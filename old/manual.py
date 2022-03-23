import sys 
import os 
import os.path as op
import multiprocessing
import subprocess
import numpy as np
import subprocess
import nibabel 
import itertools
import scipy.io
import functools

sys.path.append('/mnt/hgfs/trin2636/Drive/DPhil/pvtools')
import toblerone
from toblerone.commandline import estimate_all_cmd

ROOT = '/home/ibmeuser/LocalData/pcasl_reproducibility/raw'
repeat = 1
SUBJECT_SESSIONS = {
    "01" : ["A", "B"],
    "02" : ["A", "B", "C"],
    "03" : ["A", "B", "C"],
    "04" : ["A", "B"],
    "05" : ["B", "C"],
    "06" : ["A", "B", "C"],
    "07" : ["B", "C"]
}

n  = {
    "01" : ["A"],# "B"],
    "02" : ["C"],#"A", "B", 
    "03" : ["A"],# "B", "C"],
    "04" : ["A"],# "B"],
    "05" : ["C"],#"B", 
    "06" : ["B"],#"A", , "C"
#  #"07" : ["B", "C"]
}


subject_directory = lambda s: op.join(ROOT, 's' + s)
reference_image = lambda sub, sess, img: op.join(subject_directory(sub), sess, 
    '%s_reference/native_space/perfusion.nii.gz' % img)

images = []
images += itertools.product([('%02i'%(r), s) for r in range(1,len(SUBJECT_SESSIONS)+1) for s in SUBJECT_SESSIONS['%02i' % r]], 
    ["rest1", "task1"])
# images += itertools.product([('%02i'%(r), s) for r in range(1,7) for s in SUBJECT_SESSIONS2['%02i' % r]], 
#     ["rest2", "task2"]) ERROR: update the repeat variable (1) for all calibration image paths below. 
NIMAGES = 5 #(images)
images = [(sub,sess,img) for (sub,sess), img in images[0:NIMAGES] ]

def motion_correct_worker(subject,session,img):

    subdir = subject_directory(subject)

    ref = reference_image(subject, session, img)
    func = op.join(subdir, session, img + '.nii.gz')

    cmd = 'mcflirt -in %s -r %s' % (func, ref)
    # shell(cmd)



def motion_correct():

    # with multiprocessing.Pool(7) as p:
    #     p.starmap(motion_correct_worker, images)

    for subject, session, img in images:
        subdir = subject_directory(subject)
        ref = reference_image(subject, session, img)
        func = op.join(subdir, session, img + '.nii.gz')
        print("Presenting", subject, session, img)
        mcfunc = op.join(subdir, session, img + '_mcf.nii.gz')
        cmd = 'fsleyes %s %s -cm blue-lightblue -a 30' % (mcfunc, ref)
        shell(cmd)

def register_all(method):
        
    register_struct_ref('normcorr', '01', 'A', 'rest1')
    tomap = functools.partial(register_struct_ref, method)
    with multiprocessing.Pool(7) as p: 
        p.starmap(tomap, images)

    for subject, session, img in images:

        # View structural, calibration registered onto structural, calibration, and motion corrected functional. 
        subdir = subject_directory(subject)
        ref2structroot = op.join(subdir, session, '%s_2_struct_%s' % (img, method))
        ref_reg_struct = ref2structroot + '.nii.gz'
        struct = op.join(subdir, 'struc.nii.gz')
        print("Presenting", subject, session, method)
        cmd = 'fsleyes %s %s' % (struct, ref_reg_struct)
        # subprocess.run(cmd, shell=True)



def register_struct_ref(method, subject, session, image):

    subdir = subject_directory(subject)

    # Paths to the ref, structural
    struct = op.join(subdir, 'struc.nii.gz')
    struct_brain = op.join(subdir, 'struc_brain.nii.gz')
    ref = reference_image(subject, session, image)

    # Paths for the calibration -> structural transform, and then the 
    ref2structroot = op.join(subdir, session, '%s_2_struct_%s' % (image,method))
    ref2struct = ref2structroot + '.mat'
    struct2ref = op.join(subdir, session, 'struct_2_%s_%s.mat' % (image,method))
    struct2ref_world = op.join(subdir, session, 'struct_2_%s_%s_world.mat' % (image,method))
    ref_reg_struct = ref2structroot + '.nii.gz'
    wm = op.join(subdir, 'struc_brain_pve_2.nii.gz')

    # EPI reg from calib -> struct
    if method == 'bbr':
        cmd = 'epi_reg --epi=%s --t1=%s --t1brain=%s --out=%s --wmseg=%s' % (ref, struct, struct_brain, ref2structroot, wm)
    elif method == 'normcorr':
        cmd = 'flirt -cost %s -dof 6 -in %s -ref %s -omat %s -out %s' % (method, ref, struct, ref2struct, ref_reg_struct)
    shell(cmd)

    # Inverse FLIRT struct -> calibration
    cmd = 'convert_xfm -omat %s -inverse %s' % (struct2ref, ref2struct)
    shell(cmd)

    # Invert and transform into world-world coords for structural -> calibration
    cmd = 'wb_command -convert-affine -from-flirt %s %s %s -to-world %s -inverse' % (ref2struct, ref, struct, struct2ref_world)
    shell(cmd)


def prepare_fast_pvs(method): 
    for subject, session, img in images:
        
        subdir = subject_directory(subject)
        ref = reference_image(subject, session, img)
        struct2ref = op.join(subdir, session, 'struct_2_%s_%s.mat' % (img, method))

        # Warp FAST onto the ref image. 
        fastcsf = op.join(subdir, 'struc_brain_pve_0.nii.gz')
        fastgm = op.join(subdir, 'struc_brain_pve_1.nii.gz')
        fastwm = op.join(subdir, 'struc_brain_pve_2.nii.gz')
        wm = op.join(subdir, session, '%s_fast_%s_wm.nii.gz' % (img, method))
        gm = op.join(subdir, session, '%s_fast_%s_gm.nii.gz' % (img, method))
        csfthresh = op.join(subdir, session, '%s_fast_%s_csf_thr.nii.gz' % (img, method))
        csf = op.join(subdir, session, '%s_fast_%s_csf.nii.gz' % (img, method))
        
        for (fast, out) in zip([fastgm, fastwm, fastcsf], [gm, wm, csf]):
            cmd = 'applywarp -i %s -r %s -o %s --premat=%s' % (fast, ref, out, struct2ref)
            subprocess.run(cmd, shell=True)

        refmask = op.join(subdir, session, '%s_mask.nii.gz' % img)
        refmasktmp = op.join(subdir, session, '%s_mask_tmp.nii.gz' % img)
        structbet = op.join(subdir, 'struc_brain.nii.gz')
        cmd = 'applywarp -i %s -r %s -o %s --premat=%s' % (structbet, ref, refmasktmp, struct2ref)
        shell(cmd)

        # Rebinarise in this space
        cmd = 'fslmaths %s -bin %s' % (refmasktmp, refmask)
        shell(cmd)       
        shell('rm %s' % refmasktmp)

        refmaskero = op.join(subdir, session, '%s_mask_ero.nii.gz' % img)
        cmd = 'fslmaths %s -roi 0 -1 0 -1 1 22 0 -1 -ero -ero -ero -ero %s' % (refmask, refmaskero)
        shell(cmd)
        cmd = 'fslmaths %s -mas %s -thr 0.9 -bin %s' % (csf, refmaskero, csfthresh)
        shell(cmd)

        print("Presenting", subject, session, img)
        cmd = 'fsleyes %s %s -cm red-yellow' % (ref, gm)
        # subprocess.run(cmd, shell=True)



def caller(string):
    subprocess.run(string, shell=True)


def estimate_surf_pvs(method): 
    for subject, session, img in images:
        
        subdir = subject_directory(subject)
        ref = reference_image(subject, session, img)
        struct2ref_world = op.join(subdir, session, 'struct_2_%s_%s_world.mat' % (img,method))

        pvdir = ROOT + '/pvtools_s%s' % subject 
        outname = '%s_%s_%s' % (session, img, method)
        estimate_all_cmd(['-ref', ref, '-pvdir', pvdir, '-struct2ref', 
            struct2ref_world, '-savesurfs', '-out', outname, '-cores', '7'])


def run_all_oxasl(method):
    commands = []

    for subject, session, img in images:
        print("Registering calibrations for", subject, session, img)

        subdir = subject_directory(subject)
        ref = reference_image(subject, session, img)

        struct = op.join(subdir, 'struc.nii.gz')
        brain = op.join(subdir, 'struc_brain.nii.gz')

        # Paths for the calibration -> structural transform, and then the 
        ref2structroot = op.join(subdir, session, '%s_2_struct_%s' % (img,method))
        ref2struct = ref2structroot + '.mat'
        struct2ref = op.join(subdir, session, 'struct_2_%s_%s.mat' % (img,method))
        struct2ref_world = op.join(subdir, session, 'struct_2_%s_%s_world.mat' % (img,method))
        ref_reg_struct = ref2structroot + '.nii.gz'
        wm = op.join(subdir, 'struc_brain_pve_2.nii.gz')
        mcfunc = op.join(subdir, session, img + '_mcf.nii.gz')

        # Register the calibration onto the functional. 
        cal = op.join(subdir, session, 'calib1.nii.gz')
        calib = op.join(subdir, session, 'calib1_%s.nii.gz' % img)
        cmd = 'flirt -in %s -ref %s -o %s -dof 6' % (cal, ref, calib)
        # shell(cmd)

        cr = op.join(subdir, session, 'calibration_body1.nii.gz')
        crroi = op.join(subdir, session, 'calibbody1.nii.gz')
        cmd = 'fslroi %s %s 0 -1 0 -1 0 -1 1 1' % (cr, crroi)
        shell(cmd)
        cref = op.join(subdir, session, 'calibbody1_%s.nii.gz' % img)
        cmd = 'flirt -in %s -ref %s -o %s -dof 6' % (crroi, ref, cref)
        # shell(cmd)

        # # Prepare a brain mask from the BETed structual
        structbet = op.join(subdir, 'struc_brain.nii.gz')
        # structmask = op.join(subdir, 'struc_mask.nii.gz')
        # cmd = 'fslmaths %s -bin %s' % (structbet, structmask)
        # subprocess.run(cmd, shell=True)

        # Resample mask to calibration space
        refmask = op.join(subdir, session, '%s_mask.nii.gz' % img)

        fwm = op.join(subdir, session, '%s_fast_%s_wm.nii.gz' % (img, method))
        fgm = op.join(subdir, session, '%s_fast_%s_gm.nii.gz' % (img, method))
        fcsfthresh = op.join(subdir, session, '%s_fast_%s_csf_thr.nii.gz' % (img, method))

        # Load surface estimates (mask them all first)
        sgm, swm, scsf = list(map(lambda t: op.join(ROOT, 
            'pvtools_s%s/%s_%s_%s_%s.nii.gz' % (subject, session, img, method, t)), ['GM', 'WM', 'nonbrain']))
        sgmmask, swmmask = list(map(lambda t: op.join(ROOT, 
            'pvtools_s%s/%s_%s_%s_%s_mask.nii.gz' % (subject, session, img, method, t)), ['GM', 'WM']))
        for unmask, masked in zip([sgm, swm], [sgmmask, swmmask]):
            cmd = 'fslmaths %s -mas %s %s' % (unmask, refmask, masked)
            subprocess.run(cmd, shell=True)

        for path in [fgm, sgmmask, fwm, swmmask]:
            assert op.isfile(path), '%s does not exist' % path

        # We use the FAST-derived CSF mask for both (the surface CSF mask includes cerebellum)

        for gm,wm,source in zip([fgm,sgmmask], [fwm,swmmask], ['fast', 'surf']):

            # And prepare the oxford asl call
            # --cref for sensitivity correction use calibration_body1
            # Distortion correction requires structural image to be used
            # So we will also manually specify the registrations we have derived earlier - 
            # we do not want oxford_asl to perform its own registration 
            fmap = op.join(subdir, session, 'other', 'fieldmap_ero_prepared.nii.gz')
            fmapmag = op.join(subdir, session, 'other', 'images_009_fieldmapbrain11001.nii.gz')
            fmapmagbrain = op.join(subdir, session, 'other', 'fieldmap_mag_brainext_ero.nii.gz')
            cmd = ('/home/ibmeuser/UserApps/oxford_asl/oxford_asl '
                + '-i %s -c %s -m %s ' % (mcfunc, calib, refmask)
                + '--ibf=rpt --tis=1.65,1.9,2.15,2.4,2.65,2.9 --bolus=1.4 '
                + '--slicedt=0.0452 --fixbolus --spatial --te=13 --artoff '
                + '--casl --iaf=tc --tr=6 --te=13 '
                + '--cref=%s ' % (cref) #Calibration, sens correction 
                + '--pvcorr --pvgm=%s --pvwm=%s ' % (gm, wm)
                + '--fmap=%s --fmapmag=%s --fmapmagbrain=%s ' % (fmap, fmapmag, fmapmagbrain)
                + '--echospacing=0.0004 --pedir=-y '
                + '--asl2struc=%s -s=%s --sbrain=%s ' % (ref2struct, struct, brain))

            outdir = op.join(subdir, session, '%s_%s_%s_voxel' % (img,source,method))
            suf = '-o %s --cmethod=voxel' % outdir
            commands.append(cmd + suf)

            outdir = op.join(subdir, session, '%s_%s_%s_reference' % (img,source,method))
            suf = '-o %s --cmethod=single --csf=%s' % (outdir, fcsfthresh)
            commands.append(cmd + suf)

    with multiprocessing.Pool(7) as p:
        p.map(shell, commands)

def shell(cmd): 
    subprocess.run(cmd, shell=True)

def compare_registrations():
    for subject, session, img in images:

        print("Presenting", subject, session, img)
        subdir = subject_directory(subject)
        ref = reference_image(subject, session, img)
        gm1 = op.join(subdir, session, '%s_fast_bbr_gm.nii.gz' % img)
        # gm2 = op.join(subdir, session, '%s_fast_normcorr_gm.nii.gz' % img)
        cmd = 'fsleyes %s %s -cm blue-lightblue' % (ref, gm1)
        shell(cmd)


def generate_references():

    commands = []

    for subject, session, img in images:
        subdir = subject_directory(subject)

        # Extract a calibration frame and save to calibroi 
        calib = op.join(subdir, session, 'calibration_head1.nii.gz')
        calibroi = op.join(subdir, session, 'calib1.nii.gz')
        cmd = 'fslroi %s %s 0 -1 0 -1 0 -1 1 1' % (calib, calibroi)
        shell(cmd)

        # And get a mask in this space. 
        brain = op.join(subdir, 'struc_brain.nii.gz')

        # And prepare the oxford asl call
        # --cref for sensitivity correction use calibration_body1
        # Distortion correction requires structural image to be used
        # So we will also manually specify the registrations we have derived earlier - 
        # we do not want oxford_asl to perform its own registration 
        func = op.join(subdir, session, img + '.nii.gz')
        cref = op.join(subdir, session, 'calibration_body%d.nii.gz' % repeat)
        fmap = op.join(subdir, session, 'other', 'fieldmap_ero_prepared.nii.gz')
        fmapmag = op.join(subdir, session, 'other', 'images_009_fieldmapbrain11001.nii.gz')
        fmapmagbrain = op.join(subdir, session, 'other', 'fieldmap_mag_brainext_ero.nii.gz')
        outdir = op.join(subdir, session, '%s_reference' % (img))
        cmd = ('/home/ibmeuser/UserApps/oxford_asl/oxford_asl '
            + '-i %s -o %s -c %s --mc ' % (func, outdir, calibroi)
            + '--ibf=rpt --tis=1.65,1.9,2.15,2.4,2.65,2.9 --bolus=1.4 '
            + '--slicedt=0.0452 --fixbolus --spatial --te=13 '
            + '--casl --iaf=tc --tr=6 --te=13 --cmethod=voxel '
            + '--artoff')

        commands.append(cmd)

    with multiprocessing.Pool(7) as p:
        p.map(shell, commands)

if __name__ == '__main__':

    for subject, session, img in images:
        path = op.join(subject_directory(subject), session, '%s.nii.gz' % img)
        assert op.isfile(path), '%s does not exist' % path

    # print("Generating references")
    # generate_references()

    method = 'bbr'

    # print("Structural - reference registration")
    # register_all(method)

    # print("Preparing FAST PVs")
    # prepare_fast_pvs(method)

    # compare_registrations()

    # print("Motion correction")
    # motion_correct()

    # print("Preparing surface PVs")
    # estimate_surf_pvs(method)

    run_all_oxasl(method)
