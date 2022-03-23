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
import toblerone
import glob
import oxasl.oxford_asl

ROOT = '/home/ibmeuser/LocalData/pcasl_reproducibility/raw'
OUTROOT = '/home/ibmeuser/LocalData/pcasl_reproducibility/processing'

SUBJECT_SESSIONS = {
    "01" : ["A", "B"],
    "02" : ["A", "B", "C"],
    "03" : ["A", "B", "C"],
    "04" : ["A", "B"],
    "05" : ["B", "C"],
    "06" : ["A", "B", "C"],
    "07" : ["B", "C"]
}

SUBJECT_SESSIONS2 = {
    "01" : ["A"],# "B"],
    "02" : ["C"],#"A", "B", 
    "03" : ["A"],# "B", "C"],
    "04" : ["A"],# "B"],
    "05" : ["C"],#"B", 
    "06" : ["B"],#"A", , "C"
#  #"07" : ["B", "C"]
}




IMAGES = []
for sub in SUBJECTS
images += itertools.product([('%02i'%(r), s) for r in range(1,len(SUBJECT_SESSIONS)+1) for s in SUBJECT_SESSIONS['%02i' % r]], 
    ["rest1", "task1"])
images += itertools.product([('%02i'%(r), s) for r in range(1,7) for s in SUBJECT_SESSIONS2['%02i' % r]], 
    ["rest2", "task2"])
 
images = [(sub,sess,img) for (sub,sess), img in images]

subject_directory = lambda s: op.join(ROOT, 's' + s)


def anatdir_forsub(n):
    return op.join(OUTROOT, 's%d.anat' % int(n))

def prepare_anat_dir(s):
    struct = op.join(ROOT, 's0%d' % s, 'struc.nii.gz')
    anatdir = anatdir_forsub(s)
    if not op.isdir(anatdir):
        make_surf_anat_dir(struct=struct, out=anatdir)

def sess_dir(sub,sess):
    return op.join(ROOT, 's%s' % sub, sess)

def fmap_forsess(sub,sess):
    return op.join(ROOT, 'fieldmaps', 's%s' % sub, sess, 'fieldmap_ero_prepared.nii.gz')

def fmapmagbrain_forsess(sub,sess):
    return op.join(ROOT, 'fieldmaps', 's%s' % sub, sess, 'fieldmap_mag_brainext_ero.nii.gz')

def outdir_forimg(sub,sess,img):
    return op.join(OUTROOT, sub, sess, img)

def worker(args):
    # sys.argv = []
    # sys.argv = args 
    # oxasl.oxford_asl.main()
    a = ['oxasl'] + args[1:]
    subprocess.run(a)

def shell(arg):
    subprocess.run(arg, shell=True)

SIGMA = 1.7

if __name__ == '__main__':

    # with multiprocessing.Pool(7) as p: 
    #     p.map(prepare_anat_dir, range(1,8))
    init_runs = [] 
    pv_runs = [] 
    estimate_runs = [] 

    for sub,sess,img in images:
        print(sub,sess,img)

        if img in ["task1", "task2"]:
            cmd = "./chop_chop.sh %s %s %s" % (sub, sess, img)
            # subprocess.run(cmd, shell=True)
            jobs = [ img + suff for suff in ["_task", "_rest"] ]
        else:
            jobs = [ img ]

        N_JOBS = len(jobs)
        for job in jobs[0:N_JOBS]:

            outdir = outdir_forimg(sub, sess, job)
            rpt = img[-1]

            # FAST and Tob without PSF
            anat = anatdir_forsub(sub)
            sessdir = sess_dir(sub,sess)
            asl = op.join(sessdir, '%s.nii.gz' % job)
            cmethod = 'single'
            calib = op.join(sessdir, 'calibration_head%s.nii.gz' % rpt)
            cref = op.join(sessdir, 'calibration_body%s.nii.gz' % rpt)
            if not op.isfile(calib):
                calib = op.join(sessdir, 'calibration_head1.nii.gz')
                cref = op.join(sessdir, 'calibration_body1.nii.gz')
            fmap = fmap_forsess(sub,sess)
            fmapmagbrain = fmapmagbrain_forsess(sub,sess)
            fmapmag = op.join(sessdir, 'other', 'images_009_fieldmapbrain11001.nii.gz')

            if not op.isfile(op.join(outdir, 'output_surf_pvcorr', 'native', 'perfusion_calib.nii.gz')):
                args = ['', '-i', asl, '--fslanat', anat, '--overwrite',
                '-c', calib, '--cref', cref, '--fmap', fmap, 
                '--fmapmag', fmapmag, '--fmapmagbrain', fmapmagbrain,
                '--output', outdir, '--cmethod', cmethod,
                '--tis=1.65,1.9,2.15,2.4,2.65,2.9', '--ibf=rpt', '--bolus=1.4', 
                '--slicedt=0.0452', '--fixbolus', '--cores=3',
                '--casl', '--iaf=tc', '--tr=6', '--te=13', '--debug',
                '--basil-options', '', '--mc', '--pvcorr', '--surf-pvcorr',
                '--no-report', '--pedir=-y', '--echospacing=0.0004']
                init_runs.append(args)

            # FAST, with PSF
            fast_psf_dir = outdir + '_fast_psf'
            subprocess.run(['mkdir', '-p', op.abspath(fast_psf_dir + '/..')])
            subprocess.run(['mkdir', '-p', op.abspath(fast_psf_dir + '/custom_pvs')])

            ref = op.join(outdir, 'reg', 'regfrom.nii.gz')
            struct = op.join(outdir, 'reg', 'regto.nii.gz')
            struct2asl = op.join(outdir, 'reg', 'struc2asl.mat')

            fastgm = op.join(outdir, 'structural', 'gm_pv_asl.nii.gz')
            fastwm = op.join(outdir, 'structural', 'wm_pv_asl.nii.gz')
            fastgm_psf = op.join(fast_psf_dir, 'custom_pvs', 'fast_gm_psf.nii.gz')
            fastwm_psf = op.join(fast_psf_dir, 'custom_pvs', 'fast_wm_psf.nii.gz')

            cmd = 'fslmaths %s -s %1.1f %s' % (fastgm, SIGMA, fastgm_psf)
            estimate_runs.append(cmd)
            cmd = 'fslmaths %s -s %1.1f %s' % (fastwm, SIGMA, fastwm_psf)
            estimate_runs.append(cmd)

            if not op.isfile(op.join(fast_psf_dir, 'output_pvcorr', 'native', 'perfusion.nii.gz')):
                args = ['', '-i', asl, '--fslanat', anat, '--overwrite',
                '-c', calib, '--cref', cref, '--fmap', fmap, 
                '--fmapmag', fmapmag, '--fmapmagbrain', fmapmagbrain,
                '--output', fast_psf_dir, '--cmethod', cmethod,
                '--tis=1.65,1.9,2.15,2.4,2.65,2.9', '--ibf=rpt', '--bolus=1.4', 
                '--slicedt=0.0452', '--fixbolus', '--artoff',  '--cores=3',
                '--casl', '--iaf=tc', '--tr=6', '--te=13', '--debug',
                '--basil-options', '', '--mc', '--pvgm', fastgm_psf, '--pvwm', fastwm_psf,
                '--no-report', '--pedir=-y', '--echospacing=0.0004']
                pv_runs.append(args)

            # Tob with PSF
            tob_psf_dir = outdir + '_tob_psf'
            pvdir = tob_psf_dir + '/custom_pvs'

            subprocess.run(['mkdir', '-p', op.abspath(tob_psf_dir + '/..')])
            subprocess.run(['mkdir', '-p', op.abspath(pvdir)])

            tobgm = op.join(outdir, 'surf_pvs', 'GM.nii.gz')
            tobwm = op.join(outdir, 'surf_pvs', 'WM.nii.gz')

            tobgm_psf = op.join(tob_psf_dir, 'custom_pvs', 'GM_psf.nii.gz')
            tobwm_psf = op.join(tob_psf_dir, 'custom_pvs', 'WM_psf.nii.gz')
            cmd = 'fslmaths %s -s %1.1f %s' % (tobgm, SIGMA, tobgm_psf)
            estimate_runs.append(cmd)
            cmd = 'fslmaths %s -s %1.1f %s' % (tobwm, SIGMA, tobwm_psf)
            estimate_runs.append(cmd)

            if not op.isfile(op.join(tob_psf_dir, 'output_pvcorr', 'native', 'perfusion.nii.gz')):
                args = ['', '-i', asl, '--fslanat', anat, '--overwrite',
                '-c', calib, '--cref', cref, '--fmap', fmap, 
                '--fmapmag', fmapmag, '--fmapmagbrain', fmapmagbrain,
                '--output', tob_psf_dir, '--cmethod', cmethod,
                '--tis=1.65,1.9,2.15,2.4,2.65,2.9', '--ibf=rpt', '--bolus=1.4', 
                '--slicedt=0.0452', '--fixbolus', '--artoff',  '--cores=3',
                '--casl', '--iaf=tc', '--tr=6', '--te=13', '--debug',
                '--basil-options', '', '--mc', '--pvgm', tobgm_psf, '--pvwm', tobwm_psf,
                '--no-report', '--pedir=-y', '--echospacing=0.0004']
                pv_runs.append(args)
            
    with multiprocessing.Pool(5) as p:
        p.map(worker, init_runs) 

    with multiprocessing.Pool(12) as p:
        p.map(shell, estimate_runs)

    with multiprocessing.Pool(12) as p: 
        p.map(worker, pv_runs)

    # for r in init_runs:
    #     sys.argv = r
    #     oxasl.oxford_asl.main()