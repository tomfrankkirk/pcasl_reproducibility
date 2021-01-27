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
import oxasl.oxford_asl
import tempfile
import shutil 
import toblerone.utils as tutils
import toblerone
import glob

RAW = '/mnt/hgfs/Data/pcasl2/raw'
OUTROOT = '/mnt/hgfs/Data/pcasl2/processed'

SUBJECTS = list(range(1,9))
SESSIONS = ['A', 'B', 'C']

RAW_SESSIONS = { s: SESSIONS for s in SUBJECTS }

RAW_SESSIONS2 = {
    1 : ["A"],
    2 : ["C"],
    3 : ["A"],
    4 : ["A"],
    5 : ["C"],
    6 : ["B"],
    7 : ["A"],
    8 : ["B"]
}

IMAGES = ['rest', 'task']

ALL = []
for (sub,sess) in [(r, s) for r in SUBJECTS for s in RAW_SESSIONS[r]]:
    for img in IMAGES:
        ALL.append((1,sub,sess,img))

for (sub,sess) in [(r, s) for r in SUBJECTS for s in RAW_SESSIONS2[r]]:
    for img in IMAGES:
        ALL.append((2,sub,sess,img))

def shell(cmd):
    subprocess.run(cmd, shell=True)

def generate_path_dicts(path_func): 
    RESTS = {}
    TASKS = {} 
    ALL = {}
    for sub in range(1,8):
        RESTS[sub] = {}
        TASKS[sub] = {}
        ALL[sub] = [] 

        for sess in RAW_SESSIONS[sub]:
            r = path_func(sub,sess,'rest1')
            RESTS[sub][sess] = [r] 
            t = path_func(sub,sess,'task1')
            TASKS[sub][sess] = [t]
            # ALL[sub].append(r,t)

        if sub in RAW_SESSIONS2: 
            for sess in RAW_SESSIONS2.get(sub):
                r = path_func(sub,sess,'rest2')
                RESTS[sub][sess].append(r) 
                t = path_func(sub,sess,'task2')
                TASKS[sub][sess].append(t)
                # ALL[sub].append(r,t)


    return RESTS, TASKS

def raw_path(rpt, sub, sess, img):
    if img.count('_'):
        img = img.replace('_', '%d_' % rpt)
        return op.join(RAW, 's%02d' % sub, sess, '%s.nii.gz' % (img))
    else:
        return op.join(RAW, 's%02d' % sub, sess, '%s%d.nii.gz' % (img, rpt))

def processed_path(rpt, sub, sess, img): 
    p = op.join(OUTROOT, '%d' % sub, '%d%s_%s.nii.gz' % (rpt, sess, img))
    d = op.join(OUTROOT, str(sub))
    if not op.isdir(d): os.mkdir(d)
    return p 

def oxasl_dir(rpt, sub, sess, img, common):
    if common: 
        return op.join(OUTROOT, str(sub), '%d%s_%s_oxasl_common' % (rpt, sess, img))
    else: 
        return op.join(OUTROOT, str(sub), '%d%s_%s_oxasl' % (rpt, sess, img))

def raw_calib_name(rpt, sub, sess):
    return op.join(RAW, 's%02d' % sub, sess, 'calibration_head%d.nii.gz' % (rpt))

def processed_calib_name(rpt, sub, sess):
    return op.join(OUTROOT, str(sub), '%d%s_calib_head.nii.gz' % (rpt,sess))

def raw_cref_name(rpt, sub, sess):
    return op.join(RAW, 's%02d' % sub, sess, 'calibration_body%d.nii.gz' % rpt)

def processed_cref_name(rpt, sub, sess):
    return op.join(OUTROOT, str(sub), '%d%s_calib_body.nii.gz' % (rpt,sess))

def anatdir_name(n):
    return op.join('/home/ibmeuser/LocalData/reproducibility_anat', '%d.anat' % n)

def reg_name(sub, ref_sess, dest_sess):
    return op.join(OUTROOT, str(sub), '1%s_hc_2_1%s_hc.mat' % (ref_sess,dest_sess))

def mcflirt_extref(series, ref, outpath):
    if sys.platform.startswith('win'):
        raise RuntimeError("MCFLIRT requires UNIX system")

    try: 
        d = tempfile.mkdtemp()
        new_series = op.join(d, op.split(series)[1])
        shutil.copy(series, new_series)
        indata = nibabel.load(new_series).get_fdata()

        target = int(indata.shape[-1] / 2)
        target_frame = op.join(d, 'target.nii.gz')
        cmd = 'fslroi %s %s %d 1' % (new_series, target_frame, target)
        shell(cmd)

        premat = op.join(d, 'target2ref.mat')
        flirt_out = op.join(d, 'target2ref.nii.gz')
        cmd = 'flirt -in %s -ref %s -omat %s -out %s' % (target_frame, ref, premat, flirt_out)
        shell(cmd)
        premat = np.loadtxt(premat)
        
        mcf_dest = op.join(d, 'mcf.nii.gz')
        cmd = 'mcflirt -in %s -refvol %d -mats -out %s' % (new_series, target, mcf_dest)
        shell(cmd)
        mats = [ np.loadtxt(m) for m in glob.glob(op.join(d, 'mcf.nii.gz.mat/*')) ]

        outdata = np.zeros(indata.shape, dtype=np.float32)
        srcspace = toblerone.ImageSpace(target_frame)
        destspace = toblerone.ImageSpace(ref)
        for fn,mat in enumerate(mats): 
            frame = op.join(d, 'frame%d.nii.gz' % fn)
            frameout = op.join(d, 'frame%d_out.nii.gz' % fn)
            cmd = 'fslroi %s %s %d 1' % (new_series, frame, fn)
            shell(cmd)
            overall = op.join(d, 'mat%d.txt' % fn)
            np.savetxt(overall, premat @ mat) 
            cmd = 'flirt -in %s -ref %s -applyxfm -init %s -out %s' % (
                frame, ref, overall, frameout)
            shell(cmd)
            outdata[:,:,:,fn] = nibabel.load(frameout).get_fdata()

        destspace.save_image(outdata, outpath)

    except Exception as e: 
        print("Error in MCFLIRT with premat")
        raise e
    finally:
        shutil.rmtree(d)

def oxasl_main(argstring):
    cmds = argstring.split()
    sys.argv[0] = ''
    sys.argv += cmds[1:]
    oxasl.oxford_asl.main()

if __name__ == "__main__":
    
    anat_cmds = []
    for sub in SUBJECTS:
        d = anatdir_name(sub)
        if not op.isdir(d):
            struct = op.join(RAW, 's%02d' % sub, 'struc.nii.gz')
            cmd = 'toblerone -fsl_fs_anat -struct %s -out %s' % (struct,d)
            anat_cmds.append(cmd)

    # with multiprocessing.Pool() as p: 
    #     p.map(shell, anat_cmds)

    for common in [False]:

        preproc_cmds = [] 
        if common: 
            for sub in SUBJECTS:

                cmd = ''
                ref_sess = RAW_SESSIONS[sub][0]
                first_calib = raw_calib_name(1, sub, ref_sess)
                first_cref = raw_cref_name(1, sub, ref_sess)
                if op.exists(first_calib):
                    out = processed_calib_name(1, sub, ref_sess)
                    outcref = processed_cref_name(1, sub, ref_sess)
                    cmd += 'fslmaths %s -Tmean %s;' % (first_calib, out)
                    cmd += 'fslmaths %s -Tmean %s;' % (first_cref, outcref)
                    ref_calib = out 

                if sub in RAW_SESSIONS: 
                    for sess in RAW_SESSIONS[sub]:
                        calib = raw_calib_name(1, sub, sess)
                        cref = raw_cref_name(1, sub, sess)
                        if calib != first_calib:
                            calibm = calib.replace('.nii.gz', '_mean.nii.gz')
                            crefm = cref.replace('.nii.gz', '_mean.nii.gz')
                            cmd += 'fslmaths %s -Tmean %s;' % (calib, calibm)
                            cmd += 'fslmaths %s -Tmean %s;' % (cref, crefm)
                            outcalib = processed_calib_name(1, sub, sess)
                            outcref = processed_cref_name(1, sub, sess)
                            reg_mat = reg_name(sub, ref_sess, sess)
                            cmd += 'flirt -in %s -ref %s ' % (calibm, ref_calib)
                            cmd += '-out %s -dof 6 -omat %s;' % (outcalib, reg_mat) 
                            cmd += 'flirt -in %s -ref %s -applyxfm -init %s -out %s;' % (
                                crefm, ref_calib, reg_mat, outcref)

                preproc_cmds.append(str(cmd))

            # [ shell(cmd) for cmd in preproc_cmds ]
            # with multiprocessing.Pool() as p:
            #     p.map(shell, preproc_cmds)

        mc_cmds = []
        oxasl_jobs = [] 
        for rpt,sub,sess,img in ALL: 

            # # chop up the tasks if needed 
            # if img.count('task') and not img.count('_'):
            #     if not op.exists(raw_path(rpt,sub,sess,img+'_rest')):
            #         cmd = './chop_chop.sh %02d %s %s%d' % (sub,sess,img,rpt)
            #         shell(cmd)

            # Motion correct timeseries to first calibration 
            if common: 
                out = processed_path(rpt,sub,sess,img)
                if not op.exists(out):
                    first_calib = processed_calib_name(1, sub, RAW_SESSIONS[sub][0])
                    asl = raw_path(rpt,sub,sess,img)
                    mc_cmds.append((asl,first_calib,out))

            # Oxasl job using the pre-processed outputs
            if common: 
                asl = processed_path(rpt,sub,sess,img)
                calib = processed_calib_name(1,sub,sess)
                cref = processed_cref_name(1,sub,sess)
            else: 
                asl = raw_path(rpt, sub, sess, img)
                calib = raw_calib_name(1, sub, sess)
                cref = raw_calib_name(1, sub, sess)

            if not op.exists(asl):
                print(asl, 'does not exist')
                continue
            
            anat = anatdir_name(sub)
            od = oxasl_dir(rpt, sub, sess, img, common)
            if (not op.exists(op.join(
                od, 'output_surf_pvcorr/native/perfusion_calib.nii.gz'))): 
                cmd = 'oxasl -i %s --fslanat %s --overwrite ' % (asl, anat)
                cmd += '-c %s --cref %s --output-mni ' % (calib, cref)
                cmd += '--tis=1.65,1.9,2.15,2.4,2.65,2.9 --ibf=rpt --bolus=1.4 '
                cmd += '--slicedt=0.0452 --fixbolus --cmethod=voxel --vars '
                cmd += '--casl --iaf=tc --tr=6 --te=13 --debug --mc '
                cmd += '--pedir=-y --echospacing=0.0004 --no-report '
                cmd += '--pvcorr --surf-pvcorr --cores=1 --output %s ' % (od)

                oxasl_jobs.append(str(cmd))

        # with multiprocessing.Pool() as p: 
        #     p.starmap(mcflirt_extref, mc_cmds)
        
        with multiprocessing.Pool(14) as p: 
           p.map(shell, oxasl_jobs)

        for sub in SUBJECTS:
            od = op.join(OUTROOT, str(sub), 'stuct_pvs')
            if not op.isdir(od):
                os.makedirs(od, exist_ok=True)
                anat = anatdir_name(sub)
                struct = op.join(anat, 'T1_biascorr_brain.nii.gz')
                cmd = 'toblerone -estimate_all -ref %s -struct2ref I ' % (struct)
                cmd += '-anat %s -cores 14 -out %s ' % (anat, od)
                shell(cmd)
