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
import oxasl.oxford_asl
import glob

ROOT = '/home/ibmeuser/LocalData/pcasl_reproducibility/raw'
OUTROOT = '/mnt/hgfs/Data/pcasl_reproducibility/processed'

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




images = []
images += itertools.product([('%02i'%(r), s) for r in range(1,len(SUBJECT_SESSIONS)+1) for s in SUBJECT_SESSIONS['%02i' % r]], 
    ["rest1", "task1"])
images += itertools.product([('%02i'%(r), s) for r in range(1,7) for s in SUBJECT_SESSIONS2['%02i' % r]], 
    ["rest2", "task2"])

NIMAGES = len(images)
images = [(sub,sess,img) for (sub,sess), img in images[0:NIMAGES] ]

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

if __name__ == '__main__':

    # with multiprocessing.Pool(7) as p: 
    #     p.map(prepare_anat_dir, range(1,8))

    for sub,sess,img in images:
        print(sub,sess,img)

        if img in ["task1", "task2"]:
            cmd = "./chop_chop.sh %s %s %s" % (sub, sess, img)
            subprocess.run(cmd, shell=True)
            jobs = [ img + suff for suff in ["_task", "_rest"] ]
        else:
            jobs = [ img ]

        for job in jobs:

            cmdargs = []
            out = outdir_forimg(sub, sess, job)

            for pvcorr in [True, False]:
                if pvcorr: 
                    outdir = out + '_pvcorr'
                else: 
                    outdir = out  
                rpt = img[-1]
                subprocess.run(['mkdir', '-p', op.abspath(outdir + '/..')])
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

                if not op.isfile(op.join(outdir, 'output', 'native', 'perfusion.nii.gz')):
                    args = ['', '-i', asl, '--fslanat', anat, '--overwrite',
                    '-c', calib, '--cref', cref, '--fmap', fmap, 
                    '--fmapmag', fmapmag, '--fmapmagbrain', fmapmagbrain,
                    '--output', outdir, '--cmethod', cmethod,
                    '--tis=1.65,1.9,2.15,2.4,2.65,2.9', '--ibf=rpt', '--bolus=1.4', 
                    '--slicedt=0.0452', '--fixbolus', '--spatial', '--artoff', 
                    '--casl', '--iaf=tc', '--tr=6', '--te=13', '--debug',
                    '--basil-options', '', '--mc',
                    '--no-report', '--pedir=-y', '--echospacing=0.0004']
                    if pvcorr:
                        args.append('--surf-pvcorr')
                        args.append('--pvcorr')
                    cmdargs.append(args)

            for args in cmdargs:
                sys.argv = []
                sys.argv = args 
                oxasl.oxford_asl.main()