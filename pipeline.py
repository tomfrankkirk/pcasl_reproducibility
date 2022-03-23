import functools
import multiprocessing
import os.path as op 
import numpy as np 
import subprocess 
import regtricks as rt 
import itertools
import glob 
from fsl import wrappers
import os 
import shutil
import nibabel as nib 
import toblerone as tob 

ROOT = '/Users/thomaskirk/Data/pcasl2/raw'
OUTROOT = '/Users/thomaskirk/Data/pcasl2/derivatives/preprocessed'
REGROOT = '/Users/thomaskirk/Data/pcasl2/derivatives/registration'
OXROOT = '/Users/thomaskirk/Data/pcasl2/derivatives/fitted'

SUBJECTS = [1,2,3,4,5,6,7,8]
SESSIONS = ['A','B','C']
RPT_SESSIONS = {
    1 : "A",
    2 : "C",
    3 : "A",
    4: "A",
    5: "C",
    6: "B",
    7: "A",
    8: "B"
}

SESS_RUNS = { sub : [ (sess, 1) for sess in SESSIONS ] for sub in SUBJECTS } 
for sub in SUBJECTS: 
    SESS_RUNS[sub] += [(RPT_SESSIONS[sub], 2)]

FLIRT = { 'dof': 6, 'finesearch': 6, 'coarsesearch': 30, 'bins': 512, 
          'noresampblur': True, 'forcescaling': True, }

get_anat_path = lambda sub: op.join(ROOT, f'../derivatives/anat/s{sub:02d}.anat')

def bids_path(**kwargs): 
    suffix = kwargs.pop('suffix', None)
    ext = kwargs.pop('ext', '.nii.gz')
    stem = f"sub-{kwargs['sub']}/"
    stem += "_".join([ f"{k}-{v}" for k,v in kwargs.items() ])
    if suffix: 
        stem += f"_{suffix}"
    return stem + ext

def reg_path(**kwargs):
    p = bids_path(**kwargs) 
    p = op.join(REGROOT, p)
    os.makedirs(op.dirname(p), exist_ok=True)
    return p 

def out_path(**kwargs):
    p = bids_path(**kwargs) 
    p = op.join(OUTROOT, p)
    os.makedirs(op.dirname(p), exist_ok=True)
    return p 

def oxasl_path(**kwargs):
    p = bids_path(**kwargs, ext='') 
    p = op.join(OXROOT, p + '/')
    os.makedirs(p, exist_ok=True)
    return p 

def do_registrations(): 

    REFRESH = False 

    for sub in SUBJECTS: 
        anat_path = get_anat_path(sub)
        struct_path = op.join(anat_path, 'T1_biascorr.nii.gz')
        wm_seg = op.join(anat_path, 'T1_fast_pve_2.nii.gz')

        for sess,run in SESS_RUNS[sub]: 
            try:
                print('Registration', sub, sess, run)

                raw_asl_path = op.join(ROOT, f's{sub:02d}/{sess}/rest{run}.nii.gz')
                raw_m0h_path = op.join(ROOT, f's{sub:02d}/{sess}/calibration_head1.nii.gz')
                raw_m0b_path = op.join(ROOT, f's{sub:02d}/{sess}/calibration_body1.nii.gz')

                for path in [raw_asl_path, raw_m0h_path, raw_m0b_path]: 
                    if not op.exists(path): 
                        raise ValueError("could not find", path)
                
                # ROI vol0 of calib and ASL 
                # FLIRT calib 2 struct 
                # FLIRT ASL 2 calib
                # MCFLIRT ASL 
                asl_spc = rt.ImageSpace(raw_m0h_path)
                m0h_vol0 = nib.load(raw_m0h_path).dataobj[:,:,:,0]
                m0b_vol0 = nib.load(raw_m0b_path).dataobj[:,:,:,0]
                asl_vol0 = nib.load(raw_asl_path).dataobj[:,:,:,0]

                # Save first volume of ASL and M0 for sanity checking and registration targets 
                m0h_vol0_path = reg_path(sub=sub, sess=sess, suffix='m0h_vol0')
                m0b_vol0_path = reg_path(sub=sub, sess=sess, suffix='m0b_vol0')
                asl_vol0_path = reg_path(sub=sub, sess=sess, suffix='asl_vol0')
                asl_spc.save_image(m0h_vol0, m0h_vol0_path)
                asl_spc.save_image(m0b_vol0, m0b_vol0_path)
                asl_spc.save_image(asl_vol0, asl_vol0_path)

                # DC stuff 
                try: 
                    allfmaps = sorted(glob.glob(op.join(ROOT, f's{sub:02d}/{sess}/other/images_*_fieldmapbrain*.nii.gz')))
                    fmapphase = allfmaps[-1]
                    fmapmag = allfmaps[0]

                    epi_dc = reg_path(sub=sub, sess=sess, suffix='epi_dc_warp')
                    if not op.exists(epi_dc) or REFRESH: 
                        fmapmagbrain_noisy = reg_path(sub=sub, sess=sess, suffix='fmapmagbrain_noisy')
                        fmapmagbrain = reg_path(sub=sub, sess=sess, suffix='fmapmagbrain')
                        fmap = reg_path(sub=sub, sess=sess, suffix='fmaprads')
                        smap = reg_path(sub=sub, sess=sess, suffix='shiftmap')
                        fmap_spc = rt.ImageSpace(fmapmag)
                        fmapmag_vol0 = nib.load(fmapmag).dataobj
                        if fmapmag_vol0.ndim > 3: 
                            fmapmag_vol0 = fmapmag_vol0[:,:,:,0]

                        # BET and erode the fieldmap magnitude image 
                        wrappers.bet(fmap_spc.make_nifti(fmapmag_vol0), fmapmagbrain_noisy, f=0.6)
                        wrappers.fslmaths(fmapmagbrain_noisy).ero().run(fmapmagbrain)

                        # Prepare a fieldmap in rad/s
                        # fsl_prepare_fieldmap SIEMENS images_010_fieldmapbrain12001.nii.gz fmapmagbrain.nii.gz fmap.nii.gz 2.46 
                        cmd = ["fsl_prepare_fieldmap", "SIEMENS", fmapphase, fmapmagbrain, fmap, "2.46"] 
                        subprocess.run(cmd, check=True)

                        # Prepare the DC warp (using shiftmap output)
                        # fugue -i calibration_body1.nii.gz --dwell=0.0004 --loadfmap=fmap.nii.gz  --unwarpdir=y- --phaseconj
                        cmd = (f"fugue -i {m0h_vol0_path} --dwell=0.0004 --loadfmap={fmap} " 
                            f"--unwarpdir=y- --phaseconj --saveshift={smap}")
                        subprocess.run(cmd, shell=True, check=True)
                
                        # Convert shiftmap into a general FSL warpfield 
                        # convertwarp -s shift.nii.gz -d y- -r calibration_head1.nii.gz -o epi_dc.nii.gz --relout
                        cmd = f"convertwarp -s {smap} -d y- -r {m0h_vol0_path} -o {epi_dc} --relout"
                        subprocess.run(cmd, shell=True, check=True)

                    # Apply DC to the body M0 to assist registration onto structural data 
                    epi_dc = rt.NonLinearRegistration.from_fnirt(epi_dc, asl_spc, asl_spc, intensity_correct=True)

                except: 
                    print("Could not generate DC warp for", sub, sess, run)
                    epi_dc = rt.Registration.identity()

                m0_vol0_dc = epi_dc.apply_to_image(asl_spc.make_nifti(m0b_vol0), asl_spc)
                m0_vol0_dc_path = reg_path(sub=sub, sess=sess, suffix='m0_vol0_dc')
                m0_vol0_dc.to_filename(m0_vol0_dc_path)

                # FLIRT between M0b and struct. M0b has better intensity uniformity than M0h. 
                m0_2_struct = reg_path(sub=sub, sess=sess, suffix='flirt_m0_2_struct', ext='.mat')
                struct_2_m0_path = reg_path(sub=sub, sess=sess, suffix='flirt_struct_in_m0')
                if not op.exists(m0_2_struct) or REFRESH:
                    rt.flirt(src=m0_vol0_dc, ref=struct_path, omat=m0_2_struct,
                            cost='bbr', wmseg=wm_seg, **FLIRT)
                m0_2_struct = rt.Registration.from_flirt(m0_2_struct, asl_spc, struct_path)
                m0_2_struct.inverse().apply_to_image(struct_path, asl_spc).to_filename(struct_2_m0_path)

                # FLIRT between M0h and ASL 
                asl_2_m0 = reg_path(sub=sub, sess=sess, run=run, suffix='flirt_asl_2_m0', ext='.mat')
                if not op.exists(asl_2_m0) or REFRESH:
                    rt.flirt(src=asl_spc.make_nifti(asl_vol0), ref=asl_spc.make_nifti(m0h_vol0), omat=asl_2_m0, **FLIRT)
                asl_2_m0 = rt.Registration.from_flirt(asl_2_m0, asl_spc, asl_spc)
                m0_2_asl_path = reg_path(sub=sub, sess=sess, suffix='flirt_m0_in_asl')
                asl_2_m0.inverse().apply_to_image(asl_spc.make_nifti(m0h_vol0), asl_spc).to_filename(m0_2_asl_path)

                # MCFLIRT ASL 
                asl_mc = reg_path(sub=sub, sess=sess, run=run, suffix='mcflirt_asl', ext='.mat')
                if not op.exists(asl_mc) or REFRESH: 
                    shutil.rmtree(asl_mc, ignore_errors=True)
                    rt.mcflirt(raw_asl_path, refvol=0, out=asl_mc[:-4], mats=True)
                asl_mc = rt.MotionCorrection.from_mcflirt(asl_mc, asl_spc, asl_spc)

                # Save the struct 2 asl matrix for native space 
                str2asl_path = reg_path(sub=sub, sess=sess, run=run, space='native', suffix='flirt_struct_2_asl', ext='.mat')
                str2asl = rt.chain(m0_2_struct.inverse(), asl_2_m0.inverse())
                str2asl.save_fsl(str2asl_path, struct_path, asl_spc)

                # And for common space 
                str2asl_path = reg_path(sub=sub, sess=sess, run=run, space='common', suffix='flirt_struct_2_asl', ext='.mat')
                rt.Registration.identity().save_fsl(str2asl_path, struct_path, asl_spc)

            except Exception as e: 
                print('Skipping registration', sub, sess, run, e)


def prepare_data(): 

    for sub in SUBJECTS: 
        anat_path = get_anat_path(sub)
        fsdir = op.join(ROOT, f'../derivatives/fs/s{sub:02d}')
        firstdir = op.join(anat_path, 'first_results')
        struct_path = op.join(anat_path, 'T1_biascorr_brain.nii.gz')
        mask_path = op.join(anat_path, 'T1_biascorr_brain_mask.nii.gz')
        struct_spc = rt.ImageSpace(struct_path)

        for sess,run in SESS_RUNS[sub]: 
                print('Data preprocessing', sub, sess, run)
                try: 

                    raw_asl_path = op.join(ROOT, f's{sub:02d}/{sess}/rest{run}.nii.gz')
                    raw_m0h_path = op.join(ROOT, f's{sub:02d}/{sess}/calibration_head1.nii.gz')
                    raw_m0b_path = op.join(ROOT, f's{sub:02d}/{sess}/calibration_body1.nii.gz')

                    asl_spc = rt.ImageSpace(raw_m0h_path)
                    common_spc = struct_spc.resize_voxels(asl_spc.vox_size / struct_spc.vox_size)

                    m0_2_struct = reg_path(sub=sub, sess=sess, suffix='flirt_m0_2_struct', ext='.mat')
                    m0h_2_struct = rt.Registration.from_flirt(m0_2_struct, asl_spc, struct_path)

                    asl_2_m0h = reg_path(sub=sub, sess=sess, run=run, suffix='flirt_asl_2_m0', ext='.mat')
                    asl_2_m0h = rt.Registration.from_flirt(asl_2_m0h, asl_spc, asl_spc)

                    asl_mc = reg_path(sub=sub, sess=sess, run=run, suffix='mcflirt_asl', ext='.mat')
                    asl_mc = rt.MotionCorrection.from_mcflirt(asl_mc, asl_spc, asl_spc)

                    try: 
                        epi_dc = reg_path(sub=sub, sess=sess, suffix='epi_dc_warp')
                        epi_dc = rt.NonLinearRegistration.from_fnirt(epi_dc, asl_spc, asl_spc, intensity_correct=True)
                    except: 
                        epi_dc = rt.Registration.identity()

                    gm_struct = op.join(anat_path, 'T1_fast_pve_1.nii.gz')
                    wm_struct = op.join(anat_path, 'T1_fast_pve_2.nii.gz') 

                except: 
                    print('Data preprocessing skip', sub, sess, run)
                    continue 

                # Create TI volume (needed for slice time correction in common space)
                tivol_native = np.indices(asl_spc.size, dtype=np.float32)[2]
                tivol_native *= 0.0452
                plds_dense = 1.4 + np.repeat([1.65,1.9,2.15,2.4,2.65,2.9], 8)
                tivol_native = tivol_native[...,None] + plds_dense[None,None,None,:]
                tivol_native[:,:,0,:] = 0 
                tivol_native[:,:,-1,:] = 0 #clip out top and bottom slices because they show edge effects after MCFLRIT 

                # Native space output
                asl2nat = rt.chain(asl_mc, asl_2_m0h, epi_dc, asl_2_m0h.inverse()) 
                m02nat = rt.chain(epi_dc, asl_2_m0h.inverse())
                m02nat_lin = asl_2_m0h.inverse()
                struct2nat = rt.chain(m0h_2_struct.inverse(), asl_2_m0h.inverse())

                # Common space output 
                asl2cmn = rt.chain(asl_mc, asl_2_m0h, epi_dc, m0h_2_struct)
                m02cmn = rt.chain(epi_dc, m0h_2_struct)
                m02cmn_lin = m0h_2_struct
                struct2cmn = rt.Registration.identity() 

                for sname, spc, (asltrans, m0trans, m0trans_lin, structtrans) in zip(
                    ['native', 'common'], [asl_spc, common_spc],
                    [(asl2nat, m02nat, m02nat_lin, struct2nat), (asl2cmn, m02cmn, m02cmn_lin, struct2cmn)]): 

                    asl_path = out_path(sub=sub, sess=sess, run=run, space=sname, acq='asl_rest')
                    asl_out = asltrans.apply_to_image(raw_asl_path, spc, order=1)
                    asl_out.to_filename(asl_path)

                    m0_path = out_path(sub=sub, sess=sess, run=run, space=sname, acq='m0h')
                    m0_out = m0trans.apply_to_image(raw_m0h_path, spc, order=1)
                    m0_out.to_filename(m0_path)

                    m0_path = out_path(sub=sub, sess=sess, run=run, space=sname, acq='m0b')
                    m0_out = m0trans.apply_to_image(raw_m0b_path, spc, order=1)
                    m0_out.to_filename(m0_path) 

                    ti_path = out_path(sub=sub, sess=sess, run=run, space=sname, suffix='tivol')
                    ti_out = m0trans_lin.apply_to_array(tivol_native, asl_spc, spc, order=1)
                    spc.save_image(ti_out, ti_path)

                    strpath = out_path(sub=sub, sess=sess, run=run, space=sname, acq='anat')
                    str_out = structtrans.apply_to_image(struct_path, spc, order=1)
                    str_out.to_filename(strpath)

                    maskpath = out_path(sub=sub, sess=sess, run=run, space=sname, suffix='mask')
                    mask_out = structtrans.apply_to_image(mask_path, spc, order=1)
                    mask_out = (mask_out.dataobj > 0.8) & (ti_out[...,0] > 0.5)
                    spc.save_image(mask_out, maskpath)

                    # FAST PV estimates 
                    suff = '' if sname == 'native' else 'naive'
                    gmpath = out_path(sub=sub, sess=sess, run=run, space=sname, pvgm='fast', suffix=suff)
                    gm_out = structtrans.apply_to_image(gm_struct, spc, order=1)
                    gm_out.to_filename(gmpath)

                    wmpath = out_path(sub=sub, sess=sess, run=run, space=sname, pvwm='fast', suffix=suff)
                    wm_out = structtrans.apply_to_image(wm_struct, spc, order=1)
                    wm_out.to_filename(wmpath)

                    # Double resampled PVs - we map from the native space back to common 
                    if sname == 'native': 
                        gmpath = out_path(sub=sub, sess=sess, run=run, space='common', pvgm='fast', suffix='double')
                        gm_out = structtrans.inverse().apply_to_image(gm_out, common_spc, order=1)
                        gm_out.to_filename(gmpath)

                        wmpath = out_path(sub=sub, sess=sess, run=run, space='common', pvwm='fast', suffix='double')
                        wm_out = structtrans.inverse().apply_to_image(wm_out, common_spc, order=1)
                        wm_out.to_filename(wmpath)

                    # Toblerone PV estimates 
                    gmpath = out_path(sub=sub, sess=sess, run=run, space=sname, pvgm='tob', suffix=suff)
                    if not op.exists(gmpath):
                        tob_pvs = tob.pvestimation.complete(spc, struct2ref=structtrans, fslanat=anat_path, 
                                                            firstdir=firstdir, fsdir=fsdir, ones=False)

                        spc.save_image(tob_pvs['GM'], gmpath)
                        wmpath = out_path(sub=sub, sess=sess, run=run, space=sname, pvwm='tob', suffix=suff)
                        spc.save_image(tob_pvs['WM'], wmpath)

                        if sname == 'native': 
                            gmpath = out_path(sub=sub, sess=sess, run=run, space='common', pvgm='tob', suffix='double')
                            gm_out = structtrans.inverse().apply_to_array(tob_pvs['GM'], spc, common_spc, order=1)
                            common_spc.save_image(gm_out, gmpath)

                            wmpath = out_path(sub=sub, sess=sess, run=run, space='common', pvwm='tob', suffix='double')
                            wm_out = structtrans.inverse().apply_to_array(tob_pvs['WM'], spc, common_spc, order=1)
                            common_spc.save_image(wm_out, wmpath)



def run_oxasl():

    # 
    oxasl_args = [ '--casl', '--fixbolus', '--te', '13', '--tr', '6', '--iaf', 
                   'tc', '--debug', '--ibf', 'rpt', '--rpts', '8', '--overwrite', 
                   '--cmethod', 'voxel', '--bolus', '1.4', '--calib-aslreg', 
                   '--no-report', '--plds', '1.65,1.9,2.15,2.4,2.65,2.9']

    jobs = [] 
    for sub in SUBJECTS: 
        for (sess,run),space,pvec in itertools.product(SESS_RUNS[sub], ['common', 'native'], ['tob', 'fast']): 
            print('Fitting', sub, sess, run, space)
        
            suff = '' if space == 'native' else 'naive'
            anat = get_anat_path(sub)
            asl_path = out_path(sub=sub, sess=sess, run=run, space=space, acq='asl_rest')
            m0_path = out_path(sub=sub, sess=sess, run=run, space=space, acq='m0h')
            m0b_path = out_path(sub=sub, sess=sess, run=run, space=space, acq='m0b')
            gmpath = out_path(sub=sub, sess=sess, run=run, space=space, pvgm=pvec, suffix=suff)
            wmpath = out_path(sub=sub, sess=sess, run=run, space=space, pvwm=pvec, suffix=suff)
            ti_path = out_path(sub=sub, sess=sess, run=run, space=space, suffix='tivol')
            maskpath = out_path(sub=sub, sess=sess, run=run, space=space, suffix='mask')
            str2asl_path = reg_path(sub=sub, sess=sess, run=run, space=space, suffix='flirt_struct_2_asl', ext='.mat')
            odir = oxasl_path(sub=sub, sess=sess, run=run, space=space, pvec=pvec, suffix=suff)

            cmd = ["oxasl", '-i', asl_path, '-c', m0_path, '--pvcorr', '-o', odir, 
                    '-m', maskpath, '--pvgm', gmpath, '--pvwm', wmpath, '--tiimg',
                    ti_path, '--cref', m0b_path, *oxasl_args]
        
            if not op.exists(op.join(odir, 'output_pvcorr')) or True:
                jobs.append(" ".join(cmd))

            if space == 'native': 
                continue 

            gmpath = out_path(sub=sub, sess=sess, run=run, space=space, pvgm=pvec, suffix='double')
            wmpath = out_path(sub=sub, sess=sess, run=run, space=space, pvwm=pvec, suffix='double')
            odir = oxasl_path(sub=sub, sess=sess, run=run, space=space, pvec=pvec, suffix='double')

            cmd = ["oxasl", '-i', asl_path, '-c', m0_path, '--pvcorr', '-o', odir, 
                    '-m', maskpath, '--pvgm', gmpath, '--pvwm', wmpath, '--tiimg', 
                    ti_path, '--cref', m0b_path, *oxasl_args]
            
            if not op.exists(op.join(odir, 'output_pvcorr')) or True:
                jobs.append(" ".join(cmd))

    worker = functools.partial(subprocess.run, check=True, shell=True)

    print(*jobs, sep='\n\n')
    # for j in jobs: print(j)
    # worker(jobs[0])
    with multiprocessing.Pool() as p: 
        p.map(worker, jobs)

    for sub in SUBJECTS: 
        anat_path = get_anat_path(sub)
        struct_path = op.join(anat_path, 'T1_biascorr_brain.nii.gz')
        struct_spc = rt.ImageSpace(struct_path)

        for (sess,run),pvec in itertools.product(SESS_RUNS[sub], ['tob','fast']): 
            try: 
                print('Post transform', sub, sess, run, pvec)    

                indir = oxasl_path(sub=sub, sess=sess, run=run, space='native', pvec=pvec)
                outdir = oxasl_path(sub=sub, sess=sess, run=run, space='common', pvec=pvec, suffix='post')

                raw_m0h_path = op.join(ROOT, f's{sub:02d}/{sess}/calibration_head1.nii.gz')
                asl_spc = rt.ImageSpace(raw_m0h_path)
                common_spc = struct_spc.resize_voxels(asl_spc.vox_size / struct_spc.vox_size)
                m0_2_struct = reg_path(sub=sub, sess=sess, suffix='flirt_m0_2_struct', ext='.mat')
                m0h_2_struct = rt.Registration.from_flirt(m0_2_struct, asl_spc, struct_path)
                asl_2_m0h = reg_path(sub=sub, sess=sess, run=run, suffix='flirt_asl_2_m0', ext='.mat')
                asl_2_m0h = rt.Registration.from_flirt(asl_2_m0h, asl_spc, asl_spc)
                asl2cmn = rt.chain(asl_2_m0h, m0h_2_struct)

                gmcbf = op.join(indir, 'output/native/perfusion_calib.nii.gz')
                gpath = op.join(outdir, 'output/native/perfusion_calib.nii.gz')
                os.makedirs(op.dirname(gpath), exist_ok=True) 
                gout = asl2cmn.apply_to_image(gmcbf, common_spc, order=1)
                gout.to_filename(gpath)
                
                gmcbf = op.join(indir, 'output_pvcorr/native/perfusion_calib.nii.gz')
                gpath = op.join(outdir, 'output_pvcorr/native/perfusion_calib.nii.gz')
                os.makedirs(op.dirname(gpath), exist_ok=True) 
                gout = asl2cmn.apply_to_image(gmcbf, common_spc, order=1)
                gout.to_filename(gpath)
                wpath = op.join(outdir, 'output_pvcorr/native/perfusion_wm_calib.nii.gz')
                wmcbf = op.join(indir, 'output_pvcorr/native/perfusion_wm_calib.nii.gz')
                wout = asl2cmn.apply_to_image(wmcbf, common_spc, order=1)
                wout.to_filename(wpath)

            except Exception as e: 
                print('Skipping post transform', sub, sess, run, e)    

if __name__ == '__main__':

    # do_registrations() 
    # prepare_data() 
    run_oxasl()
