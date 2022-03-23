# This code is used to call oxford_asl consistently 
# and repeatedly for a dataset which can include
# multiple subjects, each with multiple scan sessions
# and with multiple sets of scans acquired in each session
# here it is being applied to ASL task-rest data.

import os
import os.path as op
import subprocess
import multiprocessing
import nibabel
from glob import glob
import copy

NCORES = 2

# Adding the FSLDEVDIR directory to the path to access
# the development versions of FSL ASL analysis tools
os.environ["PATH"] = os.path.join(os.environ["FSLDEVDIR"], "bin") + ":" + os.environ["PATH"]


def run_oxasl_analysis(argslist):
    input_data, output_file, fph, oph = argslist

    oxasl_call = ("/home/ibmeuser/UserApps/oxford_asl/oxford_asl -i " + input_data +
        " -o " + fph + output_file + " --casl --iaf=tc " +
        "--ibf=rpt --tis=1.65,1.9,2.15,2.4,2.65,2.9 --bolus" # using while loops
        + "=1.4 --slicedt=0.0452 --fixbolus " +
        "--spatial --mc --fslanat=" + op.abspath(op.join(oph, '..')) + 
        "/struc.anat --pvcorr -c " + fph + "calibration_head1" +
        ".nii.gz --cref=" + fph + "calibration_body1" +
        ".nii.gz --tr=6 --cmethod=single --te=13 ") 
        
        
        # --fmap=" + fph + 
        # "other/fieldmap_ero_prepared.nii.gz --fmapmag=" + fph +
        # "other/images_009_fieldmapbrain11001.nii.gz --" +
        # "fmapmagbrain=" + fph + "other/fieldmap_mag_" +
        # "brainext_ero.nii.gz --echospacing=0.0004 --pedir=-y")# --debug   "--pvgm=" + pvg + " --pvwm=" + pvw +     --artoff  +
        # #"--pvgm=" + pvg + " --pvwm=" + pvw

    print(oxasl_call)
    # print()
    ### if not os.path.exists((oph + output_file)):struc.anat
    ###     os.mkdir((oph + output_file))
    os.system(oxasl_call)
    records_file = os.path.join(fph + output_file + "/command.txt")
    with open(records_file, "w+") as text_file:
        print("Oxford_asl command: {}".format(oxasl_call), file=text_file)

# specifying scan sessions for each subject - could be done
# using while loops
# N.B. there may be more data than this but I am avoiding 
# scan sessions which don't include fieldmap acquisitions
subject_sessions = {
    "01" : ["A", "B"],
    "02" : ["A", "B", "C"],
    "03" : ["A", "B", "C"],
    "04" : ["A", "B"],
    "05" : ["B", "C"],
    "06" : ["A", "B", "C"],
    "07" : ["B", "C"]
 }


subject_sessions2 = {
    "01" : ["A"],# "B"],
    "02" : ["C"],#"A", "B", 
    "03" : ["A"],# "B", "C"],
    "04" : ["A"],# "B"],
    "05" : ["C"],#"B", 
    "06" : ["B"],#"A", , "C"
#  #"07" : ["B", "C"]
}


path_root = "/home/ibmeuser/LocalData/pcasl_reproducibility/raw/"
# path_to_trans = "/mnt/hgfs/Postdoc_data_files/Andy_data/pcasl_reproducibility/raw/"
dirname = "_oxasl"

jobs = []
for i in ['%02i'%(r+1) for r in range(len(subject_sessions))]:
    print(i)
    for j in subject_sessions[i]:
        path_to_sub = (path_root + "s" + i + "/" + j + "/")
        path_to_outputs = (path_root + "s" + i + "/" + j + "/")

        filename = (path_to_sub + "rest1.nii.gz")
        out_dir = "rest1" + dirname 
        jobs.append([filename, out_dir, path_to_sub, path_to_outputs])
        # run_oxasl_analysis(filename, out_dir, path_to_sub, path_to_outputs)

        filename = (path_to_outputs + "task1_all_rest.nii.gz")
        out_dir = "task1" + dirname 
        jobs.append([filename, out_dir, path_to_sub, path_to_outputs])
        # run_oxasl_analysis(filename, out_dir, path_to_sub, path_to_outputs)

        filename = (path_to_outputs + "task1_all_task.nii.gz")
        out_dir = "task1" + dirname 
        jobs.append([filename, out_dir, path_to_sub, path_to_outputs])
        # run_oxasl_analysis(filename, out_dir, path_to_sub, path_to_outputs) 

for i in subject_sessions2:
    for j in subject_sessions2[i]:
        path_to_sub = (path_root + "s" + i + "/" + j + "/")
        path_to_outputs = ("/home/ibmeuser/pcasl_reproducibility/raw/" + "s" + i + "/" + j + "/")

        filename = (path_to_sub + "rest2.nii.gz")
        out_dir = "rest2" + dirname 
        jobs.append([filename, out_dir, path_to_sub, path_to_outputs])
        # run_oxasl_analysis(filename, out_dir, path_to_sub, path_to_outputs, mask)

        filename = (path_to_outputs + "task2_all_rest.nii.gz")
        out_dir = "task2" + dirname 
        jobs.append([filename, out_dir, path_to_sub, path_to_outputs])
        # run_oxasl_analysis(filename, out_dir, path_to_sub, path_to_outputs, mask)

        filename = (path_to_outputs + "task2_all_task.nii.gz")
        out_dir = "task2" + dirname 
        jobs.append([filename, out_dir, path_to_sub, path_to_outputs])
        # run_oxasl_analysis(filename, out_dir, path_to_sub, path_to_outputs)         

# print(jobs)
with multiprocessing.Pool(NCORES) as p:
    p.map(run_oxasl_analysis, jobs)

# # Single core mode
# for job in jobs:
#     run_oxasl_analysis(job)

