#!/bin/bash

ORIG=$(pwd)
SUB=$1
SESS=$2
SCAN=$3

scandir=/mnt/hgfs/Data/pcasl_reproducibility/raw/s$SUB/$SESS
echo "Chopping ${scandir}/${SCAN}.nii.gz"

cd $scandir

outdir=$scandir/${SCAN}_chopped
mkdir -p $outdir

fslroi ${SCAN}.nii.gz ${outdir}/${SCAN}_part_1b.nii.gz 2 10
fslroi ${SCAN}.nii.gz ${outdir}/${SCAN}_part_1a.nii.gz 12 2

fslroi ${SCAN}.nii.gz ${outdir}/${SCAN}_part_2b.nii.gz 32 4
fslroi ${SCAN}.nii.gz ${outdir}/${SCAN}_part_2a.nii.gz 36 8

fslroi ${SCAN}.nii.gz ${outdir}/${SCAN}_part_3b.nii.gz 62 10
fslroi ${SCAN}.nii.gz ${outdir}/${SCAN}_part_3a.nii.gz 72 2

cd $outdir

fslmerge -t ../${SCAN}_rest.nii.gz ${SCAN}_part_1a.nii.gz ${SCAN}_part_1b.nii.gz \
 ${SCAN}_part_2a.nii.gz ${SCAN}_part_2b.nii.gz ${SCAN}_part_3a.nii.gz ${SCAN}_part_3b.nii.gz

cd $scandir

fslroi ${SCAN}.nii.gz ${outdir}/${SCAN}_part_1b 18 6
fslroi ${SCAN}.nii.gz ${outdir}/${SCAN}_part_1a 24 6

fslroi ${SCAN}.nii.gz ${outdir}/${SCAN}_part_2 48 12

fslroi ${SCAN}.nii.gz ${outdir}/${SCAN}_part_3b 78 6
fslroi ${SCAN}.nii.gz ${outdir}/${SCAN}_part_3a 84 6

cd $outdir

fslmerge -t ../${SCAN}_task.nii.gz ${SCAN}_part_1a.nii.gz ${SCAN}_part_1b.nii.gz \
 ${SCAN}_part_2.nii.gz ${SCAN}_part_3a.nii.gz ${SCAN}_part_3b.nii.gz

cd $ORIG
