#!/bin/bash

#Parameter 1 is input mat file

DEST_DIR=./res/`date +%F-%R`
MMD_RES_FILE=${DEST_DIR}/mmd_res.cml
C_RES_FILE=${DEST_DIR}/c_res.mat
VIS_FILE=${DEST_DIR}/c_vis.cml
ORIG_VIS_FILE=${DEST_DIR}/orig.cml
TIME=/usr/bin/time

mkdir -p ${DEST_DIR}
cp $1 ${DEST_DIR}
GA_OUTPUT=`$TIME -f "time %E" ./ga $1 ${C_RES_FILE}`
MMD_OUTPUT=`$TIME -f "time %E" ./sc_qdm_n $1`
cp out_f.cml ${DEST_DIR}/
./visualize_cluster ${C_RES_FILE} ${VIS_FILE}
./visualize_cluster $1 ${ORIG_VIS_FILE}

touch ${DEST_DIR}/summary.txt

echo "GA output:" >> ${DEST_DIR}/summary.txt
echo "" >> ${DEST_DIR}/summary.txt
echo "$GA_OUTPUT" >> ${DEST_DIR}/summary.txt
echo "" >> ${DEST_DIR}/summary.txt
echo "MMD output:" >> ${DEST_DIR}/summary.txt
echo "" >> ${DEST_DIR}/summary.txt
echo "$MMD_OUTPUT" >> ${DEST_DIR}/summary.txt
