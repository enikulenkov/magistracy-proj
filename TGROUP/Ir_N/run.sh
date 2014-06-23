#!/bin/bash

#Parameter 1 is input mat file

DEST_DIR=./res/`date +%F-%R`
MMD_RES_FILE=${DEST_DIR}/mmd_res.cml
C_RES_FILE=${DEST_DIR}/c_res.mat
VIS_FILE=${DEST_DIR}/c_vis.cml
ORIG_VIS_FILE=${DEST_DIR}/orig.cml
GA_INI_FILE=./gen_alg_prefs.ini
TIME=/usr/bin/time
RUN_MMD=1

mkdir -p ${DEST_DIR}
cp $1 ${DEST_DIR}
echo "Running GA..."
GA_OUTPUT=`$TIME -f "time %E" 2>&1 ./ga $1 ${C_RES_FILE}`
cp ${GA_INI_FILE} ${DEST_DIR}/
if [ $RUN_MMD -gt 0 ] ; then
  echo "Running MMD..."
  MMD_OUTPUT=`$TIME -f "time %E" 2>&1 ./sc_qdm_n $1`
  cp out_f.cml ${DEST_DIR}/
fi
echo "Visualizing results"
./visualize_cluster ${C_RES_FILE} ${VIS_FILE} > /dev/null
./visualize_cluster $1 ${ORIG_VIS_FILE} > /dev/null

touch ${DEST_DIR}/summary.txt

echo "GA output:" >> ${DEST_DIR}/summary.txt
echo "" >> ${DEST_DIR}/summary.txt
echo "$GA_OUTPUT" >> ${DEST_DIR}/summary.txt
echo "" >> ${DEST_DIR}/summary.txt
if [ $RUN_MMD -gt 0 ] ; then
  echo "MMD output:" >> ${DEST_DIR}/summary.txt
  echo "" >> ${DEST_DIR}/summary.txt
  echo "$MMD_OUTPUT" >> ${DEST_DIR}/summary.txt
fi
