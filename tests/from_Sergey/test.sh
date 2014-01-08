SIZES="10 100 1000 2000"
OUTDIR=results

rm -rf $OUTDIR

mkdir -p $OUTDIR

for i in $SIZES ; do 
  echo "Size: $i"
  echo "CBLAS..."
  C_BLAS/cblas_main $i > $OUTDIR/cblas_$i.txt
  echo "CUBLAS..."
  CUBLAS/cublas_main $i > $OUTDIR/cublas_$i.txt
  echo "----------"
done
