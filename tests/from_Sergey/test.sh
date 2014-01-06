SIZES="10 100 1000 2000"

for i in $SIZES ; do 
  echo "Size: $i"
  echo "--CBLAS--"
  C_BLAS/cblas_main $i
  echo "--CUBLAS--"
  CUBLAS/cublas_main $i
  echo ""
  echo "----------"
done
