#include <stdio.h>

// For the CUDA runtime routines (prefixed with "cuda_")
#include <cuda_runtime.h>
#include <cusparse.h>

/**
 * CUDA Kernel Device code
 *
 * Computes dense matrix on vector multiplication. "Naive" implementation.
 */
__global__ void
matr_vector_mult(const double *A, const double *B, double *C, int numElements)
{
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  int j;

  if (i < numElements)
  {
    double sum = 0;

    for (j = 0; j < numElements; j++)
    {
      sum += A[i * numElements + j] * B[j];
    }

    C[i] = sum;
  }
}

/* Wrapper to unify calling from Fortran and C */
extern "C" void matr_vector_mult_(double *h_A, double *h_B, double *h_C, int *numElements)
{
    // Error code to check return values for CUDA calls
    cudaError_t err = cudaSuccess;

    size_t size = *numElements * sizeof(double);

    // Allocate the device input matrix A
    double *d_A = NULL;
    err = cudaMalloc((void **)&d_A, (*numElements)*(*numElements)*sizeof(double));

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device vector A (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Allocate the device input vector B
    double *d_B = NULL;
    err = cudaMalloc((void **)&d_B, size);

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device vector B (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Allocate the device output vector C
    double *d_C = NULL;
    err = cudaMalloc((void **)&d_C, size);

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device vector C (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Copy the host input vectors A and B in host memory to the device input vectors in
    // device memory
    printf("Copy input data from the host memory to the CUDA device\n");
    err = cudaMemcpy(d_A, h_A, (*numElements)*(*numElements)*sizeof(double), cudaMemcpyHostToDevice);

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy matrix A from host to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(d_B, h_B, size, cudaMemcpyHostToDevice);

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy vector B from host to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Launch the Vector Add CUDA Kernel
    int threadsPerBlock = 256;
    int blocksPerGrid =(*numElements + threadsPerBlock - 1) / threadsPerBlock;
    printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);
    matr_vector_mult<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_B, d_C, *numElements);
    err = cudaGetLastError();

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to launch matr_vector_mult kernel (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Copy the device result vector in device memory to the host result vector
    // in host memory.
    printf("Copy output data from the CUDA device to the host memory\n");
    err = cudaMemcpy(h_C, d_C, size, cudaMemcpyDeviceToHost);

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy vector C from device to host (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Free device global memory
    err = cudaFree(d_A);

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to free device matrix A (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaFree(d_B);

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to free device vector B (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaFree(d_C);

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to free device vector C (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

extern "C" void cuda_deinit_()
{
    cudaError_t err;

    // Reset the device and exit
    err = cudaDeviceReset();

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to deinitialize the device! error=%s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

extern "C" void matr_vector_mult_sparse(int *h_A_row_indeces, int *h_A_col_indeces, double *h_A_values, int nnz,
    double *h_B, double *h_C, int *numElements)
{
  cusparseHandle_t cusparse_hndl;
  cusparseStatus_t cusparse_ret;
  cusparseMatDescr_t mat_descr;
  cudaError_t err = cudaSuccess;

  int *d_A_row_indices = NULL;
  int *d_A_col_indices = NULL;
  double *d_A_values = NULL;
  int *csrRowPtr;
  double done = 1;
  double complex_null = 0;
  cudaMalloc((void **)&d_A_row_indices, nnz*sizeof(int));
  cudaMalloc((void **)&d_A_col_indices, nnz*sizeof(int));
  cudaMalloc((void **)&d_A_values, nnz*sizeof(double));

  if (err != cudaSuccess)
  {
      fprintf(stderr, "Failed to allocate device vector A (error code %s)!\n", cudaGetErrorString(err));
      exit(EXIT_FAILURE);
  }

  // Allocate the device input vector B
  double *d_B = NULL;
  err = cudaMalloc((void **)&d_B, *numElements*sizeof(double));

  if (err != cudaSuccess)
  {
      fprintf(stderr, "Failed to allocate device vector B (error code %s)!\n", cudaGetErrorString(err));
      exit(EXIT_FAILURE);
  }

  // Allocate the device output vector C
  double *d_C = NULL;
  err = cudaMalloc((void **)&d_C, *numElements*sizeof(double));

  if (err != cudaSuccess)
  {
      fprintf(stderr, "Failed to allocate device vector C (error code %s)!\n", cudaGetErrorString(err));
      exit(EXIT_FAILURE);
  }

  // Copy the host input vectors A and B in host memory to the device input vectors in
  // device memory
  printf("Copy input data from the host memory to the CUDA device\n");
  cudaMemcpy(d_A_row_indices, h_A_row_indeces, nnz*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_A_col_indices, h_A_col_indeces, nnz*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_A_values, h_A_values, nnz*sizeof(double), cudaMemcpyHostToDevice);

  if (err != cudaSuccess)
  {
      fprintf(stderr, "Failed to copy matrix A from host to device (error code %s)!\n", cudaGetErrorString(err));
      exit(EXIT_FAILURE);
  }

  err = cudaMemcpy(d_B, h_B, *numElements*sizeof(double), cudaMemcpyHostToDevice);

  if (err != cudaSuccess)
  {
      fprintf(stderr, "Failed to copy vector B from host to device (error code %s)!\n", cudaGetErrorString(err));
      exit(EXIT_FAILURE);
  }

  cusparse_ret = cusparseCreate(&cusparse_hndl);

  if (cusparse_ret != CUSPARSE_STATUS_SUCCESS)
  {
    fprintf(stderr, "Failed to create cuSparse context (error code %d)!\n", cusparse_ret);
    exit(EXIT_FAILURE);
  }

  cusparse_ret = cusparseCreateMatDescr(&mat_descr);

  if (cusparse_ret != CUSPARSE_STATUS_SUCCESS)
  {
    fprintf(stderr, "Failed to create cuSparse matrix description (error code %d)!\n", cusparse_ret);
    exit(EXIT_FAILURE);
  }

  cusparseSetMatType(mat_descr, CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparseSetMatIndexBase(mat_descr, CUSPARSE_INDEX_BASE_ZERO);

  /* exercise conversion routines (convert matrix from COO 2 CSR format) */
  cudaMalloc((void**)&csrRowPtr,(*numElements+1)*sizeof(csrRowPtr[0]));

  cusparse_ret = cusparseXcoo2csr(cusparse_hndl, d_A_col_indices, nnz, *numElements,
                             csrRowPtr,CUSPARSE_INDEX_BASE_ZERO); 

  if (cusparse_ret != CUSPARSE_STATUS_SUCCESS)
  {
    fprintf(stderr, "Failed to make a conversion COO->CSR (error code %d)!\n", cusparse_ret);
    exit(EXIT_FAILURE);
  }

  cudaMemset(d_C, 0, nnz*sizeof(double));
  cusparse_ret = cusparseDcsrmv(cusparse_hndl,CUSPARSE_OPERATION_NON_TRANSPOSE,*numElements, *numElements,
                           done, mat_descr, d_A_values, csrRowPtr, d_A_row_indices, 
                           d_B, complex_null, d_C);

  if (cusparse_ret != CUSPARSE_STATUS_SUCCESS)
  {
    fprintf(stderr, "Failed to make cusparse mult (error code %d)!\n", cusparse_ret);
    exit(EXIT_FAILURE);
  }

  printf("Copy output data from the CUDA device to the host memory\n");
  err = cudaMemcpy(h_C, d_C, *numElements*sizeof(double), cudaMemcpyDeviceToHost);
}
