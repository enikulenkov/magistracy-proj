#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include "f77blas.h"
#include "timer.h"
#include "utils.h"

int g_sizes[] = {1000, 2000, 3500, 5000, 6500};

#define SIZES_COUNT sizeof(g_sizes)/sizeof(g_sizes[0])

void calc_blas(double *A, double *B, double *C, int n)
{
  double alpha = 1.0;
  double beta = 0;
  char trans = 'N';
 
  dgemm_(&trans, &trans, &n, &n, &n, &alpha, A, &n, B, &n, &beta, C, &n);
}


void calc_cublas(double *A, double *B, double *C, int n)
{
  cublasStatus_t stat;
#define CUDA_CALL(_call) \
  { \
    cudaError_t cudaStat = (_call); \
    if ( cudaStat != cudaSuccess) { \
      printf("CUDA Error: %s (error %d)\n", cudaGetErrorString(cudaStat), cudaStat); \
      goto error; \
    } \
  }

#define CUBLAS_CALL(_call, _msg) \
  { \
    cublasStatus_t stat = (_call); \
    if ( stat != CUBLAS_STATUS_SUCCESS) { \
      printf("CUBLAS Error:" _msg ", (error %d)\n", stat); \
      goto error; \
    } \
  }
 
  cublasHandle_t handle;
  cudaEvent_t copyStart, copyEnd;
  float elapsedTime;

  double *dA;
  double *dB;
  double *dC;

  CUDA_CALL( cudaMalloc((void **)&dA, n*n*sizeof(double)));
  CUDA_CALL( cudaMalloc((void **)&dB, n*n*sizeof(double)));
  CUDA_CALL( cudaMalloc((void **)&dC, n*n*sizeof(double)));

  CUBLAS_CALL( cublasCreate(&handle), "Cublas init error");

  cudaEventCreate(&copyStart);
  cudaEventCreate(&copyEnd);

  cudaEventRecord(copyStart, 0);

  CUBLAS_CALL( cublasSetMatrix(n, n, sizeof(double), A, n, dA, n), "Cublas setMatrix A error");

  CUBLAS_CALL( cublasSetMatrix(n, n, sizeof(double), B, n, dB, n), "Cublas setMatrix B error");

  cudaEventRecord(copyEnd, 0);
  cudaEventSynchronize(copyEnd);

  cudaEventElapsedTime(&elapsedTime, copyStart, copyEnd);

  printf("cuBlas copy to device: %3.1f ms\n", elapsedTime);

  {
    cublasOperation_t trans = CUBLAS_OP_N;
    double alpha = 1.0;
    double beta = 0.0;
    cudaError_t err = 0;

    CUBLAS_CALL( cublasDgemm(handle, trans, trans, n, n, n, &alpha, dA, n, dB, n, &beta, dC, n), "DGEMM error");

    err = cudaDeviceSynchronize();

    cudaEventRecord(copyStart, 0);

    CUBLAS_CALL( cublasGetMatrix(n, n, sizeof(double), dC, n, C, n), "Cublas getMatrix C error");

    cudaEventRecord(copyEnd, 0);
    cudaEventSynchronize(copyEnd);

    cudaEventElapsedTime(&elapsedTime, copyStart, copyEnd);

    printf("cuBlas copy from device: %3.1f ms\n", elapsedTime);
  }

error:
  cudaEventDestroy(copyStart);
  cudaEventDestroy(copyEnd);
  cudaFree(dA);
  cudaFree(dB);
  cudaFree(dC);
}


void run_test(int n)
{
  double *A = (double *)malloc(n*n*sizeof(double));
  double *B = (double *)malloc(n*n*sizeof(double));
  double *C = (double *)malloc(n*n*sizeof(double));
  timer_ctx_t timer;

  timer_init(&timer);

  random_matr(A, n, n);
  random_matr(B, n, n);

  printf("Size: %d\n", n);

  timer_start(&timer);
  calc_blas(A, B, C, n);
  timer_stop(&timer);
  printf("Blas time: %s\n", timer_diff_as_str(&timer));

  timer_start(&timer);
  calc_cublas(A, B, C, n);
  timer_stop(&timer);
  printf("cuBlas time: %s\n", timer_diff_as_str(&timer));

  printf("-----------\n", n);

  free(A);
  free(B);
  free(C);
}


int main(void)
{
  int i;

  for (i=0; i < SIZES_COUNT; i++)
  {
    run_test(g_sizes[i]);
  }
}
