#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include "f77blas.h"
#include "timer.h"
#include "utils.h"

int g_sizes[] = {1000, 2000, 3500, 5000};

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
  cudaError_t cudaStat;
  cublasStatus_t stat;
  cublasHandle_t handle;

  double *dA;
  double *dB;
  double *dC;

  cudaStat = cudaMalloc((void **)dA, n*n*sizeof(*dA));
  cudaStat = cudaMalloc((void **)dB, n*n*sizeof(*dA));
  cudaStat = cudaMalloc((void **)dC, n*n*sizeof(*dA));

  stat = cublasCreate(&handle);

  if (stat != CUBLAS_STATUS_SUCCESS)
  {
    printf("Cublas init error!!!\n");
    goto error;
  }

  stat = cublasSetMatrix(n, n, sizeof(double), A, n, dA, n);
  if (stat != CUBLAS_STATUS_SUCCESS)
  {
    printf("Cublas setMatrix error!!!\n");
    goto error;
  }

  stat = cublasSetMatrix(n, n, sizeof(double), B, n, dB, n);
  if (stat != CUBLAS_STATUS_SUCCESS)
  {
    printf("Cublas setMatrix error!!!\n");
    goto error;
  }

  {
    cublasOperation_t trans = CUBLAS_OP_N;
    double alpha = 1.0;
    double beta = 0.0;

    cublasDgemm(handle, trans, trans, n, n, n, &alpha, A, n, B, n, &beta, C, n);

    stat = cublasGetMatrix(n, n, sizeof(double), dC, n, C, n);

    if (stat != CUBLAS_STATUS_SUCCESS)
    {
      printf("Cublas data upload error!!!\n");
      goto error;
    }
  }

error:
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
  calc_blas(A, B, C, n);
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
