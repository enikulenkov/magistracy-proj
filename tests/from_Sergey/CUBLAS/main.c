#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "timer.h"
#include "cublas.h"
#include "cblas.h"

//double Z[N][N] = {{10.0, 1.0, 2.0, 3.0, 4.0},
//          {1.0, 9.0, -1.0, 2.0, -3.0},
//          {2.0, -1.0, 7.0, 3.0, -5.0},
//          {3.0, 2.0, 3.0, 12.0, -1.0},
//          {4.0, -3.0, -5.0, -1.0, 15.0}};

double Trace(double *matrix, int n)
{
  int i;

  double trace = 0;
  for (i = 0; i < n; i++)
  {
    trace += matrix[i*n + i];
  }
  return trace;
}

void RamdomMatr(double *matrix, int n)
{
  int i,j;

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      if (i == j)
      {
        double z = (i + 3.0) / (j + 1);
        matrix[i*n + j] = z;
      }
      else
      {
        matrix[i*n + j] = 7;
      }
    }
  }
}

int main(int argc, char *argv[])
{
  timer_ctx_t timer;
  timer_ctx_t timer2;
  int i,j,k;

  if (argc < 2)
  {
    return -1;
  }

  int N = atoi(argv[1]);
  printf("N=%d\n", N);

  double *A = (double*) malloc(N*N*sizeof(double));
  RamdomMatr(A, N);
  //A = *Z;

  double *A0 = (double*) malloc(N*N*sizeof(double));
  double *ut = (double*) malloc(N*sizeof(double));
  double *p = (double*) malloc(N*sizeof(double));
  double *q = (double*) malloc(N*sizeof(double));
  double sigma, trace = Trace(A, N);
  int m;
  double *sub1 = (double*) malloc(N*N*sizeof(double));
  double *sub2 = (double*) malloc(N*N*sizeof(double));

  cublasStatus status;
  double *cuda_q, *cuda_ut, *cuda_p, *cuda_sub1, *cuda_sub2, *cuda_A0;

  printf("Spur=%f\n", trace);

  // CUDA
  status = cublasAlloc(N * 1, sizeof(double), (void **) &cuda_q);
  status = cublasAlloc(N * 1, sizeof(double), (void **) &cuda_ut);
  status = cublasAlloc(N * 1, sizeof(double), (void **) &cuda_p);
  status = cublasAlloc(N * N, sizeof(double), (void **) &cuda_sub1);
  status = cublasAlloc(N * N, sizeof(double), (void **) &cuda_sub2);
  //status = cublasAlloc(N * 1, sizeof(double), (void **) &cuda_A0);
  status = cublasAlloc(N * N, sizeof(double), (void **) &cuda_A0);
  // CUDA

  timer_init(&timer);
  timer_start(&timer);
  for (i = 0; i < N - 2; i++)
  {
    
    memcpy(A0, A, N*N*sizeof(double)); // A0=A
    m = N - i - 1; //m = N - i + 1
    sigma = 0;
    for (j = 0; j < m; j++)
    {
      sigma += A[m + j*N] * A[m + j*N]; // For[j = 1, j <= m - 1, j++, sigma = sigma + A[[m, j]] ^ 2];
    }
    for (j = 0; j < N; j++)
    {
      ut[j] = j <= m - 2 ? A[m + j*N] : 0; // For[j = 1, j <= N, j++, ut[[1, j]] = If[j <= m - 2, A[[m, j]], 0.0]];
    }
    ut[m - 1] = A[m + (m - 1)*N] + sqrt(sigma);  //ut[[1, m - 1]] = A[[m, m - 1]] + Sqrt[sigma];
    
    // !CUDA very slow!
    //status = cublasSetMatrix(N, 1, sizeof(double), ut, N, cuda_ut, N);
    //double H = cublasSdot(N, cuda_ut, 1, cuda_ut, 1) / 2; // {{H}} = 0.5*ut.u;
    double H = cblas_ddot(N, ut, 1, ut, 1) / 2; // {{H}} = 0.5*ut.u;

    //status = cublasSetMatrix(N, N, sizeof(double), A0, N, cuda_A0, N);
    for (j = 0; j < N; j++)
    {
      // !CUDA very slow!
      //status = cublasSetMatrix(N, 1, sizeof(double), &A0[j*N], N, cuda_A0, N);
      //p[j] = cublasDdot(N, &cuda_A0[j*N], 1, cuda_ut, 1) / H; //p = A0.u/H;                // TODO replace cuda_A0 -> &cuda_A0[j*N]
      p[j] = cblas_ddot(N, &A0[j*N], 1, ut, 1) / H; //p = A0.u/H;
    }
    
    // !CUDA very slow!
    //status = cublasSetMatrix(N, 1, sizeof(double), p, N, cuda_p, N);
    //double K = cublasDdot(N, cuda_ut, 1, cuda_p, 1) / (2.0 * H); //{{ K }} = ut.p / (2.0*H);
    double K = cblas_ddot(N, ut, 1, p, 1) / (2.0 * H); //{{ K }} = ut.p / (2.0*H);
    for (j = 0; j < N; j++)
    {
      q[j] = p[j] - K*ut[j]; //q = p - K*u;
    }
    memset(sub1, 0, N*N*sizeof(double));
    memset(sub2, 0, N*N*sizeof(double));
    
    status = cublasSetMatrix(N, 1, sizeof(double), q, N, cuda_q, N);
    status = cublasSetMatrix(N, 1, sizeof(double), ut, N, cuda_ut, N);

    timer_init(&timer2);
    timer_start(&timer2);
    cublasDgemm('n', 'n', N, N, 1, 1.0, cuda_ut, N, cuda_q, 1, 0.0, cuda_sub1, N); // u.Transpose[q]
    timer_stop(&timer2);
    printf("Time u.Transpose[q]: %s\n", timer_diff_as_str(&timer2));
    timer_start(&timer2);
    cublasDgemm('n', 'n', N, N, 1, 1.0, cuda_q, N, cuda_ut, 1, 0.0, cuda_sub2, N); // q.ut
    timer_stop(&timer2);
    printf("Time q.ut: %s\n", timer_diff_as_str(&timer2));

    for (j = 0; j < N; j++)
    {
      for (k = 0; k < N; k++)
      {
        A[j + k*N] = A0[j + k*N] - sub1[j + k*N] - sub2[j + k*N];
      }
    }
  }

  timer_stop(&timer);
  printf("Total time: %s\n", timer_diff_as_str(&timer));

  Trace(A, N);
  printf("\nSpur=%f\n", trace);

  return 0;
}
