#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "timer.h"
#include "cblas.h"

//double Z[N][N] = {{10.0, 1.0, 2.0, 3.0, 4.0},
//				  {1.0, 9.0, -1.0, 2.0, -3.0},
//				  {2.0, -1.0, 7.0, 3.0, -5.0},
//				  {3.0, 2.0, 3.0, 12.0, -1.0},
//				  {4.0, -3.0, -5.0, -1.0, 15.0}};

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

	printf("Spur=%f\n", trace);

	timer_init(&timer);
	timer_start(&timer);
	for (i = 0; i < N - 2; i++)
	{
		memcpy(A0, A, N*N*sizeof(double)); // A0=A
		m = N - i - 1; //m = N - i + 1
		sigma = 0;
		for (j = 0; j < m; j++)
		{
			sigma += A[m+ j*N] * A[m + j*N]; // For[j = 1, j <= m - 1, j++, sigma = sigma + A[[m, j]] ^ 2];
		}
		for (j = 0; j < N; j++)
		{
			ut[j] = j <= m - 2 ? A[m + j*N] : 0; // For[j = 1, j <= N, j++, ut[[1, j]] = If[j <= m - 2, A[[m, j]], 0.0]];
		}
		ut[m - 1] = A[m+ (m - 1)*N] + sqrt(sigma);  //ut[[1, m - 1]] = A[[m, m - 1]] + Sqrt[sigma];
		
		double H = cblas_ddot(N, ut, 1, ut, 1) / 2; // {{H}} = 0.5*ut.u;
		for (j = 0; j < N; j++)
		{
			p[j] = cblas_ddot(N, A0, 1, ut, 1) / H; //p = A0.u/H;
		}

		double K = cblas_ddot(N, ut, 1, p, 1) / (2.0 * H); //{{ K }} = ut.p / (2.0*H);
		for (j = 0; j < N; j++)
		{
			q[j] = p[j] - K*ut[j]; //q = p - K*u;
		}
		memset(sub1, 0, N*N*sizeof(double));
		memset(sub2, 0, N*N*sizeof(double));

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, 1, 1.0, q, 1, ut, N, 0.0, sub1, N); // u.Transpose[q]
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, 1, 1.0, ut, 1, q, N, 0.0, sub2, N); // q.ut

		for (j = 0; j < N; j++)
		{
			for (k = 0; k < N; k++)
			{
				A[j + k*N] = A0[j + k*N] - sub1[j + k*N] - sub2[j + k*N];
			}
		}
	}
	timer_stop(&timer);
  printf("Elapsed time: %s", timer_diff_as_str(&timer));

	Trace(A, N);
	printf("\nSpur=%f\n", trace);

	return 0;
}
