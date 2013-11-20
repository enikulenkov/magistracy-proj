#include <stdlib.h>
#include <stdio.h>
#include <f77blas.h>
#include "timer.h"

int g_sizes[] = {1000, 2000, 3500, 5000};

#define SIZES_COUNT sizeof(g_sizes)/sizeof(g_sizes[0])

void run_test(int n)
{
  double *A = (double *)malloc(n*n*sizeof(double));
  double *B = (double *)malloc(n*n*sizeof(double));
  double *C = (double *)malloc(n*n*sizeof(double));
  double alpha = 1.0;
  double beta = 0;
  char trans = 'N';
  timer_ctx_t timer;
 
  timer_init(&timer);

  random_matr(A, n, n);
  random_matr(B, n, n);

  timer_start(&timer);
  dgemm_(&trans, &trans, &n, &n, &n, &alpha, A, &n, B, &n, &beta, C, &n);
  timer_stop(&timer);

  printf("Size: %d\n", n);
  printf("Time: %s\n", timer_diff_as_str(&timer));

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
