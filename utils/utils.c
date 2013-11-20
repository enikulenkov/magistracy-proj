#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "utils_internal.h"
#include "utils.h"

void random_matr(double *A, int n, int m)
{
  int i,j;

  for (i=0; i<n; i++)
  {
    for (j=0; j < m; j++)
    {
      A[ IDX(i,j,m) ] = (double)rand() / (double)RAND_MAX;
    }
  }
}


int check_vectors_for_eq(double *v1, double *v2, int n)
{
  int i;
  int ret = 1;
  double eps=0.0000001;

  for (i=0; i<n; i++)
  {
    if (fabs(v1[i] - v2[i]) > eps)
    {
      ret = 0;
      break;
    }
  }

  return ret;
}
