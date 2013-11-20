#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "utils_internal.h"

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
