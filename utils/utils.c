#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "utils_internal.h"

/*#define MAX_PREAMBLE_FLD_NAME_LEN    sizeof("columns")*/
/*#define MAX_PREAMBLE_FLD_SSCANF_STR  4*/

/*typedef struct preamble_field_s*/
/*{*/
  /*char name[MAX_PREAMBLE_FLD_NAME_LEN];*/
  /*int  name_len;*/
  /*char sscanf_str[MAX_PREAMBLE_FLD_SSCANF_STR];*/
  /*int  sscanf_str_len;*/
/*}*/
/*preamble_field_t;*/


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
