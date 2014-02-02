#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <utils.h>
#include "utils_internal.h"

static void write_matr_preambule(FILE *f, int rows, int cols)
{
  time_t t = time(NULL);
  struct tm tm = *localtime(&t);

  fprintf(f, "# Created on %d-%d-%d %d:%d:%d\n", tm.tm_year + 1900, tm.tm_mon + 1,
      tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
  /* TODO: Fix hardcoded value */
  fprintf(f, "# name: B\n");
  fprintf(f, "# type: matrix\n");
  fprintf(f, "# rows: %d\n", rows);
  fprintf(f, "# columns: %d\n", cols);
}


int save_matr_to_mat_file(char *file_name, double *matr, int rows, int cols)
{
  int i,j;
  FILE *f = fopen(file_name, "w");

  if (f == NULL)
  {
    DBG_LOG("Error opening file %s\n", file_name);
    return -1;
  }

  write_matr_preambule(f, rows, cols);

  for (i=0; i < rows; i++)
  {
    for (j=0; j < cols-1; j++)
    {
      fprintf(f, "%lf ", matr[i*cols+j]);
    }

    fprintf(f, "%lf\n", matr[i*cols+cols-1]);
  }

  close (f);

  return 0;
}
