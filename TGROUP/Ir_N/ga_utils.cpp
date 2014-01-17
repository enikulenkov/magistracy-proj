#include <stdlib.h>
#include <stdio.h>

extern "C" {
#include "utils.h"
#include "mat_file_parser.h"
}

/* Every atom have 3 coordinates */
#define COORDS_CNT    3

int read_input(char *test_file, double **A)
{
  int ret = 0;
  int atoms_cnt;
  mat_file_ctx_t mat_ctx;

  init_mat_file_ctx(&mat_ctx, test_file);

  ret = parse_mat_preamble(&mat_ctx);

  if ((ret != 0) || (mat_ctx.preamble.obj_type != DENSE_MATRIX)
      || (mat_ctx.preamble.cols_cnt != COORDS_CNT))
  {
    DBG_LOG("Error loading matrix \n");
    ret = -1;
  }
  else
  {
    atoms_cnt = mat_ctx.preamble.rows_cnt;
    *A = (double *)malloc(atoms_cnt*COORDS_CNT*sizeof(double));
    load_matr(&mat_ctx, *A);
  }

  deinit_mat_file_ctx(&mat_ctx);

  return (ret == 0) ? atoms_cnt : ret;
}