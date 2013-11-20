#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "mat_file_parser.h"
#include "matr_mult.h"

#define  TEST_FILE  "./test1/full.mat"

#ifdef USE_DENSE_MATRIX
int load_data(char *test_file, double **A, double **B, double **res)
#else
int load_data(char *test_file, int **A_row_indeces, int **A_col_indeces,
    double **A_values, int *nnz, double **B, double **res)
#endif
{
  int ret = 0;
  int size;
  *A_row_indeces = *A_col_indeces = NULL;
  *A_values = *B = *res = NULL;

  mat_file_ctx_t mat_ctx;

  init_mat_file_ctx(&mat_ctx, test_file);

  /* Load matrix A */
  ret = parse_mat_preamble(&mat_ctx);

  if ((ret != 0)
      || (mat_ctx.preamble.obj_type != SPARSE_MATRIX)
      || (mat_ctx.preamble.rows_cnt != mat_ctx.preamble.cols_cnt))
  {
    DBG_LOG("Error loading matrix A\n");
    ret = -1;
    goto done;
  }
  else
  {
    size = mat_ctx.preamble.rows_cnt;
    DBG_LOG("matrix size = %d\n", size);
#ifdef USE_DENSE_MATRIX
    *A = (double *)malloc(size*size*sizeof(double));
    load_sparse_matrix_as_dense(&mat_ctx, *A);
#else
    *A_row_indeces = (int *)malloc(mat_ctx.preamble.nnz*sizeof(double));
    *A_col_indeces = (int *)malloc(mat_ctx.preamble.nnz*sizeof(double));
    *A_values = (double *)malloc(mat_ctx.preamble.nnz*sizeof(double));
    *nnz = mat_ctx.preamble.nnz;
    load_sparse_matrix_coo(&mat_ctx, *A_row_indeces, *A_col_indeces, *A_values);
#endif
  }

  /* Load vector B */
  ret = parse_mat_preamble(&mat_ctx);

  if ((ret != 0)
      || (mat_ctx.preamble.obj_type != DENSE_MATRIX)
      || (mat_ctx.preamble.cols_cnt != 1)
      || (mat_ctx.preamble.rows_cnt != size))
  {
    DBG_LOG("Error loading vector B\n");
    ret = -1;
    goto done;
  }
  else
  {
    *B = (double *)malloc(size*sizeof(double));
    load_vector(&mat_ctx, *B);
  }

  /* Load result vector */
  ret = parse_mat_preamble(&mat_ctx);

  if ((ret != 0)
      || (mat_ctx.preamble.obj_type != DENSE_MATRIX)
      || (mat_ctx.preamble.cols_cnt != 1)
      || (mat_ctx.preamble.rows_cnt != size))
  {
    DBG_LOG("Error loading result vector\n");
    ret = -1;
    goto done;
  }
  else
  {
    *res = (double *)malloc(size*sizeof(double));
    load_vector(&mat_ctx, *res);
  }

  deinit_mat_file_ctx(&mat_ctx);

done:
  if (ret != 0)
  {
    deinit_mat_file_ctx(&mat_ctx);
    free(*A_row_indeces);
    free(*A_col_indeces);
    free(*A_values);
    free(*B);
    free(*res);
  }

  return (ret == 0) ? size : ret;
}


int main()
{
  int size, nnz;
  int check_res;
  double *vector, *res, *answer;

#ifdef USE_DENSE_MATRIX
  double *matr;
  size = load_data(TEST_FILE, &matr, &vector, &answer);
#else
  int *A_row_indeces, *A_col_indeces;
  double *A_values;
  size = load_data(TEST_FILE, &A_row_indeces, &A_col_indeces, &A_values, &nnz, &vector, &answer);
#endif

  if (size < 0)
  {
    printf("Error: can't load input data!\n");
    return 1;
  }

  res = (double *)malloc(size*sizeof(double));

#ifdef USE_DENSE_MATRIX
  matr_vector_mult_(matr, vector, res, &size);
#else
  matr_vector_mult_sparse(A_row_indeces, A_col_indeces, A_values, nnz, vector, res, &size);
#endif

  check_res = check_vectors_for_eq(res, answer, size);

  if (check_res)
  {
    printf("Result is correct\n");
  }
  else
  {
    printf("Result is incorrect!!!\n");
  }

  return 0;
}
