#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mat_file_parser.h"
#include "utils_internal.h"

int init_mat_file_ctx(mat_file_ctx_t *ctx, char *file_name)
{
  memset(ctx, 0, sizeof(mat_file_ctx_t));

  ctx->f = fopen(file_name, "r");

  if (ctx->f == NULL)
  {
    DBG_LOG("Error opening file %s\n", file_name);
    return -1;
  }

  return 0;
}


void deinit_mat_file_ctx(mat_file_ctx_t *ctx)
{
  if (ctx->cur_line != NULL)
  {
    free(ctx->cur_line);
  }

  if (ctx->f != NULL)
  {
    fclose(ctx->f);
  }
}


int parse_mat_preamble(mat_file_ctx_t *ctx)
{
  int          ret = 0;
  ssize_t      read;
  size_t       len = 0;
  prmbl_fld_t  exp_prmbl_fld; 

  exp_prmbl_fld = PRMBL_FLD_TYPE;

  while (((read = getline(&ctx->cur_line, &len, ctx->f)) != -1)
           && (ret == 0)
           && (!strncmp(ctx->cur_line, "# ", 2)))
  {
    ctx->cur_line_no++;

    /* ctx->cur_line is a comment, skip hash and space after it */
    read -= 2;

    if (read > 0)
    {
      char *comment = &ctx->cur_line[2];

      switch(exp_prmbl_fld)
      {
        case PRMBL_FLD_TYPE:
        {
          if (!strncmp(comment, MAT_PRMBL_TYPE_STR, strlen(MAT_PRMBL_TYPE_STR)))
          {
            if (!strcmp(&comment[strlen(MAT_PRMBL_TYPE_STR)], PRMBL_FLD_DELIM"matrix\n"))
            {
              ctx->preamble.obj_type = DENSE_MATRIX;
              exp_prmbl_fld = PRMBL_FLD_ROWS;
            }
            else if (!strcmp(&comment[strlen(MAT_PRMBL_TYPE_STR)], PRMBL_FLD_DELIM"sparse matrix\n"))
            {
              ctx->preamble.obj_type = SPARSE_MATRIX;
              exp_prmbl_fld = PRMBL_FLD_NNZ;
            }
            else
            {
              DBG_LOG("Error: Unknown math object type! (ctx->cur_line %d)\n", ctx->cur_line_no);
              ret = -1;
            }
          }
        }
        break;

        case PRMBL_FLD_NNZ:
        {
          if (!strncmp(comment, MAT_PRMBL_NNZ_STR, strlen(MAT_PRMBL_NNZ_STR)))
          {
            sscanf(&comment[strlen(MAT_PRMBL_NNZ_STR)], PRMBL_FLD_DELIM"%d\n", &ctx->preamble.nnz);
            exp_prmbl_fld = PRMBL_FLD_ROWS;
          }
        }
        break;

        case PRMBL_FLD_ROWS:
        {
          if (!strncmp(comment, MAT_PRMBL_ROWS_STR, strlen(MAT_PRMBL_ROWS_STR)))
          {
            sscanf(&comment[strlen(MAT_PRMBL_ROWS_STR)], PRMBL_FLD_DELIM"%d\n", &ctx->preamble.rows_cnt);
            exp_prmbl_fld = PRMBL_FLD_COLUMNS;
          }
        }
        break;

        case PRMBL_FLD_COLUMNS:
        {
          if (!strncmp(comment, MAT_PRMBL_COLS_STR, strlen(MAT_PRMBL_COLS_STR)))
          {
            sscanf(&comment[strlen(MAT_PRMBL_COLS_STR)], PRMBL_FLD_DELIM"%d\n", &ctx->preamble.cols_cnt);
            exp_prmbl_fld = PRMBL_FLD_NONE;
          }
        }
        break;
      }
    }
  }

  if (exp_prmbl_fld != PRMBL_FLD_NONE)
  {
    DBG_LOG("Error: Not all of the expected fields are found, curr exp_field %d!\n", exp_prmbl_fld);
    ret = -1;
  }

  return ret;
}


/**
 * @brief Loads N*N sparse matrix from .mat file
 *
 * In .mat file all matrix elements have the form:
 * "x y value"
 * x and y enumeration starts with 1, not 0
 *
 * Function loads matrix in contiguous memory block, starting
 * with address of passed matrix A.
 *
 * @param A         Pointer to the square matrix
 * @param n         Size of the matrix A
 * @param file_name File to read matrix from
 *
 * @return N - size of the matrix that was loaded from the file 
 *         -1  on error
 */
int load_sparse_matrix_as_dense(mat_file_ctx_t *ctx, double *A)
{
  int     ret = 0;
  int     i = 0;
  ssize_t read;
  size_t  len = 0;

  /* First line with data is stored in ctx->cur_line.
   * This line was read during parsing preamble */
  do
  {
    int row_idx, col_idx;
    double value;

    sscanf(ctx->cur_line, "%d %d %lf\n", &row_idx, &col_idx, &value);

    if ((row_idx > ctx->preamble.rows_cnt) || (col_idx > ctx->preamble.cols_cnt))
    {
      DBG_LOG("Element index is out of range, ctx->cur_line %d\n", ctx->cur_line_no);
      ret = -1;
    }
    else
    {
      A[(row_idx-1)*ctx->preamble.rows_cnt + (col_idx-1)] = value; 
    }

    read = getline(&ctx->cur_line, &len, ctx->f);

    if (read != -1)
    {
      ctx->cur_line_no++;
    }
    else
    {
      DBG_LOG("Error: too few rows in file\n");
      ret = -1;
    }

    i++;
  }
  while ((i < ctx->preamble.nnz) && (ret == 0));

  return ret;
}

/**
 * @brief Loads vector from .mat file
 *
 * In .mat file vector elements written one by line.
 *
 * @param v         Pointer to the vector 
 * @param n         Size of the vector
 * @param file_name File to read vector from
 *
 * @return N - size of the vector that was loaded from the file 
 *         -1  on error
 */
int load_vector(mat_file_ctx_t *ctx, double *v)
{
  int     i = 0;
  int     ret = 0;
  ssize_t read;
  size_t  len = 0;

  if (ctx->preamble.cols_cnt > 1)
  {
    return -1;
  }

  /* First line with data is stored in ctx->cur_line.
   * This line was read during parsing preamble */
  do
  {
    double value;

    sscanf(ctx->cur_line, "%lf\n", &value);

    v[i] = value;

    read = getline(&ctx->cur_line, &len, ctx->f);

    if (read != -1)
    {
      ctx->cur_line_no++;
    }
    else
    {
      DBG_LOG("Error: too few rows in file\n");
      ret = -1;
    }

    i++;
  }
  while ((i < ctx->preamble.rows_cnt) && (ret == 0));

  return ret;
}


int load_sparse_matrix_coo(mat_file_ctx_t *ctx,
    int *row_indeces,
    int *col_indeces,
    double *values)
{
  int     ret = 0;
  int     i = 0;
  ssize_t read;
  size_t  len = 0;

  /* First line with data is stored in ctx->cur_line.
   * This line was read during parsing preamble */
  do
  {
    int row_idx, col_idx;
    double value;

    sscanf(ctx->cur_line, "%d %d %lf\n", &row_idx, &col_idx, &value);

    if ((row_idx > ctx->preamble.rows_cnt) || (col_idx > ctx->preamble.cols_cnt))
    {
      DBG_LOG("Element index is out of range, ctx->cur_line %d\n", ctx->cur_line_no);
      ret = -1;
    }
    else
    {
      row_indeces[i] = row_idx - 1;
      col_indeces[i] = col_idx - 1;
      values[i] = value;
    }

    read = getline(&ctx->cur_line, &len, ctx->f);

    if (read != -1)
    {
      ctx->cur_line_no++;
    }
    else
    {
      DBG_LOG("Error: too few rows in file\n");
      ret = -1;
    }

    i++;
  }
  while ((i < ctx->preamble.nnz) && (ret == 0));

  return ret;
}



