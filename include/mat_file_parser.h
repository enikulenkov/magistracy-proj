#ifndef MAT_FILE_PARSER_H
#define MAT_FILE_PARSER_H 1

#ifdef DEBUG
#define  DBG_LOG printf
#else
#define  DBG_LOG
#endif

typedef enum mat_obj_type_e
{
  SPARSE_MATRIX,
  DENSE_MATRIX,
}
mat_obj_type_t;


typedef struct mat_file_preamble_s
{
  mat_obj_type_t  obj_type;
  int             rows_cnt;
  int             cols_cnt;
  int             nnz;       /* non-zero values in sparse matrix */
}
mat_file_preamble_t;

typedef struct mat_file_ctx_s
{
  FILE *f;
  char *cur_line;
  char cur_line_no;
  mat_file_preamble_t preamble;
}
mat_file_ctx_t;

int init_mat_file_ctx(mat_file_ctx_t *ctx, char *file_name);
int parse_mat_preamble(mat_file_ctx_t *ctx);
void deinit_mat_file_ctx(mat_file_ctx_t *ctx);

int load_sparse_matrix_coo(mat_file_ctx_t *ctx,
    int *row_indeces,
    int *col_indeces,
    double *values);

int load_sparse_matrix_as_dense(mat_file_ctx_t *ctx, double *A);
int load_dense_matrix(mat_file_ctx_t *ctx, double *A, int n);
int load_vector(mat_file_ctx_t *ctx, double *v);

#endif /* MAT_FILE_PARSER_H */
