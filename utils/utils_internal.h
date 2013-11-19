#ifndef UTILS_INTERNAL_H
#define UTILS_INTERNAL_H 1

#define MAT_PRMBL_TYPE_STR "type"
#define MAT_PRMBL_ROWS_STR "rows"
#define MAT_PRMBL_COLS_STR "columns"
#define MAT_PRMBL_NNZ_STR  "nnz"

#define PRMBL_FLD_DELIM    ": "

#define PRMBL_MAX_TYPE_LEN 64

typedef enum prmbl_fld_e
{
  PRMBL_FLD_NONE,  /* all expected fields are read, parsing is finished */
  PRMBL_FLD_TYPE,
  PRMBL_FLD_NNZ,
  PRMBL_FLD_ROWS,
  PRMBL_FLD_COLUMNS,
}
prmbl_fld_t;

#endif
