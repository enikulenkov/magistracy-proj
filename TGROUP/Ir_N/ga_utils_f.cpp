#include <stdlib.h>
#include <stdio.h>
#include <iostream>

#include <openbabel/obconversion.h>
#include <openbabel/mol.h>

extern "C" {
#include "utils.h"
#include "mat_file_parser.h"
}

/* Every atom have 3 coordinates */
#define COORDS_CNT    3

extern "C" int read_input_(char *test_file, double *A, int *na)
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
    load_matr(&mat_ctx, A);
  }

  deinit_mat_file_ctx(&mat_ctx);

  *na = (ret == 0) ? atoms_cnt : ret;
}

extern "C" int write_output_(char *out_filename, double *A, int *n)
{
  std::ofstream ofs(out_filename);
  OpenBabel::OBConversion ob(NULL, &ofs);
  OpenBabel::OBAtom atom;
  OpenBabel::OBMol mol;
  int i;

  ob.SetOutFormat("CML");

  /* Atom is Iridium */
  atom.SetAtomicNum(77);

  for (i = 0; i < *n; i++)
  {
    atom.SetVector(A[i*3], A[i*3+1], A[i*3+2]);
    mol.AddAtom(atom);
  }

  //for (i=0; i < *n; i++)
  //{
    //for (int j=i; j < *n; j++)
    //{
      //mol.AddBond(i+1, j+1, 0);
    //}
  //}

  ob.Write(&mol);
  ob.CloseOutFile();
}
