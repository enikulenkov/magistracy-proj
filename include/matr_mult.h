#ifndef MATR_MULT_H
#define MATR_MULT_H 1

void matr_vector_mult_(double *h_A, double *h_B, double *h_C, int *numElements);

void matr_vector_mult_sparse(int *h_A_row_indeces, int *h_A_col_indeces, double *h_A_values, int nnz,
    double *h_B, double *h_C, int *numElements);

#endif //MATR_MULT_H
