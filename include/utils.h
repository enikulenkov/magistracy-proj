#ifndef UTILS_H
#define UTILS_H 1

#define IDX(_r, _c, _n) ((_r)*(_n) + (_c))

int check_vectors_for_eq(double *v1, double *v2, int n);
void random_matr(double *A, int n, int m);

#endif /* UTILS_H */
