#ifndef GA_UTILS_H
#define GA_UTILS_H 1

extern "C" {
#include "utils.h"
#include "mat_file_writer.h"
}

int read_input(char *test_file, double **A);

int write_output(char *out_filename, double *A, int n);

#endif /* GA_UTILS_H */
