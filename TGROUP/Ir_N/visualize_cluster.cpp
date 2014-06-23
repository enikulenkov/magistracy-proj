#include <stdio.h>
#include <stdlib.h>

#include "ga_utils.h"

int main(int argc, char *argv[])
{
  int total_atoms_cnt;
  double *coords;
  char *in_file;
  char *out_file;

  if (argc < 3)
  {
    fprintf(stderr, "Usage: %s {input_file.mat} {output_file.cml}", argv[0]);
    return 1;
  }

  in_file = argv[1];
  out_file = argv[2];

  total_atoms_cnt = read_input(in_file, &coords);

  if (total_atoms_cnt < 0)
  {
    DBG_LOG("Invalid input, exit...");
    return 1;
  }

  DBG_LOG("Loaded %d atoms' coordinates\n", total_atoms_cnt);

  write_output(out_file, coords, total_atoms_cnt);

  printf("Result is written to %s\n", out_file);

  free(coords);

  return 0;
}
