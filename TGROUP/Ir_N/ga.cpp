#include <stdio.h>
#include <stdlib.h>

#include <pagmo/pagmo.h>
#include "cga.h"
#include "ga_utils.h"

#define GA_ITERATIONS_NUM  10000

double Clength=0.529177;   //http://en.wikipedia.org/wiki/Bohr_radius
double g_a0 = 3.8344/Clength; /* Lattice constant. http://en.wikipedia.org/wiki/Lattice_constant */

void do_work(double *init_coords, int atoms_cnt, double *res_coords, double *min_value)
{
  int i;
  pagmo::algorithm::sga algo=pagmo::algorithm::sga(GA_ITERATIONS_NUM);
  cga prob = cga(atoms_cnt);
  pagmo::archipelago arch = pagmo::archipelago(algo, prob, 1, 20);
  pagmo::decision_vector init_x(atoms_cnt*3);
  pagmo::decision_vector best_decision_vec;

  for (i = 0; i < atoms_cnt*3; i++)
  {
    init_x[i] = g_a0 * init_coords[i];
  }

  /* Set initial decision vector for all individuals */
  pagmo::base_island_ptr island_copy = arch.get_island(0);
  
  for (i=0; i < island_copy->get_size(); i+=5)
  {
    island_copy->set_x(i, init_x);
  }

  arch.set_island(0, *island_copy);

  arch.evolve(1);
  arch.join();
 
  /* std::cout << arch.get_island(0)->human_readable(); */

  best_decision_vec = arch.get_island(0)->get_population().champion().x;

  for (i=0; i < best_decision_vec.size(); i++)
  {
    res_coords[i] = best_decision_vec[i];
  }

  *min_value = arch.get_island(0)->get_population().champion().f[0];
}


int main(int argc, char *argv[])
{
  int total_atoms_cnt;
  double *init_coords;
  double *res;
  double min;
  char *in_file;
  char *out_file;
  char out_fname[] = "out.cml";
  int i;

  if (argc < 3)
  {
    fprintf(stderr, "Usage: %s {input_file.mat} {output_file.mat}", argv[0]);
    return 1;
  }

  in_file = argv[1];
  out_file = argv[2];

  total_atoms_cnt = read_input(in_file, &init_coords);

  if (total_atoms_cnt < 0)
  {
    DBG_LOG("Invalid input, exit...");
    return 1;
  }

  DBG_LOG("Loaded %d atoms' coordinates\n", total_atoms_cnt);

  res = (double *)malloc(total_atoms_cnt * 3 * sizeof(double));

  do_work(init_coords, total_atoms_cnt, res, &min);

  write_output(out_fname, res, total_atoms_cnt);

  for (i=0; i < total_atoms_cnt*3; i++)
  {
    res[i] = res[i]/g_a0;
  }

  save_matr_to_mat_file(out_file, res, total_atoms_cnt, 3);

  printf("Min value: %lf\n", min);

  free(res);

  return 0;
}
