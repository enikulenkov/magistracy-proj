#include <stdio.h>
#include <stdlib.h>

#include <pagmo/pagmo.h>
#include "cga.h"
#include "ga_utils.h"

#define GA_ITERATIONS_NUM  500

double Clength=0.529177;   //http://en.wikipedia.org/wiki/Bohr_radius
double g_a0 = 3.8344/Clength; /* Lattice constant. http://en.wikipedia.org/wiki/Lattice_constant */

void do_work(double *init_coords, int atoms_cnt)
{
  int i;
  pagmo::algorithm::sga algo=pagmo::algorithm::sga(GA_ITERATIONS_NUM);
  cga prob = cga(atoms_cnt);
  pagmo::archipelago arch = pagmo::archipelago(algo, prob, 1, 20);
  pagmo::decision_vector init_x(atoms_cnt*3);

  for (i = 0; i < atoms_cnt*3; i++)
  {
    init_x[i] = g_a0 * init_coords[i];
  }

  /* Set initial decision vector for all individuals */
  pagmo::base_island_ptr island_copy = arch.get_island(0);
  
  for (i=0; i < island_copy->get_size(); i++)
  {
    island_copy->set_x(i, init_x);
  }

  arch.set_island(0, *island_copy);

  arch.evolve(1);
  arch.join();
 
  std::cout << arch.get_island(0)->human_readable();
}


int main()
{
  int total_atoms_cnt;
  double *init_coords;
  //char filename[] = "FCC.mat";
  char filename[] = "tetrahedron.mat";

  total_atoms_cnt = read_input(filename, &init_coords);

  if (total_atoms_cnt < 0)
  {
    DBG_LOG("Invalid input, exit...");
    return 1;
  }

  DBG_LOG("Loaded %d atoms' coordinates\n", total_atoms_cnt);

  do_work(init_coords, total_atoms_cnt);

  return 0;
}
