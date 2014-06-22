#include <stdio.h>
#include <stdlib.h>

#include <pagmo/pagmo.h>
#include "sutton_chen_pot.h"
#include "ga_utils.h"
#include "ini.h"

#define GA_ITERATIONS_NUM  100

double Clength=0.529177;   //http://en.wikipedia.org/wiki/Bohr_radius
double g_a0 = 3.8344/Clength; /* Lattice constant. http://en.wikipedia.org/wiki/Lattice_constant */

typedef struct
{
  int generations_number;
  int clusters_count;
  double crossover_rate;
  double binom_rate;
  double initial_rand_clusters;
  double min_atom_dist;
  double mut_move_rate;
  double mut_rotate_rate;
  double mut_replace_rate;
  int elitism_number;
  char *selection_type;
  char *crossover_type;
  double bfgs_step_size;
  double bfgs_tol;
}
app_config_t;

app_config_t config;

static int config_handler(void* user, const char* section, const char* name,
                   const char* value)
{
    app_config_t* pconfig = (app_config_t*)user;

    #define MATCH(n) strcmp(section, "common") == 0 && strcmp(name, n) == 0
    if (MATCH("generations_number")) {
        pconfig->generations_number = atoi(value);
    } else if (MATCH("clusters_count")) {
        pconfig->clusters_count = atoi(value);
    } else if (MATCH("crossover_rate")) {
        pconfig->crossover_rate = atof(value);
    } else if (MATCH("binom_rate")) {
        pconfig->binom_rate = atof(value);
    } else if (MATCH("initial_rand_clusters")) {
        pconfig->initial_rand_clusters = atof(value);
    } else if (MATCH("min_atom_dist")) {
        pconfig->min_atom_dist = atof(value);
    } else if (MATCH("mut_move_rate")) {
        pconfig->mut_move_rate = atof(value);
    } else if (MATCH("mut_rotate_rate")) {
        pconfig->mut_rotate_rate = atof(value);
    } else if (MATCH("mut_replace_rate")) {
        pconfig->mut_replace_rate = atof(value);
    } else if (MATCH("elitism_number")) {
        pconfig->elitism_number = atoi(value);
    } else if (MATCH("selection_type")) {
        pconfig->selection_type = strdup(value);
    } else if (MATCH("crossover_type")) {
        pconfig->crossover_type = strdup(value);
    } else if (MATCH("bfgs_step_size")) {
        pconfig->bfgs_step_size = atof(value);
    } else if (MATCH("bfgs_tol")) {
        pconfig->bfgs_tol = atof(value);
    } else {
        return 0;  /* unknown section/name, error */
    }
    return 1;
}


void find_min_max(double *a, int n, double *min, double *max)
{
  int i;
  *min = a[0];
  *max = a[0];

  for (i=1; i < n; i++)
  {
    if (a[i] > *max)
    {
      *max = a[i];
    }
    else if (a[i] < *min)
    {
      *min =a[i];
    }
  }
}


void do_work(double *init_coords, int atoms_cnt, double *res_coords, double *min_value)
{
  int i;
  sutton_chen_pot prob = sutton_chen_pot(atoms_cnt);
  pagmo::decision_vector init_x(atoms_cnt*3);
  pagmo::decision_vector best_decision_vec;
  pagmo::algorithm::birmingham_ga::mutation mutations[3];
  int mut_count = 0;
  pagmo::algorithm::birmingham_ga::selection::type sel_type;
  pagmo::algorithm::birmingham_ga::crossover::type cross_type;
  double min_coord;
  double max_coord;

  if (config.mut_move_rate > 0)
  {
    mutations[mut_count].type = pagmo::algorithm::birmingham_ga::mutation::MOVE;
    mutations[mut_count++].probability = config.mut_move_rate;
  }

  if (config.mut_rotate_rate > 0)
  {
    mutations[mut_count].type = pagmo::algorithm::birmingham_ga::mutation::ROTATE;
    mutations[mut_count++].probability = config.mut_rotate_rate;
  }

  if (config.mut_replace_rate > 0)
  {
    mutations[mut_count].type = pagmo::algorithm::birmingham_ga::mutation::REPLACE;
    mutations[mut_count++].probability = config.mut_replace_rate;
  }

  if (!strcmp(config.selection_type, "TOURNAMENT"))
  {
    sel_type = pagmo::algorithm::birmingham_ga::selection::TOURNAMENT;
  }
  else if (!strcmp(config.selection_type, "ROULETTE"))
  {
    sel_type = pagmo::algorithm::birmingham_ga::selection::ROULETTE;
  }

  if (!strcmp(config.crossover_type, "BINOMIAL"))
  {
    cross_type = pagmo::algorithm::birmingham_ga::crossover::BINOMIAL;
  }
  else if (!strcmp(config.crossover_type, "CUT_AND_SPLICE"))
  {
    cross_type = pagmo::algorithm::birmingham_ga::crossover::CUT_AND_SPLICE;
  }

  for (i = 0; i < atoms_cnt*3; i++)
  {
    init_x[i] = g_a0 * init_coords[i];
  }

  find_min_max(&init_x[0], atoms_cnt*3, &min_coord, &max_coord);

  if (abs(min_coord) > max_coord)
  {
    max_coord = abs(min_coord);
  }

  
  pagmo::algorithm::birmingham_ga algo = pagmo::algorithm::birmingham_ga(
      config.generations_number,
      config.crossover_rate,
      config.binom_rate,
      config.min_atom_dist,
      mutations,
      mut_count,
      config.elitism_number,
      sel_type,
      cross_type,
      max_coord,
      config.bfgs_step_size,
      config.bfgs_tol);

  prob.set_lb(-2*max_coord);
  prob.set_ub(2*max_coord);

  pagmo::archipelago arch = pagmo::archipelago(algo, prob, 1, config.clusters_count);

  /* Set initial decision vector for all individuals */
  pagmo::base_island_ptr island_copy = arch.get_island(0);
  
  for (i=0; i < island_copy->get_size(); i++)
  {
    if (i < island_copy->get_size() * config.initial_rand_clusters - 1)
    {
      static_cast<const pagmo::algorithm::birmingham_ga>(algo).randomize_cluster(init_x);
    }
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
  char out_fname_orig[] = "out_orig.cml";
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

  if (ini_parse("gen_alg_prefs.ini", config_handler, &config) < 0) {
    DBG_LOG("Can't load 'test.ini'\n");
    return 1;
  }

  DBG_LOG("Loaded %d atoms' coordinates\n", total_atoms_cnt);

  write_output(out_fname_orig, init_coords, total_atoms_cnt);

  res = (double *)malloc(total_atoms_cnt * 3 * sizeof(double));

  do_work(init_coords, total_atoms_cnt, res, &min);

  write_output(out_fname, res, total_atoms_cnt);

  for (i=0; i < total_atoms_cnt*3; i++)
  {
    res[i] = res[i]/g_a0;
  }

  save_matr_to_mat_file(out_file, res, total_atoms_cnt, 3);

  printf("Min value: %lf\n", min);

  free(init_coords);
  free(res);

  return 0;
}
