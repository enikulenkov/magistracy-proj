#include <stdio.h>
#include <stdlib.h>
#include <pagmo/pagmo.h>
#include <gtest/gtest.h>
#include "ga_utils.h"
#include "sutton_chen_pot.h"
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

double Clength=0.529177;   //http://en.wikipedia.org/wiki/Bohr_radius
double g_a0 = 3.8344/Clength; /* Lattice constant. http://en.wikipedia.org/wiki/Lattice_constant */

typedef struct numdiff_params_s
{
  sutton_chen_pot *prob;
  pagmo::decision_vector x;
  int coord;
}
numdiff_params_t;

double objfun_numdiff_wrapper(double x, void *params)
{
  numdiff_params_t *par = (numdiff_params_t *)params;
  pagmo::fitness_vector f = pagmo::fitness_vector(1);

  par->x[par->coord] = x;
  par->prob->objfun_impl_test(f, par->x);

  return f[0];
}

void gsl_calc_deriv(pagmo::decision_vector &coords, double *abs_err, double *res, sutton_chen_pot *prob)
{
  gsl_function F;
  numdiff_params_t params;
  double step_size = 1e-8;

  params.prob = prob;
  params.x = coords;

  F.function = objfun_numdiff_wrapper;
  F.params = (void *)&params;

  for (int i = 0; i < prob->get_dimension(); ++i)
  {
    params.coord = i;
    gsl_deriv_central(&F,coords[i],step_size,&res[i],&abs_err[i]);
  }
}

int main()
{
  const int atoms_cnt = 50;
  char test_file[]="50res.mat";
  //const int atoms_cnt = 2;
  //char test_file[]="2res.mat";
  double *coords;
  pagmo::decision_vector x = pagmo::decision_vector(atoms_cnt*3);
  double abs_err[atoms_cnt * 3];
  double res[atoms_cnt * 3];
  double df[atoms_cnt * 3];

  read_input(test_file, &coords);

  for (int i = 0; i < x.size(); i++)
  {
    x[i] = coords[i];
  }

  sutton_chen_pot prob = sutton_chen_pot(atoms_cnt);

  gsl_calc_deriv(x, abs_err, res, &prob);

  prob.d_objfun(coords, atoms_cnt*3, df);

  return 0;
}
