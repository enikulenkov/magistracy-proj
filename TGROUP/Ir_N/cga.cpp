/*
 * Cluster genetic algorithm
 */
#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <string>

#include <pagmo/pagmo.h>
#include "cga.h"

#define  COORDS_CNT  3
#define  SQR(_a) ((_a)*(_a))
#define  M_PARAM  6
#define  N_PARAM  13
#define  C_PARAM  224.815
//#define  EPS_PARAM  3.7674E-3
#define  CENERGY  27.2116
#define  CLENGTH  0.529177
#define  EPS_PARAM  3.7674E-3/CENERGY
#define  A_PARAM  3.8344/CLENGTH

cga::cga(int n):pagmo::problem::base(n*3)
{
  // Set bounds.
  set_lb(-10.0);
  set_ub(10.0);
}


/* Clone method */
pagmo::problem::base_ptr cga::clone() const
{
  return pagmo::problem::base_ptr(new cga(*this));
}


void cga::objfun_impl(pagmo::fitness_vector &f, const pagmo::decision_vector &x) const
{
  pagmo_assert(f.size() == 1);
  std::vector<double>::size_type n = x.size();
  pagmo_assert(n % 3 == 0);
  int atoms_cnt = n/3;
  double *dists;
  int i = 0;
  int j = 0;
  int k = 0;
  double W_pot = 0;
  double *pair_potentials, *ro;

  /* TODO: Alloc once at creation */
  pair_potentials = (double *)malloc(atoms_cnt*sizeof(double));
  ro = (double *)malloc(atoms_cnt*sizeof(double));

  memset(pair_potentials, 0, atoms_cnt*sizeof(double));
  memset(ro, 0, atoms_cnt*sizeof(double));

  for (i=0; i < atoms_cnt; i++)
  {
    for (j=i+1; j < atoms_cnt; j++)
    {
      double dist = sqrt(SQR(x[i*COORDS_CNT] - x[j*COORDS_CNT]) + 
                        SQR(x[i*COORDS_CNT+1] - x[j*COORDS_CNT+1]) +
                        SQR(x[i*COORDS_CNT+2] - x[j*COORDS_CNT+2]));
      double a_dist = A_PARAM / dist;
      double dist_pow_n = pow(a_dist, N_PARAM);
      double dist_pow_m = pow(a_dist, M_PARAM);
      pair_potentials[i] += dist_pow_n;
      pair_potentials[j] += dist_pow_n;
      ro[i] += dist_pow_m;
      ro[j] += dist_pow_m;
    }

    W_pot +=  0.5*pair_potentials[i] - C_PARAM*sqrt(ro[i]);
  }

  f[0] = EPS_PARAM * W_pot;

  free(pair_potentials);
  free(ro);
}


void cga::objfun_impl_test(pagmo::fitness_vector &f, const pagmo::decision_vector &x) const
{
  this->objfun_impl(f,x);
}

std::string cga::get_name() const
{
  return "cga";
}

BOOST_CLASS_EXPORT_IMPLEMENT(cga);
