/*
 * Cluster genetic algorithm
 */
#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <string>

#include <pagmo/pagmo.h>
//#include <pagmo.h>
#include "sutton_chen_pot.h"

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

sutton_chen_pot::sutton_chen_pot(int n):pagmo::problem::base(n*3)
{
  /* nothing */
}


/* Clone method */
pagmo::problem::base_ptr sutton_chen_pot::clone() const
{
  return pagmo::problem::base_ptr(new sutton_chen_pot(*this));
}


void sutton_chen_pot::d_objfun(const double *vec, int vec_size, double *df) const
{
  const int atoms_cnt = vec_size/COORDS_CNT;
  double *ro;
  double *fp;
  double *sp;
  double ro_sum = 0;
  ro = (double *)malloc(atoms_cnt*sizeof(double));
  fp = (double *)malloc(vec_size * 3 * sizeof(double));
  sp = (double *)malloc(vec_size * 3 * sizeof(double));

  memset(ro, 0, atoms_cnt*sizeof(double));
  memset(fp, 0, vec_size*3*sizeof(double));
  memset(sp, 0, vec_size*3*sizeof(double));

  /* Calculate sums of ro */
  for (int k = 0; k < atoms_cnt; k++)
  {
    for (int j = k+1; j < atoms_cnt; j++)
    {
      double dist = sqrt(SQR(vec[k*COORDS_CNT] - vec[j*COORDS_CNT]) +
                        SQR(vec[k*COORDS_CNT+1] - vec[j*COORDS_CNT+1]) +
                        SQR(vec[k*COORDS_CNT+2] - vec[j*COORDS_CNT+2]));
      double a_dist = A_PARAM / dist;
      //double dist_pow_n = pow(a_dist, N_PARAM);
      double dist_pow_m = pow(a_dist, M_PARAM);
      //pair_potentials[k] += dist_pow_n;
      //pair_potentials[j] += dist_pow_n;
      ro[k] += dist_pow_m;
      ro[j] += dist_pow_m;
      //dist_squared[k] += dist*dist;
      //dist_squared[j] += dist*dist;
    }
  }

  for (int k = 0; k < atoms_cnt; k++)
  {
    for (int j = k+1; j < atoms_cnt; j++)
    {
      double dist = sqrt(SQR(vec[k*COORDS_CNT] - vec[j*COORDS_CNT]) +
                        SQR(vec[k*COORDS_CNT+1] - vec[j*COORDS_CNT+1]) +
                        SQR(vec[k*COORDS_CNT+2] - vec[j*COORDS_CNT+2]));
      double a_dist = A_PARAM / dist;
      double dist_pow_n = pow(a_dist, N_PARAM);
      double dist_pow_m = pow(a_dist, M_PARAM);

      for (int i = 0; i < COORDS_CNT; i++)
      {
        double diff =  vec[j*COORDS_CNT+i] - vec[k*COORDS_CNT+i];
        double tmp1 = diff/(dist*dist) * dist_pow_n;
        double tmp2 = diff/(dist*dist) * (pow(ro[k], -0.5) + pow(ro[j], -0.5)) * dist_pow_m;
        fp[k*COORDS_CNT+i] += tmp1;
        fp[j*COORDS_CNT+i] += -tmp1;
        sp[k*COORDS_CNT+i] += tmp2;
        sp[j*COORDS_CNT+i] += -tmp2;
      }
    }
  }

  for (int k = 0; k < atoms_cnt ; k++)
  {
    for (int i = 0; i < COORDS_CNT; i++)
    {
      df[k*COORDS_CNT+i] = (-1)*EPS_PARAM*((-1)*(N_PARAM)*fp[k*COORDS_CNT+i] + C_PARAM*M_PARAM*0.5*sp[k*COORDS_CNT+i]);
    }
  }

  free(ro);
  free(fp);
  free(sp);
}


void sutton_chen_pot::objfun_impl(pagmo::fitness_vector &f, const pagmo::decision_vector &x) const
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


void sutton_chen_pot::objfun_impl_test(pagmo::fitness_vector &f, const pagmo::decision_vector &x) const
{
  this->objfun_impl(f,x);
}

std::string sutton_chen_pot::get_name() const
{
  return "Sutton-Chen potential";
}

BOOST_CLASS_EXPORT_IMPLEMENT(sutton_chen_pot);
