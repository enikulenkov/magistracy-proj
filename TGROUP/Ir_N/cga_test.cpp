#include <stdio.h>
#include <stdlib.h>
#include <pagmo/pagmo.h>
#include <gtest/gtest.h>
#include "ga_utils.h"
#include "cga.h"

double Clength=0.529177;   //http://en.wikipedia.org/wiki/Bohr_radius
double g_a0 = 3.8344/Clength; /* Lattice constant. http://en.wikipedia.org/wiki/Lattice_constant */


TEST(PotentialEnergy, SimpleTest3)
{
  const int atoms_cnt = 50;
  char test_file[]="50res.mat";
  double *coords;
  int i;
  pagmo::fitness_vector f = pagmo::fitness_vector(1);
  pagmo::decision_vector x = pagmo::decision_vector(atoms_cnt*3);

  read_input(test_file, &coords);

  for (i = 0; i < x.size(); i++)
  {
    x[i] = coords[i];
  }

  cga prob = cga(atoms_cnt);
  prob.objfun_impl_test(f, x);

  EXPECT_DOUBLE_EQ(-10.380928756912276, f[0]);
}


TEST(PotentialEnergy, SimpleTest50)
{
  const int atoms_cnt = 50;
  char test_file[]="FCC.mat";
  double *coords;
  int i;
  pagmo::fitness_vector f = pagmo::fitness_vector(1);
  pagmo::decision_vector x = pagmo::decision_vector(atoms_cnt*3);

  read_input(test_file, &coords);

  for (i = 0; i < x.size(); i++)
  {
    x[i] = g_a0 * coords[i];
  }

  cga prob = cga(atoms_cnt);
  prob.objfun_impl_test(f, x);

  EXPECT_DOUBLE_EQ(-10.299435778260022, f[0]);
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
