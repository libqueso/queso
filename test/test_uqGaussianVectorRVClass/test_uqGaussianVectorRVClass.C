#include <stdio.h>
#include <string.h>
#include <cmath>
#include <uqEnvironment.h>
#include <uqGslVector.h>
#include <uqGslMatrix.h>
#include <uqVectorSpace.h>
#include <uqVectorSubset.h>
#include <uqVectorRV.h>

using namespace std;

int main(int argc, char **argv) {
  int i, j, num_samples = 1e5, seed;
  double sigmasq = 1.0;
  double *mean;
  double *var;
  double *sumsq;
  double *delta;
  FILE *fpRand = fopen("/dev/random", "r");

  fread(&seed, sizeof(int), 1, fpRand);
  fclose(fpRand);

  uqEnvOptionsValuesClass *opts = new uqEnvOptionsValuesClass();
  opts->m_seed = seed;

#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
#endif
  uqFullEnvironmentClass *env =
#ifdef QUESO_HAS_MPI
    new uqFullEnvironmentClass(MPI_COMM_WORLD, NULL, "", opts);
#else
    new uqFullEnvironmentClass(0, NULL, "", opts);
#endif

  uqVectorSpaceClass<uqGslVectorClass, uqGslMatrixClass> *param_space;
  param_space = new uqVectorSpaceClass<uqGslVectorClass, uqGslMatrixClass>(
      *env, "param_", 4, NULL);

  mean = new double[4];
  var = new double[4];
  sumsq = new double[4];
  delta = new double[4];

  uqGslVectorClass mins(param_space->zeroVector());
  uqGslVectorClass maxs(param_space->zeroVector());
  mins.cwSet(-INFINITY);
  maxs.cwSet(INFINITY);
  uqBoxSubsetClass<uqGslVectorClass, uqGslMatrixClass> *param_domain;
  param_domain = new uqBoxSubsetClass<uqGslVectorClass, uqGslMatrixClass>(
      "param_", *param_space, mins, maxs);

  // Mean zero
  uqGslVectorClass prior_mean_vals(param_space->zeroVector());
  // Variance
  uqGslVectorClass prior_var_vals(param_space->zeroVector());
  prior_var_vals.cwSet(sigmasq);

  uqGaussianVectorRVClass<uqGslVectorClass, uqGslMatrixClass> *prior;
  prior = new uqGaussianVectorRVClass<uqGslVectorClass, uqGslMatrixClass>(
      "prior_", *(param_domain), prior_mean_vals, prior_var_vals);

  uqGslVectorClass draw(param_space->zeroVector());

  for (j = 0; j < 4; j++) {
    mean[j] = 0.0;
    sumsq[j] = 0.0;
  }

  for (i = 1; i < num_samples + 1; i++) {
    prior->realizer().realization(draw);
    for (j = 0; j < 4; j++) {
      delta[j] = draw[j] - mean[j];
      mean[j] += (double) delta[j] / i;
      sumsq[j] += delta[j] * (draw[j] - mean[j]);
    }
  }

  for (j = 0; j < 4; j++) {
    // This is the sample variance
    var[j] = sumsq[j] / (num_samples - 1);
  }

  double mean_min = -3.0 * sqrt(sigmasq) / sqrt(num_samples);
  double mean_max =  3.0 * sqrt(sigmasq) / sqrt(num_samples);

  int return_val = 0;
  for (j = 0; j < 4; j++) {
    if (mean[j] < mean_min || mean[j] > mean_max) {
      return_val = 1;
      break;
    }
  }

  double var_min = sigmasq - 3.0 * sqrt(2.0 * sigmasq * sigmasq / (num_samples - 1));
  double var_max = sigmasq + 3.0 * sqrt(2.0 * sigmasq * sigmasq / (num_samples - 1));

  // var[j] should be approximately ~ N(sigma^2, 2 sigma^4 / (num_samples - 1))
  for (j = 0; j < 4; j++) {
    if (var[j] < var_min || var[j] > var_max) {
      return_val = 1;
      break;
    }
  }

  delete env;
  delete opts;
  delete mean;
  delete var;
  delete delta;
  delete sumsq;
  delete param_space;
  delete param_domain;
  delete prior;
#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return return_val;
}
