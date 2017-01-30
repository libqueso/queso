#include <cmath>
#include <fstream>
#include <iostream>

#include <queso/GenericScalarFunction.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/UniformVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/ScalarFunction.h>
#include <queso/VectorSet.h>

#include <slope_likelihood.h>

template<class V, class M>
Likelihood<V, M>::Likelihood(const char * prefix,
	const QUESO::VectorSet<V, M> & domain)
      : QUESO::BaseScalarFunction<V, M>(prefix, domain),
	x_vals(0),
	y_obs(0),
	stdDevs(0)
{
  double xs[21], ys[21], sigmas[21];
  std::ifstream xs_file("xs.txt");
  std::ifstream ys_file("obs.txt");
  std::ifstream sig_file("sigmas.txt");

  //int i = 0;
  //std::string line[30];

  for (int i = 0; i < 21; i++) {
   
	xs_file >> xs[i];
	ys_file >> ys[i];
	sig_file >> sigmas[i];
	//std::cout << xs[i] << " " << ys[i] << " " << sigmas[i] << std::endl;
  }
  
  std::size_t const n = sizeof(xs) / sizeof(*xs);
  x_vals.assign(xs, xs + n);
  y_obs.assign(ys, ys + n);
  stdDevs.assign(sigmas, sigmas + n);
}

template<class V, class M>
Likelihood<V, M>::~Likelihood()
{
  // Deconstruct here
}


template<class V, class M>
double
Likelihood<V, M>::lnValue(const V & domainVector) const
{
  double m = domainVector[0];
  double c = domainVector[1];
  double misfit = 0.0;
  double model_y, ratio;
  
  for (unsigned int i = 0; i < x_vals.size(); i++){
	model_y = m*x_vals[i] + c;
	ratio = (model_y - y_obs[i]) / stdDevs[i];
	misfit += ratio*ratio;
  }

  return -0.5 * misfit;
}

template<class V, class M>
double
Likelihood<V, M>::actualValue(const V & domainVector,
    const V * domainDirection, V * gradVector, M * hessianMatrix,
    V * hessianEffect) const
{
  return std::exp(this->lnValue(domainVector));
}

template class Likelihood<QUESO::GslVector, QUESO::GslMatrix>;
 
