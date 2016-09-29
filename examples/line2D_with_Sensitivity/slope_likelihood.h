#ifndef SLOPE_LIKELIHOOD_H
#define SLOPE_LIKELIHOOD_H

#include <queso/ScalarFunction.h>
#include <queso/GslMatrix.h>

template<class V = QUESO::GslVector, class M = QUESO::GslMatrix>
class Likelihood : public QUESO::BaseScalarFunction<V, M>
{
public:
  Likelihood(const char * prefix, const QUESO::VectorSet<V, M> & domain);
  virtual ~Likelihood();
  virtual double lnValue(const V &domainVector, const V * domainDirection,
	V * gradVector, M * hessianMatrix, V * hessianEffect) const;
  virtual double actualValue(const V & domainVector, const V * domainDirection,
      V * gradVector, M * hessianMatrix, V * hessianEffect) const;  

private:
  std::vector<double> x_vals; // point along the x axes
  std::vector<double> y_obs; // noisy y estimates
  std::vector<double> stdDevs; // std Dev in y_obs
 };

#endif
