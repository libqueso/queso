#include <cmath>

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <sensitivity_m.h>

template<class P_V, class P_M, class Q_V, class Q_M>
Qoi_m<P_V, P_M, Q_V, Q_M>::Qoi_m(const char * prefix,
    const QUESO::VectorSet<P_V, P_M> & domainSet,
    const QUESO::VectorSet<Q_V, Q_M> & imageSet)
  : QUESO::BaseVectorFunction<P_V, P_M, Q_V, Q_M>(prefix, domainSet, imageSet),
    x_loc(3)
{
}

template<class P_V, class P_M, class Q_V, class Q_M>
Qoi_m<P_V, P_M, Q_V, Q_M>::~Qoi_m()
{
  // Deconstruct here
}

template<class P_V, class P_M, class Q_V, class Q_M>
void
Qoi_m<P_V, P_M, Q_V, Q_M>::compute(const P_V & domainVector,
    const P_V * /* domainDirection */,
    Q_V & imageVector, QUESO::DistArray<P_V *> * /* gradVectors */,
    QUESO::DistArray<P_M *> * /* hessianMatrices */,
    QUESO::DistArray<P_V *> * /* hessianEffects */) const
{
  if (domainVector.sizeLocal() != 2) {
    queso_error_msg("domainVector does not have size 2");
  }
  if (imageVector.sizeLocal() != 1) {
    queso_error_msg("imageVector does not have size 1");
  }

//  std::cout << "m = " << domainVector[0] << std::endl; 
  double m = domainVector[0];  // Sample of the RV 'line slope'
  double c = 5.0;  // Sample of the RV 'y-intercept'
  double y_obs = 0.0;
  y_obs = m*x_loc + c;

  imageVector[0] = y_obs;
}

template class Qoi_m<QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector,
                   QUESO::GslMatrix>;
