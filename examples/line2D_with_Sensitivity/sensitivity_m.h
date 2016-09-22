#ifndef QUESO_EXAMPLE_SENSITIVITY_M_H
#define QUESO_EXAMPLE_SENSITIVITY_M_H

#include <queso/VectorFunction.h>
#include <queso/DistArray.h>

template<class P_V = QUESO::GslVector, class P_M = QUESO::GslMatrix,
	 class Q_V = QUESO::GslVector, class Q_M = QUESO::GslMatrix>
class Qoi_m : public QUESO::BaseVectorFunction<P_V, P_M, Q_V, Q_M>
{
public:
	Qoi_m(const char * prefix, const QUESO::VectorSet<P_V, P_M> & domainSet,
	    const QUESO::VectorSet<Q_V, Q_M> & imageSet);
	virtual ~Qoi_m();
	virtual void compute(const P_V & domainVector, const P_V * domainDirection,
		Q_V & imageVector, QUESO::DistArray<P_V *> * gradVectors,
		QUESO::DistArray<P_M *> * hessianMatrices,
		QUESO::DistArray<P_V *> * hessianEffects) const;

private:
	double x_loc;

};

#endif
	
