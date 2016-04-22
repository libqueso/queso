#include "rwtk.h"

#include <queso/VectorSpace.h>
#include <queso/VectorRV.h>

template <class V, class M>
RWTK<V, M>::RWTK(const char * prefix,
                 const QUESO::VectorSpace<V, M> & vectorSpace,
                 const std::vector<double> & scales)
: QUESO::BaseTKGroup<V, M>(prefix, vectorSpace, scales),
  m_rv(),
  m_mean(),
  m_covariance()
{
  m_mean.reset(new QUESO::GslVector(vectorSpace.zeroVector()));
  m_mean->cwSet(0.0);

  m_covariance.reset(new QUESO::GslMatrix(vectorSpace.zeroVector()));
  (*m_covariance)(0, 0) = 1.0;
  (*m_covariance)(0, 1) = 0.0;
  (*m_covariance)(1, 0) = 0.0;
  (*m_covariance)(1, 1) = 1.0;

  m_rv.reset(new QUESO::GaussianVectorRV<V, M>("", vectorSpace, *m_mean, *m_covariance));
}

template <class V, class M>
RWTK<V, M>::~RWTK()
{
  // Do nothing?
}

template <class V, class M>
bool
RWTK<V, M>::symmetric() const
{
  return true;
}

template <class V, class M>
const QUESO::BaseVectorRV<V, M> &
RWTK<V, M>::rv(unsigned int stageId) const
{
  return *(this->m_rv);
}

template <class V, class M>
const QUESO::BaseVectorRV<V, M> &
RWTK<V, M>::rv(const std::vector<unsigned int> & stageIds)
{
  return *(this->m_rv);
}

template class RWTK<QUESO::GslVector, QUESO::GslMatrix>;
