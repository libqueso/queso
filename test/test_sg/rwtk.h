#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/TKGroup.h>
#include <queso/SharedPtr.h>
#include <queso/ScopedPtr.h>

namespace QUESO
{
  template <class V, class M> class VectorSpace;
  template <class V, class M> class BaseVectorRV;
}

template <class V = QUESO::GslVector, class M = QUESO::GslMatrix>
class RWTK : public QUESO::BaseTKGroup<V, M>
{
public:
  RWTK(const char * prefix, const QUESO::VectorSpace<V, M> & vectorSpace,
       const std::vector<double> & scales);

  virtual ~RWTK();

  virtual bool symmetric() const;
  virtual const QUESO::BaseVectorRV<V, M> & rv(unsigned int stageId) const;
  virtual const QUESO::BaseVectorRV<V, M> & rv(
      const std::vector<unsigned int> & stageIds);

private:
  //! Disallow calling the default constructor
  /*!
   * We do this because it does not contain any QUESO environemtn information
   */
  RWTK();

  //! The underlying random variable object
  typename QUESO::SharedPtr<QUESO::BaseVectorRV<V, M> >::Type m_rv;
  typename QUESO::ScopedPtr<QUESO::GslVector>::Type m_mean;
  typename QUESO::ScopedPtr<QUESO::GslMatrix>::Type m_covariance;
};
