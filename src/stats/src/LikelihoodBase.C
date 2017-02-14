//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2017 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

#include <sstream>

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorSet.h>
#include <queso/LikelihoodBase.h>
#include <queso/VectorRV.h>
#include <queso/MultiDQuadratureBase.h>

namespace QUESO {

template<class V, class M>
LikelihoodBase<V, M>::LikelihoodBase(
    const char * prefix, const VectorSet<V, M> & domainSet,
    const V & observations)
  : BaseScalarFunction<V, M>(prefix, domainSet),
    m_observations(observations)
{
}

template<class V, class M>
LikelihoodBase<V, M>::LikelihoodBase(const char * prefix,
                                     const VectorSet<V, M> & domainSet,
                                     const V & observations,
                                     typename SharedPtr<BaseVectorRV<V,M> >::Type & marg_param_pdf,
                                     typename SharedPtr<MultiDQuadratureBase<V,M> >::Type & marg_integration,
                                     bool marg_pdf_is_weight_func)
  : BaseScalarFunction<V, M>(prefix, domainSet),
  m_observations(observations),
  m_marg_param_pdf(marg_param_pdf),
  m_marg_integration(marg_integration),
  m_marg_pdf_is_weight_func(marg_pdf_is_weight_func)
{
  // The dimension of the parameter space had better match the dimension in the integration
  queso_require_equal_to_msg(this->m_marg_param_pdf->imageSet().vectorSpace().dimGlobal(),
                             this->m_marg_integration->getDomain().vectorSpace().dimGlobal(),
                             "Mismatched marginal parameter space dimension and quadrature dimension!");
}

template<class V, class M>
LikelihoodBase<V, M>::~LikelihoodBase()
{
}

template<class V, class M>
void LikelihoodBase<V, M>::evaluateModel(const V & domainVector, const V * domainDirection,
                                         V & modelOutput, V * gradVector, M * hessianMatrix,
                                         V * hessianEffect) const
{
  std::stringstream ss;
  ss << "ERROR: evaluateModel() not implemented! This interface is deprecated."
     << std::endl
     << "Prefer implementing evaluateModel(const V & domainVector, V & modelOutput)"
     << std::endl;

  queso_error_msg(ss.str());
}

template<class V, class M>
void LikelihoodBase<V, M>::evaluateModel(const V & domainVector,
                                         const V & marginalVector,
                                         V & modelOutput) const
{
  std::stringstream ss;
  ss << "ERROR: evaluateModel(const V & domainVector, const V & marginalVector, V & modelOutput)"
     << std::endl
     << "       is not implemented! Please override this function in your subclass."
     << std::endl;

  queso_error_msg(ss.str());
}

template<class V, class M>
double
LikelihoodBase<V, M>::lnValue(const V & domainVector) const
{
  V modelOutput(this->m_observations, 0, 0);

  double lnLikelihood_value = 0.0;

  // If we're not marginalizing, then we just need f(m)
  if( !m_marg_param_pdf )
    {
      this->evaluateModel(domainVector, modelOutput);

      // Note modelOutput made be modified in lnLikelihood()
      lnLikelihood_value = this->lnLikelihood(domainVector, modelOutput);
    }
  // Otherwise we're integrating over marginal parameter space
  else
    {
      queso_assert(m_marg_integration);

      const std::vector< typename QUESO::SharedPtr< V >::Type > & x =
        this->m_marg_integration->positions();

      const std::vector< double > & w = this->m_marg_integration->weights();

      unsigned int n_qpoints = x.size();

      for( unsigned int q = 0; q < n_qpoints; q++ )
        {
          this->evaluateModel(domainVector, *(x[q]), modelOutput);

          // Note modelOutput made be modified in lnLikelihood()
          double ln_pi_m_q = this->lnLikelihood(domainVector, modelOutput);

          double ln_pi_q = 0.0;

          // If we're not treating the marginal parameters pdf as
          // the weighting function in the quadrature evaluation,
          // then we need to also evaluate it.
          if(!m_marg_pdf_is_weight_func)
            this->m_marg_param_pdf->pdf().lnValue( *(x[q]) );

          /*! \todo [PB]: We might want to play games with a log(\sum) identity
                    if precision starts to become an issue. */
          lnLikelihood_value += std::exp(ln_pi_m_q + ln_pi_q)*w[q];
        }

      lnLikelihood_value = std::log(lnLikelihood_value);
    }

  return lnLikelihood_value;
}

}  // End namespace QUESO

template class QUESO::LikelihoodBase<QUESO::GslVector, QUESO::GslMatrix>;
