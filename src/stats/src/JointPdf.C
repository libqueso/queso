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

#include <queso/JointPdf.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Default constructor -----------------------------
template<class V, class M>
BaseJointPdf<V,M>::BaseJointPdf(
  const char*                  prefix,
  const VectorSet<V,M>& domainSet)
  :
  BaseScalarFunction<V,M>(((std::string)(prefix)+"pd_").c_str(), domainSet),
  m_normalizationStyle(0),
  m_logOfNormalizationFactor(0.)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering BaseJointPdf<V,M>::constructor() [3]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving BaseJointPdf<V,M>::constructor() [3]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V, class M>
BaseJointPdf<V,M>::~BaseJointPdf()
{
}

//---------------------------------------------------
template<class V,class M>
void
BaseJointPdf<V,M>::distributionMean(V & meanVector) const
{
  queso_not_implemented();
}

//---------------------------------------------------
template<class V,class M>
void
BaseJointPdf<V,M>::distributionVariance(M & covMatrix) const
{
  queso_not_implemented();
}

// Math methods -------------------------------------
template<class V,class M>
void
BaseJointPdf<V,M>::setNormalizationStyle(unsigned int value) const
{
  m_normalizationStyle = value;
  return;
}
//---------------------------------------------------
template<class V,class M>
void
BaseJointPdf<V,M>::setLogOfNormalizationFactor(double value) const
{
  m_logOfNormalizationFactor = value;
  return;
}

template <class V, class M>
void
BaseJointPdf<V, M>::print(std::ostream & os) const
{
  // Print m_env?
  // Print mean?
  // Print var?

  os << "Start printing BaseJointPdf<V, M>" << std::endl;
  os << "m_prefix:" << std::endl;
  os << this->m_prefix << std::endl;
  os << "m_domainSet:" << std::endl;
  os << this->m_domainSet << std::endl;
  os << "m_normalizationStyle:" << std::endl;
  os << this->m_normalizationStyle << std::endl;
  os << "m_logOfNormalizationFactor:" << std::endl;
  os << this->m_logOfNormalizationFactor << std::endl;
  os << "End printing BaseJointPdf<V, M>" << std::endl;
}

//---------------------------------------------------
template<class V,class M>
double
BaseJointPdf<V,M>::commonComputeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const
{
  double value = 0.;

  double volume = m_domainSet.volume();
  if ((queso_isnan(volume)) ||
      (volume == -INFINITY         ) ||
      (volume ==  INFINITY         ) ||
      (volume <= 0.                )) {
    // Do nothing
  }
  else {
    const BoxSubset<V,M>* boxSubset = dynamic_cast<const BoxSubset<V,M>* >(&m_domainSet);
    if (boxSubset == NULL) {
      // Do nothing
    }
    else {
      V tmpVec(m_domainSet.vectorSpace().zeroVector());
      double sum = 0.;
      for (unsigned int i = 0; i < numSamples; ++i) {
        tmpVec.cwSetUniform(boxSubset->minValues(),boxSubset->maxValues());
        sum += this->actualValue(tmpVec,NULL,NULL,NULL,NULL);
      }
      double avgValue = sum/((double) numSamples);
      value = -( log(avgValue) + log(volume) );
      if (updateFactorInternally) {
        m_logOfNormalizationFactor = value;
      }
    }
  }

  return value;
}

#if 0
template <class V, class M>
const BaseScalarPdf&
BaseJointPdf<V,M>::component(unsigned int componentId) const
{
  if (componentId > m_components.size()) return m_dummyComponent;
  if (m_components[componentId] == NULL) return m_dummyComponent;
  return *(m_components[componentId]);
}
#endif

}  // End namespace QUESO

template class QUESO::BaseJointPdf<QUESO::GslVector, QUESO::GslMatrix>;
