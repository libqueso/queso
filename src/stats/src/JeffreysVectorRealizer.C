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

#include <limits>
#include <queso/JeffreysVectorRealizer.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Constructor -------------------------------------
template<class V, class M>
JeffreysVectorRealizer<V,M>::JeffreysVectorRealizer(
  const char*                  prefix,
  const VectorSet<V,M>& unifiedImageSet)
  :
  BaseVectorRealizer<V,M>(((std::string)(prefix)+"gen").c_str(),unifiedImageSet,std::numeric_limits<unsigned int>::max())
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering JeffreysVectorRealizer<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving JeffreysVectorRealizer<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor --------------------------------------
template<class V, class M>
JeffreysVectorRealizer<V,M>::~JeffreysVectorRealizer()
{
}
// Realization-related methods----------------------
template<class V, class M>
void
JeffreysVectorRealizer<V,M>::realization(V& nextValues) const
{
  const BoxSubset<V,M>* imageBox = dynamic_cast<const BoxSubset<V,M>* >(&m_unifiedImageSet);

  if (imageBox == NULL) {
    queso_error_msg("For JeffreysVectorRealizer<V,M>::realization(), only box images are supported right now");
  }
  //take log of Jeffreys bounds to set uniform bounds
  GslVector logMinValues(imageBox->minValues());
  for (unsigned int i = 0; i < logMinValues.sizeLocal(); ++i) {
    if (logMinValues[i] < 0.0) {
      queso_error_msg("The minimum value for a Jeffreys distribution should be greater than or equal to zero.");
    }
    else {
      logMinValues[i] = std::log(logMinValues[i]);
    }
  }

  GslVector logMaxValues(imageBox->maxValues());
  for (unsigned int i = 0; i < logMaxValues.sizeLocal(); ++i) {
    if (logMaxValues[i] <= 0.0) {
      queso_error_msg("The maximum value for a Jeffreys distribution should be greater than zero.");
    }
    else {
      logMaxValues[i] = std::log(logMaxValues[i]);
    }
  }

  nextValues.cwSetUniform(logMinValues,logMaxValues);
  for (unsigned int i = 0; i < nextValues.sizeLocal(); ++i) {
    nextValues[i] = std::exp(nextValues[i]);
  }
  return;
}

}  // End namespace QUESO

template class QUESO::JeffreysVectorRealizer<QUESO::GslVector, QUESO::GslMatrix>;
