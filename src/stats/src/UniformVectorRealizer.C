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
#include <queso/math_macros.h>
#include <queso/UniformVectorRealizer.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Constructor -------------------------------------
template<class V, class M>
UniformVectorRealizer<V,M>::UniformVectorRealizer(
  const char*                  prefix,
  const VectorSet<V,M>& unifiedImageSet)
  :
  BaseVectorRealizer<V,M>(((std::string)(prefix)+"gen").c_str(),unifiedImageSet,std::numeric_limits<unsigned int>::max())
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering UniformVectorRealizer<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving UniformVectorRealizer<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor --------------------------------------
template<class V, class M>
UniformVectorRealizer<V,M>::~UniformVectorRealizer()
{
}
// Realization-related methods----------------------
template<class V, class M>
void
UniformVectorRealizer<V,M>::realization(V& nextValues) const
{
  const BoxSubset<V,M>* imageBox = dynamic_cast<const BoxSubset<V,M>* >(&m_unifiedImageSet);

  queso_require_msg(imageBox, "only box images are supported right now");
  queso_require_msg(queso_isfinite(imageBox->volume()), "drawing realisations from an improper uniform is not supported");

  nextValues.cwSetUniform(imageBox->minValues(),imageBox->maxValues());
}

}  // End namespace QUESO

template class QUESO::UniformVectorRealizer<QUESO::GslVector, QUESO::GslMatrix>;
