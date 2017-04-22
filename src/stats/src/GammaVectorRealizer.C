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

#include <queso/GammaVectorRealizer.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Constructor -------------------------------------
template<class V, class M>
GammaVectorRealizer<V,M>::GammaVectorRealizer(
  const char*                  prefix,
  const VectorSet<V,M>& unifiedImageSet,
  const V&                     a,
  const V&                     b)
  :
  BaseVectorRealizer<V,M>(((std::string)(prefix)+"gen").c_str(),unifiedImageSet,std::numeric_limits<unsigned int>::max()),
  m_a(a),
  m_b(b)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering GammaVectorRealizer<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving GammaVectorRealizer<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor --------------------------------------
template<class V, class M>
GammaVectorRealizer<V,M>::~GammaVectorRealizer()
{
}
// Realization-related methods----------------------
template<class V, class M>
void
GammaVectorRealizer<V,M>::realization(V& nextValues) const
{
  // nextValues.cwSetGamma(m_a,m_b);

// begin kemelli 2013-April-22 :
  const BoxSubset<V,M>* imageBox = dynamic_cast<const BoxSubset<V,M>* >(&this->m_unifiedImageSet);
//  double biggerOfMaxValues = imageBox->maxValues().getMaxValue();
  double smallerOfMaxValues = imageBox->maxValues().getMinValue();
  double smallerOfMinValues = imageBox->minValues().getMinValue();

 // Gamma dist belongs to (0,inf)
 if( smallerOfMinValues < 0 ) //(biggerOfMinValues < 0) ||
 {
   std::cerr << "In GammaVectorRealizer<V,M>::realization()\n"
       << "Gamma distribution is only defined in (0, infinity).\n"
       << "The data provided is: \n"
       << *imageBox
         << "Sampling will not cover all interval.\n"
         << std::endl;


    queso_require_greater_equal_msg(smallerOfMaxValues, 0, "invalid input: Gamma distribution is only defined in (0, infinity), and min(m_maxValues)<0. ");


 //                     m_env.worldRank(),
 //                     "GammaVectorRealizer<V,M>::realization()",
 //                     "invalid input: Gamma distribution is only defined in (0, infinity). ");

 }

  // end kemelli 2013-April-22

  // begin kemelli 2013-April-18 :
  bool outOfSupport = true;
  do {

  nextValues.cwSetGamma(m_a,m_b);
  outOfSupport = !(this->m_unifiedImageSet.contains(nextValues));

  } while (outOfSupport);

  // end kemelli 2013-April-18 :

  return;
}

}  // End namespace QUESO

template class QUESO::GammaVectorRealizer<QUESO::GslVector, QUESO::GslMatrix>;
