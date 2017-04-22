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

#include <queso/GammaVectorRV.h>
#include <queso/GammaVectorRealizer.h>
#include <queso/GammaJointPdf.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Constructor---------------------------------------
// TODO: Check: the constructor receives imageSet, but uses m_imageSet to assign m_pdf and m_realizer. (Kemelli 2013/4/22)
template<class V, class M>
GammaVectorRV<V,M>::GammaVectorRV(
  const char*                  prefix,
  const VectorSet<V,M>& imageSet,
  const V&                     a,
  const V&                     b)
  :
  BaseVectorRV<V,M>(((std::string)(prefix)+"uni").c_str(),imageSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering GammaVectorRV<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

// begin kemelli 2013-April-22 --------------------------
// better to check for the parameter values in the constructor,
// rather than in GammaVectorRealizer<V,M>::realization

  const BoxSubset<V,M>* imageBox = dynamic_cast<const BoxSubset<V,M>* >(&imageSet);
  double smallerOfMaxValues = imageBox->maxValues().getMinValue();
  double smallerOfMinValues = imageBox->minValues().getMinValue();

 // Gamma dist is defined only in (0,inf)
 if( smallerOfMinValues < 0 )
 {
   std::cerr << "In GammaVectorRV<V,M>::constructor()\n"
       << "Gamma distribution is only defined in (0, infinity).\n"
       << "The data provided is: \n"
       << *imageBox
         << "Sampling will not cover all interval.\n"
         << std::endl;

    queso_require_greater_equal_msg(smallerOfMaxValues, 0, "invalid input: Gamma distribution is only defined in (0, infinity), and min(m_maxValues)<0");
 }
// end kemelli 2013-April-22 --------------------------

  m_pdf        = new GammaJointPdf<V,M>(m_prefix.c_str(),
                                               m_imageSet,
                                               a,
                                               b);
  m_realizer   = new GammaVectorRealizer<V,M>(m_prefix.c_str(),
                                                     m_imageSet,
                                                     a,
                                                     b);

  m_subCdf     = NULL; // FIX ME: complete code
  m_unifiedCdf = NULL; // FIX ME: complete code
  m_mdf        = NULL; // FIX ME: complete code

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving GammaVectorRV<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V, class M>
GammaVectorRV<V,M>::~GammaVectorRV()
{
  delete m_mdf;
  delete m_unifiedCdf;
  delete m_subCdf;
  delete m_realizer;
  delete m_pdf;
}
// I/O methods---------------------------------------
template <class V, class M>
void
GammaVectorRV<V,M>::print(std::ostream& os) const
{
  os << "GammaVectorRV<V,M>::print() says, 'Please implement me.'" << std::endl;
  return;
}

}  // End namespace QUESO

template class QUESO::GammaVectorRV<QUESO::GslVector,QUESO::GslMatrix>;
