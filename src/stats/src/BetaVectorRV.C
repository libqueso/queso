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

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

#include <queso/BetaJointPdf.h>
#include <queso/BetaVectorRV.h>
#include <queso/BetaVectorRealizer.h>

namespace QUESO {

// Constructor---------------------------------------
template<class V, class M>
BetaVectorRV<V,M>::BetaVectorRV(
  const char*                  prefix,
  const VectorSet<V,M>& imageSet,
  const V&                     alpha,
  const V&                     beta)
  :
  BaseVectorRV<V,M>(((std::string)(prefix)+"uni").c_str(),imageSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering BetaVectorRV<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

// begin kemelli 2013-April-22 : --------------------------

  const BoxSubset<V,M>* imageBox = dynamic_cast<const BoxSubset<V,M>* >(&imageSet);

  double smallerOfMaxValues = imageBox->maxValues().getMinValue();
  double biggerOfMaxValues = imageBox->maxValues().getMaxValue();
  double smallerOfMinValues = imageBox->minValues().getMinValue();
  double biggerOfMinValues = imageBox->minValues().getMaxValue();

 // Beta dist is defined only in [0,1]
 if( (smallerOfMinValues < 0) || ( biggerOfMaxValues > 1 ) )
 {
   std::cerr << "In BetaVectorRV<V,M>::constructor()\n"
       << "Beta distribution is defined only in [0, 1].\n"
       << "The data provided is: \n"
       << *imageBox
         << "Sampling will not cover all interval.\n"
         << std::endl;

 // if at least one of the min values > 1 then exit
    queso_require_less_equal_msg(biggerOfMinValues, 1, "invalid input: Beta distribution is only defined in [0, 1], and max(m_minValues)>1");
 // if at least one of the max values < 0 then exit
    queso_require_greater_equal_msg(smallerOfMaxValues, 0, "invalid input: Beta distribution is only defined in [0, 1], and min(m_maxValues)<0");
 }
  // end kemelli 2013-April-22 --------------------------

  m_pdf        = new BetaJointPdf<V,M>(m_prefix.c_str(),
                                              m_imageSet,
                                              alpha,
                                              beta);
  m_realizer   = new BetaVectorRealizer<V,M>(m_prefix.c_str(),
                                                    m_imageSet,
                                                    alpha,
                                                    beta);
  m_subCdf     = NULL; // FIX ME: complete code
  m_unifiedCdf = NULL; // FIX ME: complete code
  m_mdf        = NULL; // FIX ME: complete code

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving BetaVectorRV<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V, class M>
BetaVectorRV<V,M>::~BetaVectorRV()
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
BetaVectorRV<V,M>::print(std::ostream& os) const
{
  os << "BetaVectorRV<V,M>::print() says, 'Please implement me.'" << std::endl;
  return;
}

}  // End namespace QUESO

template class QUESO::BetaVectorRV<QUESO::GslVector, QUESO::GslMatrix>;
