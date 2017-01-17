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

#include <queso/LogNormalVectorRV.h>
#include <queso/LogNormalJointPdf.h>
#include <queso/LogNormalVectorRealizer.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Constructor---------------------------------------
template<class V, class M>
LogNormalVectorRV<V,M>::LogNormalVectorRV(
  const char*                  prefix,
  const VectorSet<V,M>& imageSet,
  const V&                     lawExpVector,
  const V&                     lawVarVector)
  :
  BaseVectorRV<V,M>(((std::string)(prefix)+"gau").c_str(),imageSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering LogNormalVectorRV<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

// begin kemelli 2013-April-23 --------------------------
// LogNormal dist is defined only in (0,inf)

  const BoxSubset<V,M>* imageBox = dynamic_cast<const BoxSubset<V,M>* >(&imageSet);
  double smallerOfMaxValues = imageBox->maxValues().getMinValue();
  double smallerOfMinValues = imageBox->minValues().getMinValue();

 if( smallerOfMinValues < 0 )
 {
   std::cerr << "In LogNormalVectorRV<V,M>::constructor()\n"
       << "LogNormal distribution is only defined in (0, infinity).\n"
       << "The data provided is: \n"
       << *imageBox
         << "Sampling will not cover all interval.\n"
         << std::endl;


    queso_require_greater_equal_msg(smallerOfMaxValues, 0, "invalid input: LogNormal distribution is only defined in (0, infinity), and min(m_maxValues)<0");

 }
// end kemelli 2013-April-23 --------------------------

  m_pdf = new LogNormalJointPdf<V,M>(m_prefix.c_str(),
                                            m_imageSet,
                                            lawExpVector,
                                            lawVarVector);

  M lowerCholLawCovMatrix(lawVarVector);
  int iRC = lowerCholLawCovMatrix.chol();
  lowerCholLawCovMatrix.zeroUpper(false);
  if (iRC) {
    std::cerr << "In LogNormalVectorRV<V,M>::constructor() [1]: chol failed, will use svd\n";
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In LogNormalVectorRV<V,M>::constructor() [1]: chol failed; will use svd; lawVarVector contents are\n";
      *m_env.subDisplayFile() << lawVarVector; // FIX ME: might demand parallelism
      *m_env.subDisplayFile() << std::endl;
    }
    M matLaw(lawVarVector);
    M matU  (lawVarVector);
    M matVt (m_imageSet.vectorSpace().zeroVector());
    V vecS  (m_imageSet.vectorSpace().zeroVector());
    iRC = matLaw.svd(matU,vecS,matVt);
    queso_require_msg(!(iRC), "Cholesky decomposition of covariance matrix failed.");

    vecS.cwSqrt();
    m_realizer = new LogNormalVectorRealizer<V,M>(m_prefix.c_str(),
                                                         m_imageSet,
                                                         lawExpVector,
                                                         matU,
                                                         vecS, // already square rooted
                                                         matVt);
  }
  else {
    m_realizer = new LogNormalVectorRealizer<V,M>(m_prefix.c_str(),
                                                        m_imageSet,
                                                        lawExpVector,
                                                        lowerCholLawCovMatrix);
  }

  m_subCdf     = NULL; // FIX ME: complete code
  m_unifiedCdf = NULL; // FIX ME: complete code
  m_mdf        = NULL; // FIX ME: complete code

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving LogNormalVectorRV<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V, class M>
LogNormalVectorRV<V,M>::~LogNormalVectorRV()
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
LogNormalVectorRV<V,M>::print(std::ostream& os) const
{
  os << "LogNormalVectorRV<V,M>::print() says, 'Please implement me.'" << std::endl;
  return;
}


}  // End namespace QUESO

template class QUESO::LogNormalVectorRV<QUESO::GslVector,QUESO::GslMatrix>;
