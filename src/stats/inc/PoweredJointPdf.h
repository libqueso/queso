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

#ifndef UQ_POWERED_JOINT_PROB_DENSITY_H
#define UQ_POWERED_JOINT_PROB_DENSITY_H

#include <cmath>

#include <queso/JointPdf.h>
#include <queso/Environment.h>
#include <queso/ScalarFunction.h>
#include <queso/BoxSubset.h>

namespace QUESO {

class GslVector;
class GslMatrix;

//*****************************************************
// Powered probability density class [PDF-08]
//*****************************************************
/*!
 * \class PoweredJointPdf
 * \brief A class for handling a powered joint PDFs.
 *
 * This class allows the mathematical definition of a Powered Joint PDF.*/

template <class V = GslVector, class M = GslMatrix>
class PoweredJointPdf : public BaseJointPdf<V,M> {
public:
    //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*! Constructs a new object of the class, given a prefix, the domain set and the exponent of
   * the powered PDF.  */
  PoweredJointPdf(const char*                     prefix,
                         const BaseJointPdf<V,M>& srcDensity,
                               double                    exponent);
  //! Destructor
 ~PoweredJointPdf();
 //@}

  //! @name Math methods
  //@{
  //! Actual value of the powered PDF.
  /*! Finds the actual value using BaseJointPdf::actualValue() and apply it to the power of
   * \c this PDF, which given by \c exponent, and multiplies it by the normalization factor, which is
   * given by exp(m_logOfNormalizationFactor).*/
  double actualValue          (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;

  //! Logarithm of the value of the powered PDF.
    /*! Finds the logarithm of actual value using BaseJointPdf::lnValue() and multiplies by the power of
   * \c this PDF, which given by \c exponent, and then adds the normalization factor, which is
   * given by m_logOfNormalizationFactor.*/
  double lnValue              (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;

  //! Sets a value to be used in the normalization style of the powered PDF (ie, protected attribute m_srcDensity).
  void   setNormalizationStyle(unsigned int value) const;

  //! TODO: Computes the logarithm of the normalization factor.
  /*! \todo: implement me!*/
  double computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const;
  //@}

protected:
  using BaseScalarFunction<V,M>::m_env;
  using BaseScalarFunction<V,M>::m_prefix;
  using BaseScalarFunction<V,M>::m_domainSet;
  using BaseJointPdf<V,M>::m_logOfNormalizationFactor;

  const BaseJointPdf<V,M>& m_srcDensity;
  double                          m_exponent;
};

}  // End namespace QUESO

#endif // UQ_POWERED_JOINT_PROB_DENSITY_H
