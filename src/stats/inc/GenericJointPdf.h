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

#ifndef UQ_GENERIC_JOINT_PROB_DENSITY_H
#define UQ_GENERIC_JOINT_PROB_DENSITY_H

#include <cmath>

#include <queso/JointPdf.h>
#include <queso/Environment.h>
#include <queso/ScalarFunction.h>
#include <queso/BoxSubset.h>

namespace QUESO {

class GslVector;
class GslMatrix;

//*****************************************************
// Generic probability density class [PDF-01]
//*****************************************************
/*!
 * \class GenericJointPdf
 * \brief A class for handling generic joint PDFs.
 *
 * This class allows the mathematical definition of a generic Joint PDF, such as the posterior PDF.*/

template <class V = GslVector, class M = GslMatrix>
class GenericJointPdf : public BaseJointPdf<V,M> {
public:

  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  /*! Instantiates an object of this class given a prefix and a scalar function.
   * The domain of the scalar function is assigned to the protected attribute m_domainSet,
   * and the scalar function is also itself copied to the protected attribute m_scalarFunction.*/
  GenericJointPdf(const char*                           prefix,
                         const BaseScalarFunction<V,M>& scalarFunction);
  //! Destructor
 ~GenericJointPdf();
 //@}

   //! @name Math methods
  //@{
  // See base class (BaseJointPdf) for description.
  double actualValue(const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  double lnValue    (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  double computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const;
  //@}
protected:
  using BaseScalarFunction<V,M>::m_env;
  using BaseScalarFunction<V,M>::m_prefix;
  using BaseScalarFunction<V,M>::m_domainSet;
  using BaseJointPdf<V,M>::m_logOfNormalizationFactor;

  const BaseScalarFunction<V,M>& m_scalarFunction;
};

}  // End namespace QUESO

#endif // UQ_GENERIC_JOINT_PROB_DENSITY_H
