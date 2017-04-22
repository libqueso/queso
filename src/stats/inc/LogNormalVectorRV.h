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

#ifndef UQ_LOGNORM_VECTOR_RV_H
#define UQ_LOGNORM_VECTOR_RV_H

#include <queso/VectorRV.h>
#include <queso/VectorSpace.h>
#include <queso/JointPdf.h>
#include <queso/VectorRealizer.h>
#include <queso/VectorCdf.h>
#include <queso/VectorMdf.h>
#include <queso/SequenceOfVectors.h>

namespace QUESO {

class GslVector;
class GslMatrix;

//*****************************************************
// LogNormal class [RV-10]
//*****************************************************
/*!
 * \class LogNormalVectorRV
 * \brief A class representing a LogNormal vector RV.
 *
 * This class allows the user to compute the value of a LogNormal PDF and to generate realizations
 * (samples) from it.\n
 *
 * The probability density function of a log-normal distribution is:
 *   \f[ f(x|\mu,\sigma) = \frac{1}{x \sigma \sqrt{2 \pi}}\, e^{-\frac{(\ln x - \mu)^2}{2\sigma^2}}, x>0 \f]
 * where  the parameters denoted \f$ \mu \f$ and \f$ \sigma \f$ are, respectively, the mean and standard
 * deviation of the  variable's natural logarithm; and \c x>0.  */

template <class V = GslVector, class M = GslMatrix>
class LogNormalVectorRV : public BaseVectorRV<V,M> {
public:
    //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*! Construct a LogNormal vector RV with mean \c lawExpVector and diagonal covariance matrix
   * \c lawVarVector whose variates live in \c imageSet.*/
  LogNormalVectorRV(const char*                  prefix,
                           const VectorSet<V,M>& imageSet,
                           const V&                     lawExpVector,
                           const V&                     lawVarVector);

  //! Virtual destructor
  virtual ~LogNormalVectorRV();
  //@}

  //! @name I/O methods
  //@{
  //! TODO: Prints the vector RV.
  /*! \todo: implement me!*/
  void print(std::ostream& os) const;
  //@}

private:
  using BaseVectorRV<V,M>::m_env;
  using BaseVectorRV<V,M>::m_prefix;
  using BaseVectorRV<V,M>::m_imageSet;
  using BaseVectorRV<V,M>::m_pdf;
  using BaseVectorRV<V,M>::m_realizer;
  using BaseVectorRV<V,M>::m_subCdf;
  using BaseVectorRV<V,M>::m_unifiedCdf;
  using BaseVectorRV<V,M>::m_mdf;
};

}  // End namespace QUESO

#endif // UQ_LOGNORM_VECTOR_RV_H
