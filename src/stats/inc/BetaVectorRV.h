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

#ifndef UQ_BETA_VECTOR_RV_H
#define UQ_BETA_VECTOR_RV_H

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
// Beta class [RV-05]
//*****************************************************
/*!
 * \class BetaVectorRV
 * \brief A class representing a vector RV constructed via Beta distribution.
 *
 * This class allows the user to compute the value of a Beta PDF and to generate realizations
 * (samples) from it.\n
 *
 * The beta probability density function for a given value x and given pair of parameters
 * \b a and \b b is:
 *  \f[ y=f(x|a,b)= \frac{1}{B(a,b)} x^{a-1}(1-x)^{b-1}, \f]
 * where <b>B(·)</b> is the Beta function:
 * \f[  B(a,b)=\frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}=\frac{(a-1)!(b-1)!}{(a+b-1)!}.\f]
 * The parameters \b a and \b b must all be positive, and the values \c x  must lie on the
 * interval [0, 1].
 */


template <class V = GslVector, class M = GslMatrix>
class BetaVectorRV : public BaseVectorRV<V,M> {
public:
    //! @name Constructor/Destructor methods
  //@{
  //! Default Constructor
  /*! Construct a Beta vector RV with parameters \c a>0  and \c b>0, whose variates live in \c imageSet.
   * The constructor will check whether or not the data provided via \c imageSet belongs to [0,1], which
   * is a requirement imposed by the Beta distribution. If this condition is not satisfied, an error
   * message will be displayed and the program will exit. */
  BetaVectorRV(const char*                  prefix,
                      const VectorSet<V,M>& imageSet,
                      const V&                     alpha,
                      const V&                     beta);

  //! Virtual destructor
  virtual ~BetaVectorRV();
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

#endif // UQ_BETA_VECTOR_RV_H
