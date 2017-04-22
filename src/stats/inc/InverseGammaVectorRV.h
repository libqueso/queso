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

#ifndef UQ_INVGAMMA_VECTOR_RV_H
#define UQ_INVGAMMA_VECTOR_RV_H

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
// InverseGamma class [RV-07]
//*****************************************************
/*!
 * \class InverseGammaVectorRV
 * \brief A class representing a vector RV constructed via Inverse Gamma distribution.
 *
 * This class allows the user to compute the value of a Inverse Gamma PDF and to generate realizations
 * (samples) from it.\n
 *
 * The Inverse Gamma probability density function  is defined over the support x > 0 given a shape parameters
 * \b a and a scale parameter \b b is:
 *  \f[ y=f(x| a,b)= \frac{b^a}{\Gamma(a)} x^{-a - 1}\exp\left(-\frac{b}{x}\right)\f]
 * where \f$ \Gamma(.) \f$ is the Gamma function:
 * \f[  B(a,b)=\frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}=\frac{(a-1)!(b-1)!}{(a+b-1)!}.\f]
 * The parameters \b a and \b b must all be positive, and the values \c x  must lie on the
 * interval \f$ (0, \infty)\f$. */
template <class V = GslVector, class M = GslMatrix>
class InverseGammaVectorRV : public BaseVectorRV<V,M> {
public:

  //! @name Constructor/Destructor methods
  //@{
  //! Default Constructor
  /*! Construct an Inverse Gamma vector RV with parameters \c alpha>0  and \c beta>0, whose variates live in \c imageSet.
   * The constructor will check whether or not the data provided via \c imageSet belongs to
   * \f$ (0, \infty)\f$, which is a requirement imposed by the Inverse Gamma distribution. If this
   * condition is not satisfied, an error  message will be displayed and the program will exit. */
  InverseGammaVectorRV(const char*                  prefix,
                              const VectorSet<V,M>& imageSet,
                              const V&                     alpha,
                              const V&                     beta);
  //! Virtual destructor
  virtual ~InverseGammaVectorRV();
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

#endif // UQ_INVGAMMA_VECTOR_RV_H
