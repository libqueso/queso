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

#ifndef UQ_VECTOR_RV_H
#define UQ_VECTOR_RV_H

#include <queso/Environment.h>
#include <queso/VectorSet.h>
#include <queso/VectorSpace.h>
#include <queso/JointPdf.h>
#include <queso/VectorRealizer.h>
#include <queso/VectorCdf.h>
#include <queso/VectorMdf.h>
#include <queso/SequenceOfVectors.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*! \file VectorRV.h
 * \brief A templated class for handling vector random variables (RV).
 *
 * \class BaseVectorRV
 * \brief A templated base class for handling vector RV.
 *
 * This class allows two basic but quite crucial functionalities: to compute the value of the
 * PDF of a random variable (RV) at a point and to generate realizations (samples) from such PDF.
 */

template <class V = GslVector, class M = GslMatrix>
class BaseVectorRV {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*! Constructs a new instance of BaseVectorRV, given a prefix and the image set of the
   * vector RV.
   */
  BaseVectorRV(const char*                  prefix,
                      const VectorSet<V,M>& imageSet);

  //! Virtual destructor.
  virtual ~BaseVectorRV();
  //@}

  //! @name Random variable-handling methods
  //@{
  //! QUESO environment; access to private attribute m_env.
  const   BaseEnvironment&         env       () const;

  //! Image set of the vector RV; access to private attribute m_imageSet.
  const   VectorSet         <V,M>& imageSet  () const;

  //! Posterior Density Function of the vector RV; access to private attribute m_pdf.
  const   BaseJointPdf      <V,M>& pdf       () const;

  //! Returns true iff this RV has the ability to produce realizations (samples)
          bool                  has_realizer () const;

  //! Finds a realization (sample) of the PDF of this vector RV; access to private attribute m_realizer.
  const   BaseVectorRealizer<V,M>& realizer  () const;

  //! Finds the Cumulative Distribution Function of this vector RV, considering only the sub-sequence of data; access to private attribute m_subCdf.
  const   BaseVectorCdf     <V,M>& subCdf    () const;

  //! Finds the Cumulative Distribution Function of this vector RV, considering the unified sequence of data; access to private attribute m_unifiedCdf.
  const   BaseVectorCdf     <V,M>& unifiedCdf() const;

  //! Finds the Marginal Density Function of this vector RV; access to private attribute m_mdf.
  const   BaseVectorMdf     <V,M>& mdf       () const;
  //@}

  //! @name I/O methods
  //@{
  //! TODO: Prints the vector RV.
  /*! \todo: implement me!*/
  virtual void                            print     (std::ostream& os) const = 0;
  friend std::ostream& operator<<(std::ostream& os,
      const BaseVectorRV<V,M>& obj) {
    obj.print(os);
    return os;
  }
  //@}

#ifdef QUESO_HAS_ANN
  virtual double                          estimateENT_ANN() const;
  /*
  virtual double                          estimateENT_ANN( unsigned int k, double eps ) const;
  virtual double                          estimateENTSubset_ANN( const unsigned int dimSel[] ) const;
  virtual double                          estimateENTSubset_ANN( const unsigned int dimSel[], unsigned int k, double eps ) const;
  */
#endif // QUESO_HAS_ANN
  //@}
protected:
  const   BaseEnvironment&         m_env;
          std::string                     m_prefix;
  const   VectorSet         <V,M>& m_imageSet;
          BaseJointPdf      <V,M>* m_pdf;
	  	  BaseVectorRealizer<V,M>* m_realizer;
  const   BaseVectorCdf     <V,M>* m_subCdf;
  const   BaseVectorCdf     <V,M>* m_unifiedCdf;
  const   BaseVectorMdf     <V,M>* m_mdf;
};

//---------------------------------------------------
// Method declared outside class definition ---------
//---------------------------------------------------
template <class P_V, class P_M, class Q_V, class Q_M>
void
ComputeCovCorrMatricesBetweenVectorRvs(
  const BaseVectorRV<P_V,P_M>& paramRv,
  const BaseVectorRV<Q_V,Q_M>& qoiRv,
        unsigned int                  localNumSamples,
        P_M&                          pqCovMatrix,
        P_M&                          pqCorrMatrix);

}  // End namespace QUESO

#endif // UQ_VECTOR_RV_H
