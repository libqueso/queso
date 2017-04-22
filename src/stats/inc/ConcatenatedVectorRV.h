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

#ifndef UQ_CONCATENATED_VECTOR_RV_H
#define UQ_CONCATENATED_VECTOR_RV_H

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
// Concatenated class [RV-11]
//*****************************************************
/*!
 * \class ConcatenatedVectorRV
 * \brief A class representing concatenated vector RVs.
 *
 * This class allows the user to concatenate two vector RV of different types and to generate realizations
 * (samples) from this concatenated vector RV. It is used, for instance, to concatenate priors from two or
 * more RVs, where one of them has a uniform distribution whereas the other one(s) has a Gaussian distribution. */

template <class V = GslVector, class M = GslMatrix>
class ConcatenatedVectorRV : public BaseVectorRV<V,M> {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*! Concatenates two RVs: \c rv1 and \c rv2 into one vector RV, given a prefix and the image set of the vector RV.*/
  ConcatenatedVectorRV(const char*                     prefix,
                              const BaseVectorRV<V,M>& rv1,
                              const BaseVectorRV<V,M>& rv2,
                              const VectorSet<V,M>&    imageSet);

  //! Constructor
  /*! Concatenates a sequence of RVs, given by: <c> std::vector<const BaseVectorRV<V,M>* >& rvs </c>
   * into one single vector RV, given a prefix and the image set of the resulting vector RV.*/
  ConcatenatedVectorRV(const char*                                          prefix,
                              const std::vector<const BaseVectorRV<V,M>* >& rvs,
                              const VectorSet<V,M>&                         imageSet);

  //! Virtual destructor
  virtual ~ConcatenatedVectorRV();
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

  std::vector<const BaseVectorRV      <V,M>* > m_rvs;
  std::vector<const BaseJointPdf      <V,M>* > m_pdfs;
  std::vector<const BaseVectorRealizer<V,M>* > m_realizers;
};

}  // End namespace QUESO

#endif // UQ_CONCATENATED_VECTOR_RV_H
