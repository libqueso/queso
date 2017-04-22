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

#ifndef UQ_GENERIC_VECTOR_RV_H
#define UQ_GENERIC_VECTOR_RV_H

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
// Generic class [RV-01]
//*****************************************************
 /*! \class GenericVectorRV
 * \brief A templated class for handling generic vector RVs.
 *
 * This class allows the user to compute the value of the PDF of a generic random variable (RV)
 * and to generate realizations (samples) from such PDF.  This is the class used by QUESO to
 * store the solution of an statistical inverse problem. */

template <class V = GslVector, class M = GslMatrix>
class GenericVectorRV : public BaseVectorRV<V,M> {
public:
    //! @name Constructor/Destructor methods
  //@{
  //! Default constructor
  /*! Constructs a new instance, given a prefix and the image set of the vector RV.*/
  GenericVectorRV(const char*                           prefix,
                         const VectorSet         <V,M>& imageSet);

  //! Constructor
  /*! Constructs a new instance, given all the attributes that characterize the vector RV: prefix, image set, pdf, etc.*/
  GenericVectorRV(const char*                           prefix,
                         const VectorSet         <V,M>& imageSet,
                         BaseJointPdf      <V,M>& pdf,
                         BaseVectorRealizer<V,M>& realizer,
                         const BaseVectorCdf     <V,M>& subCdf,
                         const BaseVectorCdf     <V,M>& unifiedCdf,
                         const BaseVectorMdf     <V,M>& mdf);
  //! Virtual destructor
  virtual ~GenericVectorRV();
  //@}

    //! @name Random variable-handling methods
  //@{
  //! Sets the PDF of \c this vector RV  to \c pdf.
  void setPdf       (BaseJointPdf      <V,M>& pdf       );

  //! Sets the realizer of \c this vector RV  to \c realizer.
  void setRealizer  (BaseVectorRealizer<V,M>& realizer  );

  //! Sets the CDF of the sub-sequence of \c this vector RV  to \c subCdf.
  void setSubCdf    (BaseVectorCdf     <V,M>& subCdf    );

  //! Sets the CDF of the unified sequence of \c this vector RV  to \c unifiedCdf.
  void setUnifiedCdf(BaseVectorCdf     <V,M>& unifiedCdf);

  //! Sets the MDF of  \c this vector RV  to \c Mdf.
  void setMdf       (BaseVectorMdf     <V,M>& mdf       );
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

#endif // UQ_GENERIC_VECTOR_RV_H
