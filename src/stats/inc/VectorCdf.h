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

#ifndef UQ_VECTOR_CUMULATIVE_DISTRIBUTION_FUNCTION_H
#define UQ_VECTOR_CUMULATIVE_DISTRIBUTION_FUNCTION_H

#include <queso/ArrayOfOneDGrids.h>
#include <queso/ArrayOfOneDTables.h>
#include <queso/ScalarCdf.h>
#include <queso/SampledScalarCdf.h>
#include <queso/Environment.h>
#include <math.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*! \file VectorCdf.h
 * \brief Classes to accommodate a cumulative distribution function of a vector RV.
 *
 * \class BaseVectorCdf
 * \brief A templated (base) class for handling CDFs of vector functions.
 *
 * In many applications is necessary to consider the properties of two or more RVs
 * simultaneously (represented within QUESO as Random Vectors, via BaseVectorRV
 * and derived classes). When dealing simultaneously with more than one RV, ie, a vector
 * RV, the joint cumulative distribution function must also be defined. This class handles
 * the CDFs of vector RV, which are referred to as  multivariate/vector/joint CDFs.*/

template <class V = GslVector, class M = GslMatrix>
class BaseVectorCdf {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  /*! Instantiates an object of the class  given a prefix and the support (image set) of the PDF that is
   * related to this CDF (recall that the CDF of a continuous RV is the integral of the PDF of that RV).*/
  BaseVectorCdf(const char*                  prefix,
		       const VectorSet<V,M>& pdfSupport);

  //! Virtual destructor.
  virtual ~BaseVectorCdf();
  //@}

  //! @name Mathematical methods
  //@{
  //! Returns the image set (support) of the PDF; access to protected attribute \c m_pdfSupport.
  const VectorSet<V,M>&        pdfSupport      ()                                const;

  //! Finds the value of the vector CDF at each element of \c paramValue, and saves it in \c cdfVec. See template specialization.
  virtual void                                values          (const V& paramValues, V& cdfVec) const = 0;

  //!
  virtual const BaseScalarCdf<double>& cdf             (unsigned int rowId)              const = 0;
  //@}

  //! @name I/O methods
  //@{
  //! Prints the vector CDF. See template specialization.
  virtual void                                print           (std::ostream& os)                const = 0;
  friend std::ostream& operator<< (std::ostream& os,
      const BaseVectorCdf<V,M>& obj) {
    obj.print(os);
    return os;
  }

  //! Writes the CDF of an allowed sub-environment to a file.
  /*! This function does nothing and should \n not be called by the user.*/
  virtual void                                subWriteContents(const std::string&            varNamePrefix,
                                                               const std::string&            fileName,
                                                               const std::string&            fileType,
                                                               const std::set<unsigned int>& allowedSubEnvIds) const;
  //@}
protected:

  const   BaseEnvironment& m_env;
          std::string             m_prefix;
  const   VectorSet<V,M>&  m_pdfSupport;
};

//---------------------------------------------------
// Method outside either class definition------------
//---------------------------------------------------
//! It calculated the maximum horizontal distances between two vector CDFs.
//
template <class V, class M>
void
horizontalDistances(const BaseVectorCdf<V,M>& cdf1,
                    const BaseVectorCdf<V,M>& cdf2,
                    const V& epsilonVec,
                    V&       distances);

}  // End namespace QUESO

#endif // UQ_VECTOR_CUMULATIVE_DISTRIBUTION_FUNCTION_H
