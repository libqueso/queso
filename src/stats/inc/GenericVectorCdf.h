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

#ifndef UQ_GENERIC_VECTOR_CUMULATIVE_DISTRIBUTION_FUNCTION_H
#define UQ_GENERIC_VECTOR_CUMULATIVE_DISTRIBUTION_FUNCTION_H

#include <queso/VectorCdf.h>
#include <queso/ArrayOfOneDGrids.h>
#include <queso/ArrayOfOneDTables.h>
#include <queso/ScalarCdf.h>
#include <queso/Environment.h>
#include <math.h>

namespace QUESO {

class GslVector;
class GslMatrix;

//*****************************************************
// Generic cumulative distribution function class
//*****************************************************
/*!
 * \class GenericVectorCdf
 * \brief A class for handling generic vector CDFs.
 *
 * This class \b will implement a generic vector cumulative distribution function (CDF).*/

template <class V = GslVector, class M = GslMatrix>
class GenericVectorCdf : public BaseVectorCdf<V,M> {
public:
    //! @name Constructor/Destructor methods
  //@{
  //! Constructor.
  /*! Instantiates an object of the class given a prefix, the support of the related-PDF, and
   * a routine that calculates data (like a math function). */
  GenericVectorCdf(const char*                  prefix,
                          const VectorSet<V,M>& pdfSupport,
                          double (*routinePtr)(const V& paramValues, const void* routineDataPtr, V& cdfVec),
                          const void*                  routineDataPtr);
  //! Destructor
  ~GenericVectorCdf();
  //@}

    //! @name Mathematical method
  //@{
  //! TODO: Returns the values of the vector CDF at each element of \c paramValues, by calling \c m_routinePtr.
  void values(const V& paramValues, V& cdfVec) const;
  //@}

  //! @name I/O method
  //@{
  //! TODO: Prints the vector CDF.
  /*! \todo: implement me!*/
  void print (std::ostream& os) const;
  //@}

protected:
  double (*m_routinePtr)(const V& paramValues, const void* routineDataPtr, V& cdfVec);
  const void* m_routineDataPtr;

  using BaseVectorCdf<V,M>::m_env;
  using BaseVectorCdf<V,M>::m_prefix;
  using BaseVectorCdf<V,M>::m_pdfSupport;
};

}  // End namespace QUESO

#endif // UQ_GENERIC_VECTOR_CUMULATIVE_DISTRIBUTION_FUNCTION_H
