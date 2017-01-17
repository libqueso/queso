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

#ifndef UQ_GENERIC_VECTOR_MARGINAL_DENSITY_FUNCTION_H
#define UQ_GENERIC_VECTOR_MARGINAL_DENSITY_FUNCTION_H

#include <queso/VectorMdf.h>
#include <queso/ArrayOfOneDGrids.h>
#include <queso/ArrayOfOneDTables.h>
#include <queso/Environment.h>
#include <math.h>

namespace QUESO {

class GslVector;
class GslMatrix;

//*****************************************************
// Generic marginal density function
//*****************************************************
/*!\class GenericVectorMdf
 * \brief A class for handling generic MDFs of vector functions. */

template <class V = GslVector, class M = GslMatrix>
class GenericVectorMdf : public BaseVectorMdf<V,M> {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor.
  /*! Instantiates an object of the class given a prefix, the domain set, and a routine
   * (acting as a math function). */
  GenericVectorMdf(const char*                    prefix,
                          const VectorSet<V,M>& domainSet,
                          double (*routinePtr)(const V& paramValues, const void* routineDataPtr, V& mdfVec),
                          const void* routineDataPtr);
  //! Destructor
 ~GenericVectorMdf();
  //@}

    //! @name Mathematical method
  //@{
  //! Finds the values of the vector MDF at each element of \c paramValues, by calling \c m_routinePtr, and saves it at \c mdfValues.
  void values(const V& paramValues, V& mdfVec) const;
  //@}

  //! @name I/O method
  //@{
  //! TODO: Prints the vector MDF.
  /*! \todo: implement me!*/
  void print (std::ostream& os)                const;
  //@}

protected:
  double (*m_routinePtr)(const V& paramValues, const void* routineDataPtr, V& mdfVec);
  const void* m_routineDataPtr;

  using BaseVectorMdf<V,M>::m_env;
  using BaseVectorMdf<V,M>::m_prefix;
  using BaseVectorMdf<V,M>::m_domainSet;
};

}  // End namespace QUESO

#endif // UQ_GENERIC_VECTOR_MARGINAL_DENSITY_FUNCTION_H
