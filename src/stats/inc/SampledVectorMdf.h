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

#ifndef UQ_SAMPLED_VECTOR_MARGINAL_DENSITY_FUNCTION_H
#define UQ_SAMPLED_VECTOR_MARGINAL_DENSITY_FUNCTION_H

#include <queso/VectorMdf.h>
#include <queso/ArrayOfOneDGrids.h>
#include <queso/ArrayOfOneDTables.h>
#include <queso/Environment.h>
#include <math.h>

namespace QUESO {

class GslVector;
class GslMatrix;

//*****************************************************
// Sampled marginal density function class
//*****************************************************
/*!
 * \class SampledVectorMdf
 * \brief A class for handling sampled vector MDFs.
 *
 * This class implements a sampled vector marginal density function (MDF), given
 * the grid points where it will be sampled and it returns its values.*/

template <class V = GslVector, class M = GslMatrix>
class SampledVectorMdf : public BaseVectorMdf<V,M> {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  /*! Instantiates an object of the class given a prefix and the grid points
   * where it will be sampled/evaluated.*/
  SampledVectorMdf(const char*                          prefix,
                          const ArrayOfOneDGrids <V,M>& oneDGrids,
                          const ArrayOfOneDTables<V,M>& mdfValues);
  //! Destructor
  ~SampledVectorMdf();
  //@}

  //! @name Mathematical methods
  //@{
  //! TODO: Returns the values of the vector MDF at each element of \c paramValues.
  /*! \todo: implement me!*/
  void values(const V& paramValues, V& mdfVec) const;
  //@}

  //! @name I/O methods
  //@{
  //! Prints the vector MDF (values of the grid points and of the MDF at such grid points).
  void print (std::ostream& os)                const;
  //@}

protected:
  using BaseVectorMdf<V,M>::m_env;
  using BaseVectorMdf<V,M>::m_prefix;
  using BaseVectorMdf<V,M>::m_domainSet;

  const ArrayOfOneDGrids <V,M>& m_oneDGrids;
  const ArrayOfOneDTables<V,M>& m_mdfValues;
};

}  // End namespace QUESO

#endif // UQ_SAMPLED_VECTOR_MARGINAL_DENSITY_FUNCTION_H
