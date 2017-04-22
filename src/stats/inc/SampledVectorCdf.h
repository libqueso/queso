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

#ifndef UQ_SAMPLED_VECTOR_CUMULATIVE_DISTRIBUTION_FUNCTION_H
#define UQ_SAMPLED_VECTOR_CUMULATIVE_DISTRIBUTION_FUNCTION_H

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
// Sampled cumulative distribution function class
//*****************************************************
/*!
 * \class SampledVectorCdf
 * \brief A class for handling sampled vector CDFs.
 *
 * This class implements a sampled vector cumulative distribution function (CDF), given
 * the grid points where it will be sampled and it returns its values.*/

template <class V = GslVector, class M = GslMatrix>
class SampledVectorCdf : public BaseVectorCdf<V,M> {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  /*! Instantiates an object of the class given a prefix and the grid points
   * where it will be sampled/evaluated.*/
  SampledVectorCdf(const char*                          prefix,
                          const ArrayOfOneDGrids <V,M>& oneDGrids,
                          const ArrayOfOneDTables<V,M>& cdfValues);
  //! Destructor
  ~SampledVectorCdf();
  //@}

  //! @name Mathematical methods
  //@{
  //! TODO: Returns the values of the vector CDF at each element of \c paramValues.
  /*! \todo: implement me!*/
  void                          values(const V& paramValues, V& cdfVec) const;

  //! Returns a scalar CDF stored at row \c rowId of the vector CDF.
  const BaseScalarCdf<double>& cdf   (unsigned int rowId)              const;
  //@}

    //! @name I/O methods
  //@{
  //! Prints the vector CDF (values of the grid points and of the CDF at such grid points).
  void                          print (std::ostream& os)                const;

  //!Writes the CDF of an allowed sub-environment to a file.
  void                          subWriteContents(const std::string&            varNamePrefix,
             const std::string&            fileName,
             const std::string&            fileType,
             const std::set<unsigned int>& allowedSubEnvIds) const;
  //@}
protected:
  using BaseVectorCdf<V,M>::m_env;
  using BaseVectorCdf<V,M>::m_prefix;
  using BaseVectorCdf<V,M>::m_pdfSupport;

  DistArray<SampledScalarCdf<double>*> m_cdfs;
};

}  // End namespace QUESO

#endif // UQ_SAMPLED_VECTOR_CUMULATIVE_DISTRIBUTION_FUNCTION_H
