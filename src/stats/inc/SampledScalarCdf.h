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

#ifndef UQ_SAMPLED_SCALAR_CUMULATIVE_DISTRIBUTION_FUNCTION_H
#define UQ_SAMPLED_SCALAR_CUMULATIVE_DISTRIBUTION_FUNCTION_H

#include <queso/StdOneDGrid.h>
#include <queso/Environment.h>
#include <queso/ScalarCdf.h>
#include <math.h>

namespace QUESO {

//*****************************************************
// Sampled cumulative distribution function class
//*****************************************************
/*!
 * \class SampledScalarCdf
 * \brief A class for handling sampled CDFs.
 *
 * This class implements a sampled cumulative distribution function (CDF), given
 * the grid points where it will be sampled and its resulting values.*/

template <class T>
class SampledScalarCdf : public BaseScalarCdf<T> {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  /*! Instantiates an object of the class given a prefix, the environment, the grid points
   * where it will be sampled/evaluated and its resulting values.*/
  SampledScalarCdf(const BaseEnvironment& env,
                          const char*                   prefix,
                          const BaseOneDGrid<T>& cdfGrid,
                          const std::vector<double>&    cdfValues);
  //! Destructor.
 ~SampledScalarCdf();
  //@}

  //! @name Mathematical methods
  //@{
  //! Returns the value of the CDF at \c paramValue.
  double value           (T                       paramValue) const;

  //! Returns the position of a given value of CDF.
  T      inverse         (double                  cdfValue  ) const;

  //! Returns the support (image) of the CDF between two horizontal values (domain).
  void   getSupport      (T& minHorizontal, T& maxHorizontal) const;
  //@}

  //! @name I/O methods
  //@{
  //! Prints the CDF (values of the grid points and of the CDF at such grid points).
  void   print           (std::ostream&           os        ) const;


  //!Writes the CDF of an allowed sub-environment to a file.
  /*! It will write the data in  Octave/Matlab compatible format.*/
  void   subWriteContents(const std::string&            varNamePrefix,
                          const std::string&            fileName,
                          const std::string&            fileType,
                          const std::set<unsigned int>& allowedSubEnvIds) const;
 //@}
protected:
  using BaseScalarCdf<T>::m_env;
  using BaseScalarCdf<T>::m_prefix;
  using BaseScalarCdf<T>::m_minHorizontal;
  using BaseScalarCdf<T>::m_maxHorizontal;

  const BaseOneDGrid<T>& m_cdfGrid;
  const std::vector<double>&    m_cdfValues;
//std::vector<double>&          m_sortedCdfValues;
};

}  // End namespace QUESO

#endif // UQ_SAMPLED_SCALAR_CUMULATIVE_DISTRIBUTION_FUNCTION_H
