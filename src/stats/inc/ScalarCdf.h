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

#ifndef UQ_SCALAR_CUMULATIVE_DISTRIBUTION_FUNCTION_H
#define UQ_SCALAR_CUMULATIVE_DISTRIBUTION_FUNCTION_H

#include <queso/StdOneDGrid.h>
#include <queso/Environment.h>
#include <math.h>

namespace QUESO {

/*! \file ScalarCdf.h
 * \brief Classes to accommodate a cumulative distribution function.
 *
 * \class BaseScalarCdf
 * \brief A templated (base) class for handling CDFs.
 *
 * This class allows the mathematical definition of a cumulative distribution function (CDF),
 * which is a scalar function such as * \f$ f: B \subset R \rightarrow C \subset R \f$; ie a
 * function of one or more variables that has always one-dimensional range, which for the
 * specific CDF case, the image set \f$ C = [0,1]\f$. The CDF describes the probability that
 * a real-valued random variable X with a given probability distribution will be found at a
 * value less than or equal to x. In the case of a continuous distribution, it gives the area
 * under the probability density function (PDF) from minus infinity to x.*/

template <class T>
class BaseScalarCdf {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  /*! Instantiates an object of the class  given a prefix and the environment.*/
  BaseScalarCdf(const BaseEnvironment& env, const char* prefix);

  //! Virtual destructor.
  virtual ~BaseScalarCdf();
  //@}

  //! @name Environment methods
  //@{
  //! Environment.  Access to private attribute m_env.
  const BaseEnvironment&  env             () const;

  //! Access to private attribute m_prefix.
  const std::string&             prefix          () const;
  //@}

  //! @name Mathematical methods
  //@{
  //! Returns the value of the CDF at \c paramValue. See template specialization.
  virtual double                  value           (T             paramValue          ) const = 0;

  //! Returns the position of a given value of CDF. See template specialization.
  virtual T                       inverse         (double        cdfValue            ) const = 0;

  //! Returns the support (image) of the CDF between two horizontal values (domain). See template specialization.
  virtual void                    getSupport      (T& minHorizontal, T& maxHorizontal) const = 0;
  //@}
   //! @name I/O methods
  //@{
  //! Prints the CDF. See template specialization.
  virtual void                    print           (std::ostream& os                  ) const = 0;
  friend std::ostream& operator<< (std::ostream& os,
      const BaseScalarCdf<T>& obj) {
    obj.print(os);
    return os;
  }

  //! Writes the CDF of an allowed sub-environment to a file.
  /*! This function does nothing and should \n not be called by the user.*/
  virtual void                    subWriteContents(const std::string&            varNamePrefix,
                                                   const std::string&            fileName,
                                                   const std::string&            fileType,
                                                   const std::set<unsigned int>& allowedSubEnvIds) const;
  //@}
protected:
  const BaseEnvironment& m_env;
        std::string             m_prefix;
  mutable T                     m_minHorizontal;
  mutable T                     m_maxHorizontal;
};

//---------------------------------------------------
// Method outside either class definition------------
//---------------------------------------------------

//*****************************************************
// Horizontal distance
//*****************************************************

//! It calculated the maximum horizontal distance between two CDFs.
template <class T>
double
horizontalDistance(const BaseScalarCdf<T>& cdf1,
                   const BaseScalarCdf<T>& cdf2,
                   double epsilon);

}  // End namespace QUESO

#endif // UQ_SCALAR_CUMULATIVE_DISTRIBUTION_FUNCTION_H
