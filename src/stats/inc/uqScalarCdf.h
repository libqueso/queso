//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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
// 
// $Id$
//
//--------------------------------------------------------------------------

#ifndef __UQ_SCALAR_CUMULATIVE_DISTRIBUTION_FUNCTION_H__
#define __UQ_SCALAR_CUMULATIVE_DISTRIBUTION_FUNCTION_H__

#include <uqOneDGrid.h>
#include <uqEnvironment.h>
#include <math.h>

namespace QUESO {

//*****************************************************
// Classes to accommodate a cumulative distribution function
//*****************************************************

//*****************************************************
// Base class
//*****************************************************
/*! \file uqScalarCdf.h
 * \brief Classes to accommodate a cumulative distribution function.
 * 
 * \class uqBaseScalarCdfClass
 * \brief A templated (base) class for handling CDFs.
 *
 * This class allows the mathematical definition of a cumulative distribution function (CDF), 
 * which is a scalar function such as * \f$ f: B \subset R \rightarrow C \subset R \f$; ie a 
 * function of one or more variables that has always one-dimensional range, which for the 
 * specific CDF case, the image set \f$ C = [0,1]\f$. The CDF describes the probability that 
 * a real-valued random variable X with a given probability distribution will be found at a
 * value less than or equal to x. In the case of a continuous distribution, it gives the area
 * under the probability density function (PDF) from minus infinity to x.*/

template<class T>
class uqBaseScalarCdfClass {
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  /*! Instantiates an object of the class  given a prefix and the environment.*/
  uqBaseScalarCdfClass(const uqBaseEnvironmentClass& env, const char* prefix);
  
  //! Virtual destructor.
  virtual ~uqBaseScalarCdfClass();
  //@}

  //! @name Environment methods
  //@{
  //! Environment.  Access to private attribute m_env.
  const uqBaseEnvironmentClass&  env             () const;
  
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
  
  //! Writes the CDF of an allowed sub-environment to a file. 
  /*! This function does nothing and should \n not be called by the user.*/
  virtual void                    subWriteContents(const std::string&            varNamePrefix,
                                                   const std::string&            fileName,
                                                   const std::string&            fileType,
                                                   const std::set<unsigned int>& allowedSubEnvIds) const;
  //@}
protected:
  const uqBaseEnvironmentClass& m_env;
        std::string             m_prefix;
  mutable T                     m_minHorizontal;
  mutable T                     m_maxHorizontal;
};
// Default constructor -----------------------------
template<class T>
uqBaseScalarCdfClass<T>::uqBaseScalarCdfClass(
  const uqBaseEnvironmentClass& env,
  const char*                   prefix)
  :
  m_env          (env),
  m_prefix       ((std::string)(prefix)+""),
  m_minHorizontal(-INFINITY),
  m_maxHorizontal( INFINITY)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqBaseScalarCdfClass<T>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqBaseScalarCdfClass<T>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }
}
// Destructor ---------------------------------------
template<class T>
uqBaseScalarCdfClass<T>::~uqBaseScalarCdfClass()
{
}
// Environment methods ------------------------------
template <class T>
const uqBaseEnvironmentClass&
uqBaseScalarCdfClass<T>::env() const
{
  return m_env;
}
// --------------------------------------------------
template <class T>
const std::string&
uqBaseScalarCdfClass<T>::prefix() const
{
  return m_prefix;
}
// I/O methods---------------------------------------
template<class T>
void
uqBaseScalarCdfClass<T>::subWriteContents(
  const std::string&            varNamePrefix,
  const std::string&            fileName,
  const std::string&            fileType,
  const std::set<unsigned int>& allowedSubEnvIds) const
{
  std::cerr << "WARNING: uqBaseScalarCdfClass<T>::subWriteContents() being used..."
            << std::endl;
  return;
}

// --------------------------------------------------
// Operator defined outside class limits-------------
// --------------------------------------------------
//! Operator to be used with print().
template <class T>
std::ostream& operator<< (std::ostream& os, const uqBaseScalarCdfClass<T>& obj)
{
  obj.print(os);
  return os;
}

//*****************************************************
// Sampled cumulative distribution function class
//*****************************************************
/*! 
 * \class uqSampledScalarCdfClass
 * \brief A class for handling sampled CDFs.
 *
 * This class implements a sampled cumulative distribution function (CDF), given 
 * the grid points where it will be sampled and its resulting values.*/

template<class T>
class uqSampledScalarCdfClass : public uqBaseScalarCdfClass<T> {
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  /*! Instantiates an object of the class given a prefix, the environment, the grid points 
   * where it will be sampled/evaluated and its resulting values.*/
  uqSampledScalarCdfClass(const uqBaseEnvironmentClass& env,
                          const char*                   prefix,
                          const uqBaseOneDGridClass<T>& cdfGrid,
                          const std::vector<double>&    cdfValues);
  //! Destructor.
 ~uqSampledScalarCdfClass();
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
  using uqBaseScalarCdfClass<T>::m_env;
  using uqBaseScalarCdfClass<T>::m_prefix;
  using uqBaseScalarCdfClass<T>::m_minHorizontal;
  using uqBaseScalarCdfClass<T>::m_maxHorizontal;

  const uqBaseOneDGridClass<T>& m_cdfGrid;
  const std::vector<double>&    m_cdfValues;
//std::vector<double>&          m_sortedCdfValues;
};
// Default constructor -----------------------------
template<class T>
uqSampledScalarCdfClass<T>::uqSampledScalarCdfClass(
  const uqBaseEnvironmentClass& env,
  const char*                   prefix,
  const uqBaseOneDGridClass<T>& cdfGrid,
  const std::vector<double>&    cdfValues)
  :
  uqBaseScalarCdfClass<T>(env,((std::string)(prefix)+"").c_str()),
  m_cdfGrid              (cdfGrid  ),
  m_cdfValues            (cdfValues)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqSampledScalarCdfClass<T>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  //for (unsigned int i = 0; i < m_cdfValues.size(); ++i) {
  //  m_sortedCdfValues[i] = m_cdfValues[i];
  //}
  //std::sort(m_sortedCdfValues.begin(), m_sortedCdfValues.end());
 
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqSampledScalarCdfClass<T>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class T>
uqSampledScalarCdfClass<T>::~uqSampledScalarCdfClass()
{
}
// Math methods--------------------------------------
template<class T>
double
uqSampledScalarCdfClass<T>::value(T paramValue) const
{
  double result = 0.;
  if (paramValue <= m_cdfGrid[0]) {
    result = 0.;
  }
  else if (m_cdfGrid[m_cdfGrid.size()-1] <= paramValue) {
    result = 1.;
  }
  else {
    unsigned int intervalId = m_cdfGrid.findIntervalId(paramValue);
    UQ_FATAL_TEST_MACRO(intervalId >= (m_cdfGrid.size()-1),
                        m_env.worldRank(),
                        "uqSampledScalarCdfClass<T>::value()",
                        "invalid intervalId");

    double intervalLen = m_cdfGrid[intervalId+1] - m_cdfGrid[intervalId];
    double ratio = (paramValue - m_cdfGrid[intervalId])/intervalLen;
#if 0
    *m_env.subDisplayFile() << "In uqSampledScalarCdf::value()"
                            << ": paramValue = "              << paramValue
                            << ", intervalId = "              << intervalId
                            << ", cdfGrid.size() = "          << m_cdfGrid.size()
                            << ", m_cdfGrid[intervalId] = "   << m_cdfGrid[intervalId]
                            << ", m_cdfGrid[intervalId+1] = " << m_cdfGrid[intervalId+1]
                            << ", intervalLen = "             << intervalLen
                            << ", ratio = "                   << ratio
                            << std::endl;
#endif
    UQ_FATAL_TEST_MACRO(ratio < 0.,
                        m_env.worldRank(),
                        "uqSampledScalarCdfClass<T>::value()",
                        "invalid ratio");

    result = (1.-ratio)*m_cdfValues[intervalId] + ratio*m_cdfValues[intervalId+1];
  }

  return result;
}
//---------------------------------------------------
template<class T>
T
uqSampledScalarCdfClass<T>::inverse(double cdfValue) const
{
  //*m_env.subDisplayFile() << "In uqSampledScalarCdf::inverse(): cdfValue = " << cdfValue
  //                       << std::endl;
  UQ_FATAL_TEST_MACRO((cdfValue < 0.) || (1. < cdfValue),
                      m_env.worldRank(),
                      "uqSampledScalarCdfClass<T>::inverse()",
                      "invalid cdfValue");
  double result = 0.;
  unsigned int i = 0;
  unsigned int j = m_cdfValues.size()-1;
  bool searchPosition = true;
  do {
    if (cdfValue == m_cdfValues[i]) {
      while ((0 < i) && (cdfValue == m_cdfValues[i-1])) --i;
      result = m_cdfGrid[i];
      searchPosition = false;
    }
    
    if (cdfValue == m_cdfValues[j]) {
      while ((0 < j) && (cdfValue == m_cdfValues[j-1])) --j;
      result = m_cdfGrid[j];
      searchPosition = false;
    }

    if ((j-i) <= 0) {
      UQ_FATAL_TEST_MACRO(true,
                          m_env.worldRank(),
                          "uqSampledScalarCdfClass<T>::inverse()",
                          "invalid pair of values 'i' and 'j'");
    }
    else if ((j-i) == 1) {
      double ratio = (cdfValue-m_cdfValues[i])/(m_cdfValues[j]-m_cdfValues[i]);
      result = (1.-ratio)*m_cdfGrid[i] + ratio*m_cdfGrid[j];
      searchPosition = false;
    }
    else {
      unsigned int k= (unsigned int) ((i+j)*.5);
      if (cdfValue < m_cdfValues[k]) {
        j = k;
      }
      else if (cdfValue == m_cdfValues[k]) {
        while ((0 < k) && (cdfValue == m_cdfValues[k-1])) --k;
        result = m_cdfGrid[k];
        searchPosition = false;
      }
      else {
        i = k;
      }
    }
  } while (searchPosition);
  
  return result;
}
//---------------------------------------------------
template<class T>
void
uqSampledScalarCdfClass<T>::getSupport(T& minHorizontal, T& maxHorizontal) const
{
  if ((m_minHorizontal == -INFINITY) ||
      (m_maxHorizontal ==  INFINITY)) {
    UQ_FATAL_TEST_MACRO((m_minHorizontal != -INFINITY) || (m_maxHorizontal != INFINITY),
                        m_env.worldRank(),
                        "uqSampledScalarCdfClass<T>::getSupport()",
                        "unexpected values of m_minHorizontal and/or m_maxHorizontal");

    unsigned int iMax = m_cdfGrid.size();

    for (unsigned int i = 0; i < iMax; ++i) {
      if (m_cdfValues[i] > 0.) {
        if (i > 0) --i;
        m_minHorizontal = m_cdfGrid[i];
        break;
      }
    }

    UQ_FATAL_TEST_MACRO(m_minHorizontal == -INFINITY,
                        m_env.worldRank(),
                        "uqSampledScalarCdfClass<T>::getSupport()",
                        "unexpected value for m_minHorizontal");

    if (iMax == 1) {
      UQ_FATAL_TEST_MACRO(m_cdfValues[iMax - 1] != 1.,
                          m_env.worldRank(),
                          "uqSampledScalarCdfClass<T>::getSupport()",
                          "unexpected value for case 'iMax = 1'");
      m_maxHorizontal = m_cdfGrid[iMax-1];
    }
    else for (unsigned int i = 0; i < iMax; ++i) {
      if (m_cdfValues[iMax-1-i] < 1.) {
        if (i > 0) --i;
        m_maxHorizontal = m_cdfGrid[iMax-1-i];
        break;
      }
    }

    UQ_FATAL_TEST_MACRO(m_maxHorizontal == INFINITY,
                        m_env.worldRank(),
                        "uqSampledScalarCdfClass<T>::getSupport()",
                        "unexpected value for m_maxHorizontal");
  }

  minHorizontal = m_minHorizontal;
  maxHorizontal = m_maxHorizontal;

  return;
}
// I/O methods---------------------------------------
template <class T>
void
uqSampledScalarCdfClass<T>::print(std::ostream& os) const
{
  // Print values *of* grid points
  os << m_cdfGrid;

  // Print *cdf* values *at* grid points
  os << m_prefix << "values_sub" << m_env.subIdString() << " = zeros(" << m_cdfValues.size()
     << ","                                                            << 1
     << ");"
     << std::endl;
  os << m_prefix << "values_sub" << m_env.subIdString() << " = [";
  for (unsigned int j = 0; j < m_cdfValues.size(); ++j) {
    os << m_cdfValues[j] << " ";
  }
  os << "];"
     << std::endl;

  return;
}
//---------------------------------------------------
template<class T>
void
uqSampledScalarCdfClass<T>::subWriteContents(
  const std::string&            varNamePrefix,
  const std::string&            fileName,
  const std::string&            fileType,
  const std::set<unsigned int>& allowedSubEnvIds) const
{
  UQ_FATAL_TEST_MACRO(m_env.subRank() < 0,
                      m_env.worldRank(),
                      "uqSampledScalarCdfClass<T>::subWriteContents()",
                      "unexpected subRank");

  uqFilePtrSetStruct filePtrSet;
  if (m_env.openOutputFile(fileName,
                           fileType, // "m or hdf"
                           allowedSubEnvIds,
                           false,
                           filePtrSet)) {

    // Grid
    *filePtrSet.ofsVar << varNamePrefix << "grid_sub" << m_env.subIdString() << " = zeros(" << m_cdfGrid.size()
                       << ","                                                               << 1
                       << ");"
                       << std::endl;
    *filePtrSet.ofsVar << varNamePrefix << "grid_sub" << m_env.subIdString() << " = [";

    unsigned int savedPrecision = filePtrSet.ofsVar->precision();
    filePtrSet.ofsVar->precision(16);
    for (unsigned int j = 0; j < m_cdfGrid.size(); ++j) {
      *filePtrSet.ofsVar << m_cdfGrid[j] << " ";
    }
    filePtrSet.ofsVar->precision(savedPrecision);

    *filePtrSet.ofsVar << "];\n";

    // Values
    *filePtrSet.ofsVar << varNamePrefix << "values_sub" << m_env.subIdString() << " = zeros(" << m_cdfValues.size()
                       << ","                                                                 << 1
                       << ");"
                       << std::endl;
    *filePtrSet.ofsVar << varNamePrefix << "values_sub" << m_env.subIdString() << " = [";

    savedPrecision = filePtrSet.ofsVar->precision();
    filePtrSet.ofsVar->precision(16);
    for (unsigned int j = 0; j < m_cdfValues.size(); ++j) {
      *filePtrSet.ofsVar << m_cdfValues[j] << " ";
    }
    filePtrSet.ofsVar->precision(savedPrecision);

    *filePtrSet.ofsVar << "];\n";

    // Close file
    m_env.closeFile(filePtrSet,fileType);
  }

  return;
}

//*****************************************************
// Std cumulative distribution function class
//*****************************************************
/*! 
 * \class uqStdScalarCdfClass
 * \brief A class for handling standard CDFs.
 *
 * This class implements a standard cumulative distribution function (CDF). Its protected attribute
 * m_sampledCdfGrid is an object of the class uqSampledScalarCdfClass<T>, so all members of this
 * class are implemented using the members of uqSampledScalarCdfClass<T>.*/

template<class T>
class uqStdScalarCdfClass : public uqBaseScalarCdfClass<T> {
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  /*! Instantiates an object of the class given a prefix, the environment, the grid points 
   * where it will be evaluated and its resulting values.*/
  uqStdScalarCdfClass(const uqBaseEnvironmentClass& env,
                      const char*                   prefix,
                      const std::vector<T>&         cdfGrid,
                      const std::vector<double>&    cdfValues);
  //! Destructor
  ~uqStdScalarCdfClass();
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
  using uqBaseScalarCdfClass<T>::m_env;
  using uqBaseScalarCdfClass<T>::m_prefix;
  using uqBaseScalarCdfClass<T>::m_minHorizontal;
  using uqBaseScalarCdfClass<T>::m_maxHorizontal;

  const uqStdOneDGridClass<T> m_cdfGrid;
  const std::vector<double>   m_cdfValues;
  uqSampledScalarCdfClass<T>* m_sampledCdfGrid;
};
// Default constructor -----------------------------
template<class T>
uqStdScalarCdfClass<T>::uqStdScalarCdfClass(
  const uqBaseEnvironmentClass& env,
  const char*                   prefix,
  const std::vector<T>&         cdfGrid,
  const std::vector<double>&    cdfValues)
  :
  uqBaseScalarCdfClass<T>(env,((std::string)(prefix)+"").c_str()),
  m_cdfGrid              (env,prefix,cdfGrid),
  m_cdfValues            (cdfValues),
  m_sampledCdfGrid       (NULL)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqStdScalarCdfClass<T>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  m_sampledCdfGrid = new uqSampledScalarCdfClass<T>(env,prefix,m_cdfGrid,m_cdfValues);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqStdScalarCdfClass<T>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class T>
uqStdScalarCdfClass<T>::~uqStdScalarCdfClass()
{
  delete m_sampledCdfGrid;
}
// Math methods--------------------------------------
template<class T>
double
uqStdScalarCdfClass<T>::value(T paramValue) const
{
  return m_sampledCdfGrid->value(paramValue);
}
//---------------------------------------------------
template<class T>
T
uqStdScalarCdfClass<T>::inverse(double cdfValue) const
{
  return m_sampledCdfGrid->inverse(cdfValue);
}
//---------------------------------------------------
template<class T>
void
uqStdScalarCdfClass<T>::getSupport(T& minHorizontal, T& maxHorizontal) const
{
  return m_sampledCdfGrid->getSupport(minHorizontal,maxHorizontal);
}
// I/O methods---------------------------------------
template <class T>
void
uqStdScalarCdfClass<T>::print(std::ostream& os) const
{
  m_sampledCdfGrid->print(os);
  return;
}
//---------------------------------------------------
template<class T>
void
uqStdScalarCdfClass<T>::subWriteContents(
  const std::string&            varNamePrefix,
  const std::string&            fileName,
  const std::string&            fileType,
  const std::set<unsigned int>& allowedSubEnvIds) const
{
  m_sampledCdfGrid->subWriteContents(varNamePrefix,fileName,fileType,allowedSubEnvIds);
  return;
}
//---------------------------------------------------
// Method outside either class definition------------
//---------------------------------------------------

//*****************************************************
// Horizontal distance
//*****************************************************

//! It calculated the maximum horizontal distance between two CDFs.
template <class T>
double
horizontalDistance(const uqBaseScalarCdfClass<T>& cdf1,
                   const uqBaseScalarCdfClass<T>& cdf2,
                   double epsilon)
{
  double maxDistance     = 0.;
  double xForMaxDistance = 0.;

  double x1 = cdf1.inverse(epsilon*.5);
  double x2 = cdf1.inverse(1.-epsilon*.5);
  if (cdf1.env().subDisplayFile()) {
    *cdf1.env().subDisplayFile() << "In horizontalDistance()"
                                 << ", cdf1.prefix() = " << cdf1.prefix()
                                 << ", cdf2.prefix() = " << cdf2.prefix()
                                 << ", epsilon = "       << epsilon
                                 << ": x1 = "            << x1
                                 << ", x2 = "            << x2
                                 << std::endl;
  }

  //if (cdf1.env().subDisplayFile()) {
  //  *cdf1.env().subDisplayFile() << "In horizontalDistance: x1 = " << x1
  //                              << ", x2 = " << x2
  //                              << std::endl;
  //}

  double numEvaluationPoints = 1001.;
  for (double i = 0.; i < numEvaluationPoints; ++i) {
    double ratio = i/(numEvaluationPoints-1.); // IMPORTANT: Yes, '-1.'
    double x = (1.-ratio)*x1 + ratio*x2;
    double y = cdf2.inverse(cdf1.value(x));
    double d = fabs(x-y);
    if ((cdf1.env().subDisplayFile()) && (cdf1.env().displayVerbosity() >= 3)) {
      *cdf1.env().subDisplayFile() << "In horizontalDistance"
                                   << ": i = "                  << i
                                   << ", x = "                  << x
                                   << ", cdf1.value(x) = "      << cdf1.value(x)
                                   << ", y = "                  << y
                                   << ", d = "                  << d
                                   << ", currentMaxDistance = " << maxDistance
                                   << std::endl;
    }
    if (maxDistance < d) {
      maxDistance     = d;
      xForMaxDistance = x;
      if ((cdf1.env().subDisplayFile()) && (cdf1.env().displayVerbosity() >= 3)) {
        *cdf1.env().subDisplayFile() << "In horizontalDistance"
                                     << ": i = "               << i
                                     << ", NOW maxDistance = " << maxDistance
                                     << ", xForMaxDistance = " << xForMaxDistance
                                     << std::endl;
      }
    }
  }

  if (cdf1.env().subDisplayFile()) {
    *cdf1.env().subDisplayFile() << "In horizontalDistance()"
                                 << ", cdf1.prefix() = "   << cdf1.prefix()
                                 << ", cdf2.prefix() = "   << cdf2.prefix()
                                 << ", epsilon = "         << epsilon
                                 << ": maxDistance = "     << maxDistance
                                 << ", xForMaxDistance = " << xForMaxDistance
                                 << std::endl;
  }

  return maxDistance;
}

}  // End namespace QUESO

#endif // __UQ_SCALAR_CUMULATIVE_DISTRIBUTION_FUNCTION_H__
