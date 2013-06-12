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

#ifndef __UQ_VECTOR_SET_H__
#define __UQ_VECTOR_SET_H__

#include <uqEnvironment.h>
#include <uqDefines.h>

/*! \file uqVectorSetClass.h
 * \brief A templated class for handling sets.
 * 
 * \class uqVectorSetClass
 * \brief A templated class for handling sets.
 *
 * This class allows the mathematical definition of a scalar function such as:
 * \f$ \pi: B \subset R^n \rightarrow R \f$, since it requires the specification 
 * of the domain \f$ B \f$, which is a subset of the vector space \f$ R^n \f$,  
 * which is itself a set.*/


template <class V, class M>
class uqVectorSpaceClass;

template <class V, class M>
class uqVectorSetClass
{
public:
  //! @name Constructor/Destructor methods.
  //@{ 
  //! Default Constructor
  /*! It should not be used by the user.*/
  uqVectorSetClass();
  
  //! Shaped constructor.
  /*! Creates a vector set given an environment, a identifying prefix and a volume.*/
  uqVectorSetClass(const uqBaseEnvironmentClass& env, const char* prefix, double volume);
  
  //! Virtual destructor.
  virtual ~uqVectorSetClass();
  //@}
  
  //! @name Environment methods
  //@{
  //! Environment.  Access to private attribute m_env.
  const uqBaseEnvironmentClass&  env        ()                 const;
  
  //! Access to private attribute m_prefix.
  const std::string&             prefix     ()                 const;
  //@}
  
  //! @name Mathematical methods.
  //@{
  //! Set volume; access to private attribute m_volume.
  double                   volume     ()                 const;
  
  //! Vector space to which \c this set belongs to. See template specialization.
  virtual const uqVectorSpaceClass<V,M>& vectorSpace()                 const = 0;
  
  //! Checks whether a set contains vector \c vec. See template specialization.
  virtual       bool                     contains   (const V& vec)     const = 0;
  //@}
  
  //! @name I/O methods.
  //@{
  //! Prints nothing.
  virtual       void                     print      (std::ostream& os) const; 
  //@}
  
protected:
  const uqBaseEnvironmentClass& m_env;
        std::string             m_prefix;
        double                  m_volume;
};

// --------------------------------------------------
// Constructor/Destructor methods -------------------
// Default constructor-------------------------------
template <class V, class M>
uqVectorSetClass<V,M>::uqVectorSetClass()
  :
  m_env(*(new uqEmptyEnvironmentClass()))
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.worldRank(),
                      "uqVectorSetClass<V,M>::constructor(), default",
                      "should not be used by user");
}
// Shaped constructor--------------------------------
template <class V, class M>
uqVectorSetClass<V,M>::uqVectorSetClass(
  const uqBaseEnvironmentClass& env,
  const char*                   prefix,
        double                  volume)
  :
  m_env   (env),
  m_prefix(prefix),
  m_volume(volume)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqVectorSetClass<V,M>::constructor()"
                           << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqVectorSetClass<V,M>::constructor()"
                           << std::endl;
  }
}
// Destructor --------------------------------------------
template <class V, class M>
uqVectorSetClass<V,M>::~uqVectorSetClass()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqVectorSetClass<V,M>::destructor()"
                           << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqVectorSetClass<V,M>::destructor()"
                           << std::endl;
  }
}

// Environment methods------------------------------------
// -------------------------------------------------------
template <class V, class M>
const uqBaseEnvironmentClass&
uqVectorSetClass<V,M>::env() const
{
  return m_env;
}
// -------------------------------------------------------
template <class V, class M>
const std::string&
uqVectorSetClass<V,M>::prefix() const
{
  return m_prefix;
}
// Mathematical methods-----------------------------------
template <class V, class M>
double
uqVectorSetClass<V,M>::volume() const
{
  return m_volume;
}
// I/O methods--------------------------------------------
template <class V, class M>
void
uqVectorSetClass<V,M>::print(std::ostream& os) const
{
  os << "In uqVectorSetClass<V,M>::print()"
     << ": nothing to be printed" << std::endl;
  return;
}
// --------------------------------------------------
template<class V, class M>
std::ostream&
operator<<(std::ostream& os, const uqVectorSetClass<V,M>& obj)
{
  obj.print(os);

  return os;
}

#endif // __UQ_VECTOR_SET_H__

