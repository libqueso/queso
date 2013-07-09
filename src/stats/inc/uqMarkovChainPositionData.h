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

#ifndef __UQ_CHAIN_POSITION_DATA_H__
#define __UQ_CHAIN_POSITION_DATA_H__

#include <uqEnvironment.h>

/*! \file uqMarkovChainPositionData.h
 * \brief A templated class that represents a Markov Chain.
 *
 * \class uqMarkovChainPositionDataClass
 * \brief A templated class that represents a Markov Chain.
 *
 * This class implements a Markov Chain. It contains important information about a position such as 
 * whether or not it is out of the target support and its logarithmic values of the likelihood and 
 * of the target, besides the position location. A Markov chain is collection of random variables 
 * {X_t}, where the index t runs through 0, 1, ..., having the property that, given the present, 
 * the future is conditionally independent of the past. For instance, a Markov chain is passed to 
 * the Metropolis-Hastings algorithm and on it is calculated the acceptance ration 
 * uqMetropolisHastingsSGClass::alpha(). */

template <class V>
class uqMarkovChainPositionDataClass
{
public:
    //! @name Constructor/Destructor methods
  //@{
  //! Constructor 1.
  /*! It allocates a new Markov chain, given the environment. All the other private attributes are either
   * NULL, false or zero.*/
  uqMarkovChainPositionDataClass(const uqBaseEnvironmentClass& env);
  
  //! Constructor 2.
 /*! It allocates a new Markov chain, and the vector \c vecValues, given the environment, the likelihood and 
  * target values and sets whether or not it is out of target support.*/
  uqMarkovChainPositionDataClass(const uqBaseEnvironmentClass& env,
                                 const V& vecValues,
                                 bool     outOfTargetSupport,
                                 double   logLikelihood,
                                 double   logTarget);
  //! Constructor 3: copy.
  /*! The new Markov chain is a copy of \c rhs.*/
  uqMarkovChainPositionDataClass(const uqMarkovChainPositionDataClass<V>& rhs);
  
  //! Destructor
  ~uqMarkovChainPositionDataClass();
  //@}

  //! @name Set methods
  //@{
  //! Assignment operator.
  uqMarkovChainPositionDataClass<V>& operator= (const uqMarkovChainPositionDataClass<V>& rhs);
  //@}
  
  //! @name Statistical/Mathematical methods
  //@{
  //! Values of the chain (vector); access to private attribute m_vecValues.
  const V& vecValues         () const;
  
  //! Whether or not a position is out of target support; access to private attribute m_outOfTargetSupport.
  bool     outOfTargetSupport() const;
  
  //! Logarithm of the value of the likelihood; access to private attribute m_logLikelihood.
  double   logLikelihood     () const;
  
  //! Logarithm of the value of the target; access to private attribute m_logTarget.
  double   logTarget         () const;

 //! Sets the values of the chain.
  void     set               (const V& vecValues,
                              bool     outOfTargetSupport,
                              double   logLikelihood,
                              double   logTarget);
  //@}
  //! @name I/O methods
  //@{
  //! TODO: Prints the Markov chain.
  /*! \todo: implement me!*/
  void     print             (std::ostream& os) const;
  //@}
  
private:
  const uqBaseEnvironmentClass& m_env;
  V*     m_vecValues;
  bool   m_outOfTargetSupport;
  double m_logLikelihood;
  double m_logTarget;
};

// Constructor 1 -----------------------------------
template <class V>
uqMarkovChainPositionDataClass<V>::uqMarkovChainPositionDataClass(const uqBaseEnvironmentClass& env)
  :
  m_env               (env),
  m_vecValues         (NULL),
  m_outOfTargetSupport(false),
  m_logLikelihood     (0.),
  m_logTarget         (0.)
{
}
// Constructor 2 -----------------------------------
template <class V>
uqMarkovChainPositionDataClass<V>::uqMarkovChainPositionDataClass(
  const uqBaseEnvironmentClass& env,
  const V& vecValues,
  bool     outOfTargetSupport,
  double   logLikelihood,
  double   logTarget)
  :
  m_env               (env),
  m_vecValues         (new V(vecValues)),
  m_outOfTargetSupport(outOfTargetSupport),
  m_logLikelihood     (logLikelihood),
  m_logTarget         (logTarget)
{
}
// Copy constructor---------------------------------
template <class V>
uqMarkovChainPositionDataClass<V>::uqMarkovChainPositionDataClass(const uqMarkovChainPositionDataClass<V>& rhs)
  :
  m_env               (rhs.m_env               ),
  m_vecValues         (new V(*rhs.m_vecValues )),
  m_outOfTargetSupport(rhs.m_outOfTargetSupport),
  m_logLikelihood     (rhs.m_logLikelihood     ),
  m_logTarget         (rhs.m_logTarget         )
{
}
// Destructor---------------------------------------
template <class V>
uqMarkovChainPositionDataClass<V>::~uqMarkovChainPositionDataClass()
{
  if (m_vecValues) delete m_vecValues;
}
// Set methods--------------------------------------
template <class V>
uqMarkovChainPositionDataClass<V>&
uqMarkovChainPositionDataClass<V>::operator=(const uqMarkovChainPositionDataClass<V>& rhs)
{
  if (m_vecValues == NULL) m_vecValues = new V(*rhs.m_vecValues);
  else                    *m_vecValues = *rhs.m_vecValues;
  m_outOfTargetSupport = rhs.m_outOfTargetSupport;
  m_logLikelihood      = rhs.m_logLikelihood;
  m_logTarget          = rhs.m_logTarget;

  return *this;
}

// Statistical methods-------------------------------
template <class V>
const V&
uqMarkovChainPositionDataClass<V>::vecValues() const
{
  UQ_FATAL_TEST_MACRO((m_vecValues == NULL),
                      m_env.worldRank(),
                      "uqMarkovChainPositionDataClass<V>::vecValues()",
                      "m_vecValues is NULL");
  return *m_vecValues;
}

//--------------------------------------------------
template <class V>
bool
uqMarkovChainPositionDataClass<V>::outOfTargetSupport() const
{
  return m_outOfTargetSupport;
}

//--------------------------------------------------
template <class V>
double
uqMarkovChainPositionDataClass<V>::logLikelihood() const
{
  return m_logLikelihood;
}

//--------------------------------------------------
template <class V>
double
uqMarkovChainPositionDataClass<V>::logTarget() const
{
  return m_logTarget;
}
//--------------------------------------------------
template <class V>
void
uqMarkovChainPositionDataClass<V>::set(
  const V& vecValues,
  bool     outOfTargetSupport,
  double   logLikelihood,
  double   logTarget)
{
  if (m_vecValues == NULL) m_vecValues = new V(vecValues);
  else                    *m_vecValues = vecValues;
  m_outOfTargetSupport = outOfTargetSupport;
  m_logLikelihood      = logLikelihood;
  m_logTarget          = logTarget;

  return;
}
//--------------------------------------------------
// Operator declared outside class definition-------
//--------------------------------------------------
template <class V>
std::ostream& operator<<(std::ostream& os, const uqMarkovChainPositionDataClass<V>& obj)
{
  obj.print(os);

  return os;
}

#endif // __UQ_CHAIN_POSITION_DATA_H__
