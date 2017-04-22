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

#ifndef UQ_CHAIN_POSITION_DATA_H
#define UQ_CHAIN_POSITION_DATA_H

#include <queso/Environment.h>

namespace QUESO {

class GslVector;

/*! \file MarkovChainPositionData.h
 * \brief A templated class that represents a Markov Chain.
 *
 * \class MarkovChainPositionData
 * \brief A templated class that represents a Markov Chain.
 *
 * This class implements a Markov Chain. It contains important information about a position such as
 * whether or not it is out of the target support and its logarithmic values of the likelihood and
 * of the target, besides the position location. A Markov chain is collection of random variables
 * {X_t}, where the index t runs through 0, 1, ..., having the property that, given the present,
 * the future is conditionally independent of the past. For instance, a Markov chain is passed to
 * the Metropolis-Hastings algorithm and on it is calculated the acceptance ration
 * MetropolisHastingsSG::alpha(). */

template <class V = GslVector>
class MarkovChainPositionData
{
public:
    //! @name Constructor/Destructor methods
  //@{
  //! Constructor 1.
  /*! It allocates a new Markov chain, given the environment. All the other private attributes are either
   * NULL, false or zero.*/
  MarkovChainPositionData(const BaseEnvironment& env);

  //! Constructor 2.
 /*! It allocates a new Markov chain, and the vector \c vecValues, given the environment, the likelihood and
  * target values and sets whether or not it is out of target support.*/
  MarkovChainPositionData(const BaseEnvironment& env,
                                 const V& vecValues,
                                 bool     outOfTargetSupport,
                                 double   logLikelihood,
                                 double   logTarget);
  //! Constructor 3: copy.
  /*! The new Markov chain is a copy of \c rhs.*/
  MarkovChainPositionData(const MarkovChainPositionData<V>& rhs);

  //! Destructor
  ~MarkovChainPositionData();
  //@}

  //! @name Set methods
  //@{
  //! Assignment operator.
  MarkovChainPositionData<V>& operator= (const MarkovChainPositionData<V>& rhs);
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
  friend std::ostream& operator<<(std::ostream& os,
      const MarkovChainPositionData<V>& obj)
  {
    obj.print(os);

    return os;
  }
  //@}

private:
  const BaseEnvironment& m_env;
  V*     m_vecValues;
  bool   m_outOfTargetSupport;
  double m_logLikelihood;
  double m_logTarget;
};

}  // End namespace QUESO

#endif // UQ_CHAIN_POSITION_DATA_H
