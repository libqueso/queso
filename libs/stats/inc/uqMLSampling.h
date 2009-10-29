/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * $Id$
 *
 * Brief description of this file: 
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __UQ_MULTI_LEVEL_SAMPLING_H__
#define __UQ_MULTI_LEVEL_SAMPLING_H__

#include <uqMLSamplingOptions.h>
#include <uqFiniteDistribution.h>
#include <uqVectorRV.h>
#include <uqVectorSpace.h>
#include <uqMarkovChainPositionData.h>
#include <uqScalarFunctionSynchronizer.h>
#include <uqSequenceOfVectors.h>
#include <uqArrayOfSequences.h>
#include <sys/time.h>
#include <fstream>

#define UQ_ML_SAMPLING_USES_OMEGA_LN

struct uqLinkedChainControlStruct
{
  unsigned int initialPositionIndexInPreviousChain;
  unsigned int numberOfPositions;
};

struct uqLinkedChainsPerNodeStruct
{
  std::vector<uqLinkedChainControlStruct> linkedChains;
};

/*! A templated class that represents a Multi Level sampler.
 */
template <class P_V,class P_M>
class uqMLSamplingClass
{
public:
  /*! Constructor: */
  uqMLSamplingClass(/*! Prefix                  */ const char*                               prefix,
                    /*! The prior rv            */ const uqBaseVectorRVClass      <P_V,P_M>& priorRv,
                    /*! The likelihood function */ const uqBaseScalarFunctionClass<P_V,P_M>& likelihoodFunction);
  /*! Destructor: */
 ~uqMLSamplingClass();

  /*! Operation to generate the chain */
  void   generateSequence   (uqBaseVectorSequenceClass<P_V,P_M>&             workingChain,
                             uqScalarSequenceClass<double>*                  workingLogLikelihoodValues,
                             uqScalarSequenceClass<double>*                  workingLogTargetValues);

  void   print              (std::ostream& os) const;

private:
  const uqBaseEnvironmentClass&             m_env;
  const uqBaseVectorRVClass      <P_V,P_M>& m_priorRv;
  const uqBaseScalarFunctionClass<P_V,P_M>& m_likelihoodFunction;
  const uqVectorSpaceClass       <P_V,P_M>& m_vectorSpace;
};

template<class P_V,class P_M>
std::ostream& operator<<(std::ostream& os, const uqMLSamplingClass<P_V,P_M>& obj);

template<class P_V,class P_M>
uqMLSamplingClass<P_V,P_M>::uqMLSamplingClass(
  const char*                               prefix,
  const uqBaseVectorRVClass      <P_V,P_M>& priorRv,            
  const uqBaseScalarFunctionClass<P_V,P_M>& likelihoodFunction)
  :
  m_env               (priorRv.env()),
  m_priorRv           (priorRv),
  m_likelihoodFunction(likelihoodFunction),
  m_vectorSpace       (m_priorRv.imageSet().vectorSpace()),
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering uqMLSamplingClass<P_V,P_M>::constructor()"
                            << std::endl;
  }

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving uqMLSamplingClass<P_V,P_M>::constructor()"
                            << std::endl;
  }
}

template<class P_V,class P_M>
uqMLSamplingClass<P_V,P_M>::~uqMLSamplingClass()
{
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::generateSequence(
  uqBaseVectorSequenceClass<P_V,P_M>& workingChain,
  uqScalarSequenceClass<double>*      workingLogLikelihoodValues,
  uqScalarSequenceClass<double>*      workingLogTargetValues)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Entering uqMLSamplingClass<P_V,P_M>::generateSequence()..."
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving uqMLSamplingClass<P_V,P_M>::generateSequence()"
                            << std::endl;
  }

  return;
}

template<class P_V,class P_M>
std::ostream& operator<<(std::ostream& os, const uqMLSamplingClass<P_V,P_M>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_MULTI_LEVEL_SAMPLING_H__
