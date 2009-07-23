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
#include <uqVectorRV.h>
#include <uqVectorSpace.h>
#include <uqMarkovChainPositionData.h>
#include <uqScalarFunctionSynchronizer.h>
#include <uqSequenceOfVectors.h>
#include <uqArrayOfSequences.h>
#include <sys/time.h>
#include <fstream>

/*! A templated class that represents a Markov Chain generator. 'SG' stands for 'Sequence Generator'.
 */
template <class P_V,class P_M>
class uqMLSamplingClass
{
public:

  /*! Constructor: */
  uqMLSamplingClass(/*! Prefix                 */ const char*                         prefix,                  
                    /*! The source rv          */ const uqBaseVectorRVClass<P_V,P_M>& sourceRv,                
                    /*! Initial chain position */ const P_V&                          initialPosition,
                    /*! Proposal cov. matrix   */ const P_M*                          inputProposalCovMatrix);  
  /*! Destructor: */
 ~uqMLSamplingClass();

  /*! Operation to generate the chain */
  void   generateSequence(uqBaseVectorSequenceClass<P_V,P_M>& workingChain); /*! Generate the chain and store it in 'workingChain' object */

  void   print           (std::ostream& os) const;

private:
  const uqBaseEnvironmentClass&              m_env;
  const uqVectorSpaceClass <P_V,P_M>&        m_vectorSpace;
  const uqBaseJointPdfClass<P_V,P_M>&        m_targetPdf;
        P_V                                  m_initialPosition;
  const P_M*                                 m_initialProposalCovMatrix;
        bool                                 m_nullInputProposalCovMatrix;

  const uqBaseVectorRVClass       <P_V,P_M>& m_sourceRv;
        uqBaseJointPdfClass       <P_V,P_M>* m_currPdf;
        uqGenericVectorRVClass    <P_V,P_M>* m_currRv;
        uqMarkovChainSGClass     <P_V,P_M>*  m_mcSeqGenerator;
        uqBaseVectorSequenceClass<P_V,P_M>*  m_currChain;
        uqBaseVectorSequenceClass<P_V,P_M>*  m_prevChain;

        uqMLSamplingOptionsClass             m_options;
};

template<class P_V,class P_M>
std::ostream& operator<<(std::ostream& os, const uqMLSamplingClass<P_V,P_M>& obj);

template<class P_V,class P_M>
uqMLSamplingClass<P_V,P_M>::uqMLSamplingClass(
  const char*                         prefix,
  const uqBaseVectorRVClass<P_V,P_M>& sourceRv,
  const P_V&                          initialPosition,
  const P_M*                          inputProposalCovMatrix)
  :
  m_env                       (sourceRv.env()),
  m_vectorSpace               (sourceRv.imageSet().vectorSpace()),
  m_targetPdf                 (sourceRv.pdf()),
  m_initialPosition           (initialPosition),
  m_initialProposalCovMatrix  (inputProposalCovMatrix),
  m_nullInputProposalCovMatrix(inputProposalCovMatrix == NULL),
  m_sourceRv                  (sourceRv),
  m_currPdf                   (NULL),
  m_currRv                    (NULL),
  m_mcSeqGenerator            (NULL),
  m_currChain                 (NULL),
  m_prevChain                 (NULL),
  m_options                   (m_env,prefix)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering uqMLSamplingClass<P_V,P_M>::constructor()"
                            << std::endl;
  }

  m_options.scanOptionsValues();

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving uqMLSamplingClass<P_V,P_M>::constructor()"
                           << std::endl;
  }
}

template<class P_V,class P_M>
uqMLSamplingClass<P_V,P_M>::~uqMLSamplingClass()
{
  if (m_currChain) {
    m_currChain->clear();
    delete m_currChain;
  }
  if (m_prevChain) {
    m_prevChain->clear();
    delete m_prevChain;
  }
  if (m_mcSeqGenerator            ) delete m_mcSeqGenerator;
  if (m_currRv                    ) delete m_currRv;
  if (m_currPdf                   ) delete m_currPdf;
  if (m_nullInputProposalCovMatrix) delete m_initialProposalCovMatrix;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::generateSequence(uqBaseVectorSequenceClass<P_V,P_M>& workingChain)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqMLSamplingClass<P_V,P_M>::generateSequence()..."
                           << std::endl;
  }

  m_currChain = new uqSequenceOfVectorsClass<P_V,P_M>(m_sourceRv.imageSet().vectorSpace(),0,m_options.m_prefix+"curr_chain");
  m_prevChain = new uqSequenceOfVectorsClass<P_V,P_M>(m_sourceRv.imageSet().vectorSpace(),0,m_options.m_prefix+"prev_chain");

  for (unsigned int i = 0; i < m_options.m_maxNumberOfLevels; ++i) {
    if (i > 0) {
      m_prevChain->clear();
      *(m_prevChain) = *(m_currChain);
      m_currChain->clear();
      delete m_currRv;
      delete m_currPdf;
    }
#if 0
    m_currPdf = new uqPoweredPdfClass<P_V,P_M>(m_options.m_prefix.c_str(),
                                               m_sourceRv.pdf(),
                                               .5);

    m_currRv = new uqGenericRvClass<P_V,P_M>();
#endif
    m_currRv->setPdf(*m_currPdf);

    m_mcSeqGenerator = new uqMarkovChainSGClass<P_V,P_M>(m_options.m_prefix.c_str(),
                                                         *m_currRv,
                                                         m_initialPosition,
                                                         m_initialProposalCovMatrix);

    m_mcSeqGenerator->generateSequence(*m_currChain);
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
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