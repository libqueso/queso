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

#ifndef __UQ_SIP_H__
#define __UQ_SIP_H__

#include <uqStatisticalInverseProblemOptions.h>
#include <uqMarkovChainSG1.h>
#include <uqMLSampling.h>
#include <uqInstantiateIntersection.h>
#include <uqVectorRV.h>
#include <uqScalarFunction.h>
#include <hpct.h>

/*! A templated class that represents statistical inverse problems. */
/*! */
/*! Conceptually, a statistical inverse problem has two input entities and one output entity.
    The input entities are the prior rv and the likelihood function. The output entity is the posterior rv, which stores the solution according the Bayesian approach.
    A similar situation occurs e.g. in the case of a system Ax=b of linear equations, where A and b are inputs, and x is the solution of the inverse problem. */
/*! The solution of a statistical inverse problem is computed by calling 'solveWithBayesMarkovChain(...)'.
    Upon return from such operation, the posterior rv is available through the operation 'postRv()'. Such posterior rv is able to supply a joint pdf (up to a multiplicative constant) and a vector realizer. */
template <class P_V,class P_M>
class uqStatisticalInverseProblemClass
{
public:
  
  /*! Constructor: */
  uqStatisticalInverseProblemClass(/*! The prefix              */ const char*                               prefix,             
                                   /*! The prior rv            */ const uqBaseVectorRVClass      <P_V,P_M>& priorRv,            
                                   /*! The likelihood function */ const uqBaseScalarFunctionClass<P_V,P_M>& likelihoodFunction, 
                                   /*! The posterior rv        */       uqGenericVectorRVClass   <P_V,P_M>& postRv);
  /*! Destructor: */
 ~uqStatisticalInverseProblemClass();

        bool computeSolutionFlag      () const;
	/*! Operation to solve the problem */
        void solveWithBayesMarkovChain(const P_V& initialValues,
                                       const P_M* initialProposalCovMatrix);
        void solveWithBayesMLSampling (const P_V& initialValues,
                                       const P_M* initialProposalCovMatrix);
  const uqBaseVectorRVClass   <P_V,P_M>& priorRv() const;
  const uqGenericVectorRVClass<P_V,P_M>& postRv () const;

        void print                    (std::ostream& os) const;

private:
  const uqBaseEnvironmentClass&                 m_env;

  const uqBaseVectorRVClass       <P_V,P_M>&    m_priorRv;
  const uqBaseScalarFunctionClass <P_V,P_M>&    m_likelihoodFunction;
        uqGenericVectorRVClass    <P_V,P_M>&    m_postRv;

        uqVectorSetClass          <P_V,P_M>*    m_solutionDomain;
        uqBaseJointPdfClass       <P_V,P_M>*    m_solutionPdf;
        uqBaseVectorMdfClass      <P_V,P_M>*    m_subSolutionMdf;
        uqBaseVectorCdfClass      <P_V,P_M>*    m_subSolutionCdf;
        uqBaseVectorRealizerClass <P_V,P_M>*    m_solutionRealizer;

        uqMarkovChainSGClass      <P_V,P_M>*    m_mcSeqGenerator;
        uqMLSamplingClass         <P_V,P_M>*    m_mlSampler;
        uqBaseVectorSequenceClass <P_V,P_M>*    m_chain;
        uqArrayOfOneDGridsClass   <P_V,P_M>*    m_subMdfGrids;
        uqArrayOfOneDTablesClass  <P_V,P_M>*    m_subMdfValues;

        uqStatisticalInverseProblemOptionsClass m_options;
};

template<class P_V,class P_M>
std::ostream& operator<<(std::ostream& os, const uqStatisticalInverseProblemClass<P_V,P_M>& obj);

template <class P_V,class P_M>
uqStatisticalInverseProblemClass<P_V,P_M>::uqStatisticalInverseProblemClass(
  const char*                               prefix,
  const uqBaseVectorRVClass      <P_V,P_M>& priorRv,
  const uqBaseScalarFunctionClass<P_V,P_M>& likelihoodFunction,
        uqGenericVectorRVClass   <P_V,P_M>& postRv)
  :
  m_env               (priorRv.env()),
  m_priorRv           (priorRv),
  m_likelihoodFunction(likelihoodFunction),
  m_postRv            (postRv),
  m_solutionDomain    (NULL),
  m_solutionPdf       (NULL),
  m_subSolutionMdf    (NULL),
  m_subSolutionCdf    (NULL),
  m_solutionRealizer  (NULL),
  m_mcSeqGenerator    (NULL),
  m_mlSampler         (NULL),
  m_chain             (NULL),
  m_options           (m_env,prefix)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering uqStatisticalInverseProblemClass<P_V,P_M>::constructor()"
                           << ": prefix = " << m_options.m_prefix
                           << std::endl;
  }

  m_options.scanOptionsValues();

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving uqStatisticalInverseProblemClass<P_V,P_M>::constructor()"
                           << ": prefix = " << m_options.m_prefix
                           << std::endl;
  }

  return;
}

template <class P_V,class P_M>
uqStatisticalInverseProblemClass<P_V,P_M>::~uqStatisticalInverseProblemClass()
{
  if (m_chain) {
    m_chain->clear();
    delete m_chain;
  }
  if (m_mlSampler       ) delete m_mlSampler;
  if (m_mcSeqGenerator  ) delete m_mcSeqGenerator;
  if (m_solutionRealizer) delete m_solutionRealizer;
  if (m_subSolutionCdf  ) delete m_subSolutionCdf;
  if (m_subSolutionMdf  ) delete m_subSolutionMdf;
  if (m_solutionPdf     ) delete m_solutionPdf;
  if (m_solutionDomain  ) delete m_solutionDomain;
}

template <class P_V,class P_M>
void
uqStatisticalInverseProblemClass<P_V,P_M>::solveWithBayesMarkovChain(
  const P_V& initialValues,
  const P_M* initialProposalCovMatrix)
{
  //hpct_timer_begin("BayesMarkovChain");     TODO: revisit timing output

  m_env.fullComm().Barrier();
  m_env.syncPrintDebugMsg("Entering uqStatisticalInverseProblemClass<P_V,P_M>::solveWithBayesMarkovChain()",1,3000000,m_env.fullComm());

  if (m_options.m_computeSolution == false) {
    if ((m_env.subDisplayFile())) {
      *m_env.subDisplayFile() << "In uqStatisticalInverseProblemClass<P_V,P_M>::solveWithBayesMarkovChain()"
                              << ": avoiding solution, as requested by user"
                              << std::endl;
    }
    return;
  }
  if ((m_env.subDisplayFile())) {
    *m_env.subDisplayFile() << "In uqStatisticalInverseProblemClass<P_V,P_M>::solveWithBayesMarkovChain()"
                            << ": computing solution, as requested by user"
                            << std::endl;
  }

  if (m_mlSampler       ) delete m_mlSampler;
  if (m_mcSeqGenerator  ) delete m_mcSeqGenerator;
  if (m_solutionRealizer) delete m_solutionRealizer;
  if (m_subSolutionCdf  ) delete m_subSolutionCdf;
  if (m_subSolutionMdf  ) delete m_subSolutionMdf;
  if (m_solutionPdf     ) delete m_solutionPdf;
  if (m_solutionDomain  ) delete m_solutionDomain;

  P_V numEvaluationPointsVec(m_priorRv.imageSet().vectorSpace().zeroVector());
  numEvaluationPointsVec.cwSet(250.);

  // Compute output pdf up to a multiplicative constant: Bayesian approach
  m_solutionDomain = uqInstantiateIntersection(m_priorRv.pdf().domainSet(),m_likelihoodFunction.domainSet());



  m_solutionPdf = new uqBayesianJointPdfClass<P_V,P_M>(m_options.m_prefix.c_str(),
                                                       m_priorRv.pdf(),
                                                       m_likelihoodFunction,
                                                      *m_solutionDomain);


  m_postRv.setPdf(*m_solutionPdf);

  // Compute output realizer: Markov Chain approach
  m_chain = new uqSequenceOfVectorsClass<P_V,P_M>(m_postRv.imageSet().vectorSpace(),0,m_options.m_prefix+"chain");
  m_mcSeqGenerator = new uqMarkovChainSGClass<P_V,P_M>(m_options.m_prefix.c_str(),
                                                       m_postRv,
                                                       initialValues,
                                                       initialProposalCovMatrix);

  m_mcSeqGenerator->generateSequence(*m_chain);

  m_solutionRealizer = new uqSequentialVectorRealizerClass<P_V,P_M>(m_options.m_prefix.c_str(),
                                                                   *m_chain);

  m_postRv.setRealizer(*m_solutionRealizer);

  m_env.syncPrintDebugMsg("In uqStatisticalInverseProblemClass<P_V,P_M>::solveWithBayesMarkovChain(), code place 1",3,3000000,m_env.fullComm());
  //m_env.fullComm().Barrier();

  // Compute output mdf: uniform sampling approach
  m_subMdfGrids  = new uqArrayOfOneDGridsClass <P_V,P_M>((m_options.m_prefix+"Mdf_").c_str(),m_postRv.imageSet().vectorSpace());
  m_subMdfValues = new uqArrayOfOneDTablesClass<P_V,P_M>((m_options.m_prefix+"Mdf_").c_str(),m_postRv.imageSet().vectorSpace());
  m_chain->subUniformlySampledMdf(numEvaluationPointsVec, // input
                                  *m_subMdfGrids,         // output
                                  *m_subMdfValues);       // output
  m_subSolutionMdf = new uqSampledVectorMdfClass<P_V,P_M>(m_options.m_prefix.c_str(),
                                                          *m_subMdfGrids,
                                                          *m_subMdfValues);
  m_postRv.setMdf(*m_subSolutionMdf);

  if ((m_options.m_dataOutputFileName                       != UQ_SIP_FILENAME_FOR_NO_FILE) &&
      (m_options.m_dataOutputAllowedSet.find(m_env.subId()) != m_options.m_dataOutputAllowedSet.end()       )) {
    if (m_env.subRank() == 0) {
      // Write data output file
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "Opening data output file '" << m_options.m_dataOutputFileName
                                      << "' for calibration problem with problem with prefix = " << m_options.m_prefix
                                      << std::endl;
      }

      // Open file
#if 0
      // Always write over an eventual pre-existing file
      std::ofstream* ofsvar = new std::ofstream((m_options.m_dataOutputFileName+"_sub"+m_env.subIdString()+".m").c_str(), std::ofstream::out | std::ofstream::trunc);
#else
      // Always write at the end of an eventual pre-existing file
      std::ofstream* ofsvar = new std::ofstream((m_options.m_dataOutputFileName+"_sub"+m_env.subIdString()+".m").c_str(), std::ofstream::out | std::ofstream::in | std::ofstream::ate);
      if ((ofsvar            == NULL ) ||
          (ofsvar->is_open() == false)) {
        delete ofsvar;
        ofsvar = new std::ofstream((m_options.m_dataOutputFileName+"_sub"+m_env.subIdString()+".m").c_str(), std::ofstream::out | std::ofstream::trunc);
      }
#endif
      UQ_FATAL_TEST_MACRO((ofsvar && ofsvar->is_open()) == false,
                          m_env.fullRank(),
                          "uqStatisticalInverseProblem<P_V,P_M>::solveWithBayesMarkovChain()",
                          "failed to open file");

      m_postRv.mdf().print(*ofsvar);

      // Close file
      ofsvar->close();
      delete ofsvar;
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "Closed data output file '" << m_options.m_dataOutputFileName
                                      << "' for calibration problem with problem with prefix = " << m_options.m_prefix
                                      << std::endl;
      }
    }
  }
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << std::endl;
  }

  m_env.syncPrintDebugMsg("Leaving uqStatisticalInverseProblemClass<P_V,P_M>::solveWithBayesMarkovChain()",1,3000000,m_env.fullComm());
  m_env.fullComm().Barrier();
  //  hpct_timer_end("BayesMarkovChain");  TODO: revist timers
  return;
}

template <class P_V,class P_M>
void
uqStatisticalInverseProblemClass<P_V,P_M>::solveWithBayesMLSampling(
  const P_V& initialValues,
  const P_M* initialProposalCovMatrix)
{
  m_env.fullComm().Barrier();
  m_env.syncPrintDebugMsg("Entering uqStatisticalInverseProblemClass<P_V,P_M>::solveWithBayesMLSampling()",1,3000000,m_env.fullComm());

  if (m_options.m_computeSolution == false) {
    if ((m_env.subDisplayFile())) {
      *m_env.subDisplayFile() << "In uqStatisticalInverseProblemClass<P_V,P_M>::solveWithBayesMLSampling()"
                              << ": avoiding solution, as requested by user"
                              << std::endl;
    }
    return;
  }
  if ((m_env.subDisplayFile())) {
    *m_env.subDisplayFile() << "In uqStatisticalInverseProblemClass<P_V,P_M>::solveWithBayesMLSampling()"
                            << ": computing solution, as requested by user"
                            << std::endl;
  }

  if (m_mlSampler       ) delete m_mlSampler;
  if (m_mcSeqGenerator  ) delete m_mcSeqGenerator;
  if (m_solutionRealizer) delete m_solutionRealizer;
  if (m_subSolutionCdf  ) delete m_subSolutionCdf;
  if (m_subSolutionMdf  ) delete m_subSolutionMdf;
  if (m_solutionPdf     ) delete m_solutionPdf;
  if (m_solutionDomain  ) delete m_solutionDomain;

  P_V numEvaluationPointsVec(m_priorRv.imageSet().vectorSpace().zeroVector());
  numEvaluationPointsVec.cwSet(250.);

  // Compute output pdf up to a multiplicative constant: Bayesian approach
  m_solutionDomain = uqInstantiateIntersection(m_priorRv.pdf().domainSet(),m_likelihoodFunction.domainSet());

  m_solutionPdf = new uqBayesianJointPdfClass<P_V,P_M>(m_options.m_prefix.c_str(),
                                                       m_priorRv.pdf(),
                                                       m_likelihoodFunction,
                                                      *m_solutionDomain);

  m_postRv.setPdf(*m_solutionPdf);

  // Compute output realizer: ML approach
  m_chain = new uqSequenceOfVectorsClass<P_V,P_M>(m_postRv.imageSet().vectorSpace(),0,m_options.m_prefix+"chain");
  m_mlSampler = new uqMLSamplingClass<P_V,P_M>(m_options.m_prefix.c_str(),
                                               m_postRv,
                                               initialValues,
                                               initialProposalCovMatrix);

  m_mcSeqGenerator->generateSequence(*m_chain);

  m_solutionRealizer = new uqSequentialVectorRealizerClass<P_V,P_M>(m_options.m_prefix.c_str(),
                                                                   *m_chain);

  m_postRv.setRealizer(*m_solutionRealizer);

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << std::endl;
  }

  m_env.syncPrintDebugMsg("Leaving uqStatisticalInverseProblemClass<P_V,P_M>::solveWithBayesMLSampling()",1,3000000,m_env.fullComm());
  m_env.fullComm().Barrier();

  return;
}

template <class P_V,class P_M>
const uqBaseVectorRVClass<P_V,P_M>& 
uqStatisticalInverseProblemClass<P_V,P_M>::priorRv() const
{
  return m_priorRv;
}

template <class P_V,class P_M>
const uqGenericVectorRVClass<P_V,P_M>& 
uqStatisticalInverseProblemClass<P_V,P_M>::postRv() const
{
  return m_postRv;
}

template <class P_V,class P_M>
void
uqStatisticalInverseProblemClass<P_V,P_M>::print(std::ostream& os) const
{
  return;
}

template<class P_V,class P_M>
std::ostream& operator<<(std::ostream& os, const uqStatisticalInverseProblemClass<P_V,P_M>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_SIP_H__
