//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
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

#include <queso/StatisticalInverseProblem.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/GPMSA.h>
#include <queso/GslOptimizer.h>
#include <queso/OptimizerMonitor.h>

namespace QUESO {

// Default constructor -----------------------------
template <class P_V,class P_M>
StatisticalInverseProblem<P_V,P_M>::StatisticalInverseProblem(
  /*! The prefix                 */ const char*                               prefix,
  /*! Options (if no input file) */ const SipOptionsValues*            alternativeOptionsValues, // dakota
  /*! The prior RV               */ const BaseVectorRV      <P_V,P_M>& priorRv,
  /*! The likelihood function    */ const BaseScalarFunction<P_V,P_M>& likelihoodFunction,
  /*! The posterior RV           */       GenericVectorRV   <P_V,P_M>& postRv)
  :
  m_env                     (priorRv.env()),
  m_priorRv                 (priorRv),
  m_likelihoodFunction      (likelihoodFunction),
  m_postRv                  (postRv),
  m_solutionDomain          (NULL),
  m_solutionPdf             (NULL),
  m_subSolutionMdf          (NULL),
  m_subSolutionCdf          (NULL),
  m_solutionRealizer        (NULL),
  m_mhSeqGenerator          (NULL),
  m_mlSampler               (NULL),
  m_chain                   (NULL),
  m_logLikelihoodValues     (NULL),
  m_logTargetValues         (NULL),
  m_optionsObj              (alternativeOptionsValues),
  m_seedWithMAPEstimator    (false),
  m_userDidNotProvideOptions(false)
{
#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "Entering Sip" << std::endl;
#endif
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering StatisticalInverseProblem<P_V,P_M>::constructor()"
                            << ": prefix = " << prefix
                            << ", alternativeOptionsValues = " << alternativeOptionsValues
                            << ", m_env.optionsInputFileName() = " << m_env.optionsInputFileName()
                            << std::endl;
  }

  // If NULL, we create one
  if (m_optionsObj == NULL) {
    SipOptionsValues * tempOptions = new SipOptionsValues(&m_env, prefix);

    // We did this dance because scanOptionsValues is not a const method, but
    // m_optionsObj is a pointer to const
    m_optionsObj = tempOptions;

    // We set this flag so we don't delete the user-created object when it
    // comes time to deconstruct
    m_userDidNotProvideOptions = true;
  }

  if (m_optionsObj->m_help != "") {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << (*m_optionsObj) << std::endl;
    }
  }

#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "In Sip, finished scanning options" << std::endl;
#endif

  queso_require_equal_to_msg(priorRv.imageSet().vectorSpace().dimLocal(), likelihoodFunction.domainSet().vectorSpace().dimLocal(), "'priorRv' and 'likelihoodFunction' are related to vector spaces of different dimensions");

  queso_require_equal_to_msg(priorRv.imageSet().vectorSpace().dimLocal(), postRv.imageSet().vectorSpace().dimLocal(), "'priorRv' and 'postRv' are related to vector spaces of different dimensions");

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving StatisticalInverseProblem<P_V,P_M>::constructor()"
                            << ": prefix = " << m_optionsObj->m_prefix
                            << std::endl;
  }

  return;
}

template <class P_V,class P_M>
StatisticalInverseProblem<P_V,P_M>::StatisticalInverseProblem(
    const char * prefix,
    const SipOptionsValues * alternativeOptionsValues,
    const GPMSAFactory<P_V, P_M> & gpmsaFactory,
    GenericVectorRV <P_V,P_M> & postRv)
  :
  m_env                     (gpmsaFactory.m_totalPrior->env()),
  m_priorRv                 (*(gpmsaFactory.m_totalPrior)),
  m_likelihoodFunction      (gpmsaFactory.getGPMSAEmulator()),
  m_postRv                  (postRv),
  m_solutionDomain          (NULL),
  m_solutionPdf             (NULL),
  m_subSolutionMdf          (NULL),
  m_subSolutionCdf          (NULL),
  m_solutionRealizer        (NULL),
  m_mhSeqGenerator          (NULL),
  m_mlSampler               (NULL),
  m_chain                   (NULL),
  m_logLikelihoodValues     (NULL),
  m_logTargetValues         (NULL),
  m_optionsObj              (alternativeOptionsValues),
  m_seedWithMAPEstimator    (false),
  m_userDidNotProvideOptions(false)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering StatisticalInverseProblem<P_V,P_M>::constructor()"
                            << ": prefix = " << prefix
                            << ", alternativeOptionsValues = " << alternativeOptionsValues
                            << ", m_env.optionsInputFileName() = " << m_env.optionsInputFileName()
                            << std::endl;
  }

  // If NULL, we create one
  if (m_optionsObj == NULL) {
    SipOptionsValues * tempOptions = new SipOptionsValues(&m_env, prefix);

    // We did this dance because scanOptionsValues is not a const method, but
    // m_optionsObj is a pointer to const
    m_optionsObj = tempOptions;

    m_userDidNotProvideOptions = true;
  }

  if (m_optionsObj->m_help != "") {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << (*m_optionsObj) << std::endl;
    }
  }

  queso_require_equal_to_msg(m_priorRv.imageSet().vectorSpace().dimLocal(), m_likelihoodFunction.domainSet().vectorSpace().dimLocal(), "'priorRv' and 'likelihoodFunction' are related to vector spaces of different dimensions");

  queso_require_equal_to_msg(m_priorRv.imageSet().vectorSpace().dimLocal(), postRv.imageSet().vectorSpace().dimLocal(), "'priorRv' and 'postRv' are related to vector spaces of different dimensions");

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving StatisticalInverseProblem<P_V,P_M>::constructor()"
                            << ": prefix = " << m_optionsObj->m_prefix
                            << std::endl;
  }

  return;
}

// Destructor --------------------------------------
template <class P_V,class P_M>
StatisticalInverseProblem<P_V,P_M>::~StatisticalInverseProblem()
{
  if (m_chain) {
    m_chain->clear();
    delete m_chain;
  }
  if (m_logLikelihoodValues) {
    m_logLikelihoodValues->clear();
    delete m_logLikelihoodValues;
  }
  if (m_logTargetValues) {
    m_logTargetValues->clear();
    delete m_logTargetValues;
  }
  if (m_mlSampler       ) delete m_mlSampler;
  if (m_mhSeqGenerator  ) delete m_mhSeqGenerator;
  if (m_solutionRealizer) delete m_solutionRealizer;
  if (m_subSolutionCdf  ) delete m_subSolutionCdf;
  if (m_subSolutionMdf  ) delete m_subSolutionMdf;
  if (m_solutionPdf     ) delete m_solutionPdf;
  if (m_solutionDomain  ) delete m_solutionDomain;

  if (m_optionsObj && m_userDidNotProvideOptions) {
    delete m_optionsObj;
  }
}
// Statistical methods -----------------------------
template <class P_V,class P_M>
void
StatisticalInverseProblem<P_V,P_M>::solveWithBayesMetropolisHastings(
  const MhOptionsValues* alternativeOptionsValues, // dakota
  const P_V&                    initialValues,
  const P_M*                    initialProposalCovMatrix)
{
  //grvy_timer_begin("BayesMetropolisHastings"); TODO: revisit timing output
  //std::cout << "proc " << m_env.fullRank() << ", HERE sip 000" << std::endl;
  m_env.fullComm().Barrier();
  //std::cout << "proc " << m_env.fullRank() << ", HERE sip 001" << std::endl;
  m_env.fullComm().syncPrintDebugMsg("Entering StatisticalInverseProblem<P_V,P_M>::solveWithBayesMetropolisHastings()",1,3000000);

  if (m_optionsObj->m_computeSolution == false) {
    if ((m_env.subDisplayFile())) {
      *m_env.subDisplayFile() << "In StatisticalInverseProblem<P_V,P_M>::solveWithBayesMetropolisHastings()"
                              << ": avoiding solution, as requested by user"
                              << std::endl;
    }
    return;
  }
  if ((m_env.subDisplayFile())) {
    *m_env.subDisplayFile() << "In StatisticalInverseProblem<P_V,P_M>::solveWithBayesMetropolisHastings()"
                            << ": computing solution, as requested by user"
                            << std::endl;
  }

  queso_require_equal_to_msg(m_priorRv.imageSet().vectorSpace().dimLocal(), initialValues.sizeLocal(), "'m_priorRv' and 'initialValues' should have equal dimensions");

  if (initialProposalCovMatrix) {
    queso_require_equal_to_msg(m_priorRv.imageSet().vectorSpace().dimLocal(), initialProposalCovMatrix->numRowsLocal(), "'m_priorRv' and 'initialProposalCovMatrix' should have equal dimensions");
    queso_require_equal_to_msg(initialProposalCovMatrix->numCols(), initialProposalCovMatrix->numRowsGlobal(), "'initialProposalCovMatrix' should be a square matrix");
  }

  if (m_mlSampler       ) delete m_mlSampler;
  if (m_mhSeqGenerator  ) delete m_mhSeqGenerator;
  if (m_solutionRealizer) delete m_solutionRealizer;
  if (m_subSolutionCdf  ) delete m_subSolutionCdf;
  if (m_subSolutionMdf  ) delete m_subSolutionMdf;
  if (m_solutionPdf     ) delete m_solutionPdf;
  if (m_solutionDomain  ) delete m_solutionDomain;

  P_V numEvaluationPointsVec(m_priorRv.imageSet().vectorSpace().zeroVector());
  numEvaluationPointsVec.cwSet(250.);

  // Compute output pdf up to a multiplicative constant: Bayesian approach
  m_solutionDomain = InstantiateIntersection(m_priorRv.pdf().domainSet(),m_likelihoodFunction.domainSet());

  m_solutionPdf = new BayesianJointPdf<P_V,P_M>(m_optionsObj->m_prefix.c_str(),
                                                       m_priorRv.pdf(),
                                                       m_likelihoodFunction,
                                                       1.,
                                                       *m_solutionDomain);

  m_postRv.setPdf(*m_solutionPdf);
  m_chain = new SequenceOfVectors<P_V,P_M>(m_postRv.imageSet().vectorSpace(),0,m_optionsObj->m_prefix+"chain");

  // Decide whether or not to create a MetropolisHastingsSG instance from the
  // user-provided initial seed, or use the user-provided seed for a
  // deterministic optimisation instead and seed the chain with the result of
  // the optimisation
  if (this->m_seedWithMAPEstimator ||
      m_optionsObj->m_seedWithMAPEstimator) {
    // Unfortunately, I didn't have the foresight to make this an input file
    // option from the beginning, hence needing to check two flags

    // Do optimisation before sampling
    GslOptimizer optimizer(*m_solutionPdf);
    optimizer.setInitialPoint(dynamic_cast<const GslVector &>(initialValues));

    OptimizerMonitor monitor(m_env);
    monitor.set_display_output(true, true);

    // If the input file option is set, then use the monitor, otherwise don't
    if (m_optionsObj->m_useOptimizerMonitor) {
      optimizer.minimize(&monitor);
    }
    else {
      optimizer.minimize();
    }

    // Compute output realizer: Metropolis-Hastings approach
    m_mhSeqGenerator = new MetropolisHastingsSG<P_V, P_M>(
        m_optionsObj->m_prefix.c_str(), alternativeOptionsValues,
        m_postRv, optimizer.minimizer(), initialProposalCovMatrix);
  }
  else {
    // Compute output realizer: Metropolis-Hastings approach
    m_mhSeqGenerator = new MetropolisHastingsSG<P_V, P_M>(
        m_optionsObj->m_prefix.c_str(), alternativeOptionsValues, m_postRv,
        initialValues, initialProposalCovMatrix);
  }


  m_logLikelihoodValues = new ScalarSequence<double>(m_env, 0,
                                                     m_optionsObj->m_prefix +
                                                     "logLike");

  m_logTargetValues = new ScalarSequence<double>(m_env, 0,
                                                 m_optionsObj->m_prefix +
                                                 "logTarget");

  // m_logLikelihoodValues and m_logTargetValues may be NULL
  m_mhSeqGenerator->generateSequence(*m_chain, m_logLikelihoodValues,
                                     m_logTargetValues);

  m_solutionRealizer = new SequentialVectorRealizer<P_V,P_M>(m_optionsObj->m_prefix.c_str(),
                                                                    *m_chain);

  m_postRv.setRealizer(*m_solutionRealizer);

  m_env.fullComm().syncPrintDebugMsg("In StatisticalInverseProblem<P_V,P_M>::solveWithBayesMetropolisHastings(), code place 1",3,3000000);
  //m_env.fullComm().Barrier();

  // Compute output mdf: uniform sampling approach
#ifdef UQ_ALSO_COMPUTE_MDFS_WITHOUT_KDE
  m_subMdfGrids  = new ArrayOfOneDGrids <P_V,P_M>((m_optionsObj->m_prefix+"Mdf_").c_str(),m_postRv.imageSet().vectorSpace());
  m_subMdfValues = new ArrayOfOneDTables<P_V,P_M>((m_optionsObj->m_prefix+"Mdf_").c_str(),m_postRv.imageSet().vectorSpace());
  m_chain->subUniformlySampledMdf(numEvaluationPointsVec, // input
                                  *m_subMdfGrids,         // output
                                  *m_subMdfValues);       // output
  m_subSolutionMdf = new SampledVectorMdf<P_V,P_M>(m_optionsObj->m_prefix.c_str(),
                                                          *m_subMdfGrids,
                                                          *m_subMdfValues);
  m_postRv.setMdf(*m_subSolutionMdf);

  if ((m_optionsObj->m_dataOutputFileName                       != UQ_SIP_FILENAME_FOR_NO_FILE                    ) &&
      (m_optionsObj->m_dataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_dataOutputAllowedSet.end())) {
    if (m_env.subRank() == 0) {
      // Write data output file
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "Opening data output file '" << m_optionsObj->m_dataOutputFileName
                                << "' for calibration problem with problem with prefix = " << m_optionsObj->m_prefix
                                << std::endl;
      }

      // Open file
      // Always write at the end of an eventual pre-existing file
      std::ofstream* ofsvar = new std::ofstream((m_optionsObj->m_dataOutputFileName+"_sub"+m_env.subIdString()+".m").c_str(), std::ofstream::out | std::ofstream::in | std::ofstream::ate);
      if ((ofsvar            == NULL ) ||
          (ofsvar->is_open() == false)) {
        delete ofsvar;
        ofsvar = new std::ofstream((m_optionsObj->m_dataOutputFileName+"_sub"+m_env.subIdString()+".m").c_str(), std::ofstream::out | std::ofstream::trunc);
      }
      queso_require_msg((ofsvar && ofsvar->is_open()), "failed to open file");

      m_postRv.mdf().print(*ofsvar);

      // Close file
      //ofsvar->close();
      delete ofsvar;
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "Closed data output file '" << m_optionsObj->m_dataOutputFileName
                                << "' for calibration problem with problem with prefix = " << m_optionsObj->m_prefix
                                << std::endl;
      }
    }
  }
#endif
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << std::endl;
  }

  m_env.fullComm().syncPrintDebugMsg("Leaving StatisticalInverseProblem<P_V,P_M>::solveWithBayesMetropolisHastings()",1,3000000);
  m_env.fullComm().Barrier();
  // grvy_timer_end("BayesMetropolisHastings"); TODO: revisit timers
  return;
}

template <class P_V, class P_M>
void
StatisticalInverseProblem<P_V, P_M>::seedWithMAPEstimator()
{
  this->m_seedWithMAPEstimator = true;
}

template <class P_V,class P_M>
void
StatisticalInverseProblem<P_V,P_M>::solveWithBayesMLSampling()
{
  m_env.fullComm().Barrier();
  m_env.fullComm().syncPrintDebugMsg("Entering StatisticalInverseProblem<P_V,P_M>::solveWithBayesMLSampling()",1,3000000);

  if (m_optionsObj->m_computeSolution == false) {
    if ((m_env.subDisplayFile())) {
      *m_env.subDisplayFile() << "In StatisticalInverseProblem<P_V,P_M>::solveWithBayesMLSampling()"
                              << ": avoiding solution, as requested by user"
                              << std::endl;
    }
    return;
  }
  if ((m_env.subDisplayFile())) {
    *m_env.subDisplayFile() << "In StatisticalInverseProblem<P_V,P_M>::solveWithBayesMLSampling()"
                            << ": computing solution, as requested by user"
                            << std::endl;
  }

  if (m_mlSampler       ) delete m_mlSampler;
  if (m_mhSeqGenerator  ) delete m_mhSeqGenerator;
  if (m_solutionRealizer) delete m_solutionRealizer;
  if (m_subSolutionCdf  ) delete m_subSolutionCdf;
  if (m_subSolutionMdf  ) delete m_subSolutionMdf;
  if (m_solutionPdf     ) delete m_solutionPdf;
  if (m_solutionDomain  ) delete m_solutionDomain;

  P_V numEvaluationPointsVec(m_priorRv.imageSet().vectorSpace().zeroVector());
  numEvaluationPointsVec.cwSet(250.);

  // Compute output pdf up to a multiplicative constant: Bayesian approach
  m_solutionDomain = InstantiateIntersection(m_priorRv.pdf().domainSet(),m_likelihoodFunction.domainSet());

  m_solutionPdf = new BayesianJointPdf<P_V,P_M>(m_optionsObj->m_prefix.c_str(),
                                                       m_priorRv.pdf(),
                                                       m_likelihoodFunction,
                                                       1.,
                                                       *m_solutionDomain);

  m_postRv.setPdf(*m_solutionPdf);

  // Compute output realizer: ML approach
  m_chain = new SequenceOfVectors<P_V,P_M>(m_postRv.imageSet().vectorSpace(),0,m_optionsObj->m_prefix+"chain");
  m_mlSampler = new MLSampling<P_V,P_M>(m_optionsObj->m_prefix.c_str(),
                                             //m_postRv,
                                               m_priorRv,
                                               m_likelihoodFunction);
  //                                           initialValues,
  //                                           initialProposalCovMatrix);

  m_mlSampler->generateSequence(*m_chain,
                                NULL,
                                NULL);

  m_solutionRealizer = new SequentialVectorRealizer<P_V,P_M>(m_optionsObj->m_prefix.c_str(),
                                                                    *m_chain);

  m_postRv.setRealizer(*m_solutionRealizer);

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << std::endl;
  }

  m_env.fullComm().syncPrintDebugMsg("Leaving StatisticalInverseProblem<P_V,P_M>::solveWithBayesMLSampling()",1,3000000);
  m_env.fullComm().Barrier();

  return;
}

template <class P_V, class P_M>
const MetropolisHastingsSG<P_V, P_M> &
StatisticalInverseProblem<P_V, P_M>::sequenceGenerator() const
{
  return *m_mhSeqGenerator;
}

//--------------------------------------------------
template <class P_V,class P_M>
const BaseVectorRV<P_V,P_M>&
StatisticalInverseProblem<P_V,P_M>::priorRv() const
{
  return m_priorRv;
}
//--------------------------------------------------
template <class P_V,class P_M>
const GenericVectorRV<P_V,P_M>&
StatisticalInverseProblem<P_V,P_M>::postRv() const
{
  return m_postRv;
}
//--------------------------------------------------
template <class P_V,class P_M>
const BaseVectorSequence<P_V,P_M>&
StatisticalInverseProblem<P_V,P_M>::chain() const
{
  queso_require_msg(m_chain, "m_chain is NULL");
  return *m_chain;
}
//--------------------------------------------------
template <class P_V,class P_M>
const ScalarSequence<double>&
StatisticalInverseProblem<P_V,P_M>::logLikelihoodValues() const
{
  queso_require_msg(m_logLikelihoodValues, "m_logLikelihoodValues is NULL");
  return *m_logLikelihoodValues;
}
//--------------------------------------------------
template <class P_V,class P_M>
const ScalarSequence<double>&
StatisticalInverseProblem<P_V,P_M>::logTargetValues() const
{
  queso_require_msg(m_logTargetValues, "m_logTargetValues is NULL");
  return *m_logTargetValues;
}
//--------------------------------------------------
template <class P_V,class P_M>
double StatisticalInverseProblem<P_V,P_M>::logEvidence() const
{
  queso_require_msg(m_mlSampler, "m_mlSampler is NULL");
  return m_mlSampler->logEvidence();
}
//--------------------------------------------------
template <class P_V,class P_M>
double StatisticalInverseProblem<P_V,P_M>::meanLogLikelihood() const
{
  queso_require_msg(m_mlSampler, "m_mlSampler is NULL");
  return m_mlSampler->meanLogLikelihood();
}
//--------------------------------------------------
template <class P_V,class P_M>
double StatisticalInverseProblem<P_V,P_M>::eig() const
{
  queso_require_msg(m_mlSampler, "m_mlSampler is NULL");
  return m_mlSampler->eig();
}
// I/O methods--------------------------------------
template <class P_V,class P_M>
void
StatisticalInverseProblem<P_V,P_M>::print(std::ostream& os) const
{
  return;
}

}  // End namespace QUESO

template class QUESO::StatisticalInverseProblem<QUESO::GslVector, QUESO::GslMatrix>;
