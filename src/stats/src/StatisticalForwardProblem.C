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

#include <queso/StatisticalForwardProblem.h>
#include <queso/SequentialVectorRealizer.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Default constructor -----------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>::StatisticalForwardProblem(
  /*! The prefix                 */ const char*                                       prefix,
  /*! Options (if no input file) */ const SfpOptionsValues*                    alternativeOptionsValues, // dakota
  /*! The input RV               */ const BaseVectorRV      <P_V,P_M>&         paramRv,
  /*! The QoI function           */ const BaseVectorFunction<P_V,P_M,Q_V,Q_M>& qoiFunction,
  /*! The QoI RV                 */       GenericVectorRV   <Q_V,Q_M>&         qoiRv)
  :
  m_env                     (paramRv.env()),
  m_paramRv                 (paramRv),
  m_qoiFunction             (qoiFunction),
  m_qoiRv                   (qoiRv),
  m_paramChain              (NULL),
  m_qoiChain                (NULL),
  m_mcSeqGenerator          (NULL),
  m_solutionRealizer        (NULL),
#ifdef UQ_ALSO_COMPUTE_MDFS_WITHOUT_KDE
  m_subMdfGrids             (NULL),
  m_subMdfValues            (NULL),
#endif
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
  m_subSolutionMdf          (NULL),
  m_subCdfGrids             (NULL),
  m_subCdfValues            (NULL),
  m_subSolutionCdf          (NULL),
  m_unifiedCdfGrids         (NULL),
  m_unifiedCdfValues        (NULL),
  m_unifiedSolutionCdf      (NULL),
#endif
  m_solutionPdf             (NULL),
  m_optionsObj              (alternativeOptionsValues),
  m_userDidNotProvideOptions(false)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = " << prefix
                            << ", alternativeOptionsValues = " << alternativeOptionsValues
                            << ", m_env.optionsInputFileName() = " << m_env.optionsInputFileName()
                            << std::endl;
  }

  // If NULL, we create one
  if (m_optionsObj == NULL) {
    SfpOptionsValues * tempOptions = new SfpOptionsValues(&m_env, prefix);

    // We did this dance because scanOptionsValues is not a const method, but
    // m_optionsObj is a pointer to const
    m_optionsObj = tempOptions;

    // We do this so we don't delete the user's object in the dtor
    m_userDidNotProvideOptions = true;
  }

  if (m_optionsObj->m_help != "") {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << (*m_optionsObj) << std::endl;
    }
  }

  queso_require_equal_to_msg(paramRv.imageSet().vectorSpace().dimLocal(), qoiFunction.domainSet().vectorSpace().dimLocal(), "'paramRv' and 'qoiFunction' are related to vector spaces of different dimensions");

  queso_require_equal_to_msg(qoiFunction.imageSet().vectorSpace().dimLocal(), qoiRv.imageSet().vectorSpace().dimLocal(), "'qoiFunction' and 'qoiRv' are related to vector spaces of different dimensions");

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = " << m_optionsObj->m_prefix
                            << std::endl;
  }
}

// Destructor --------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>::~StatisticalForwardProblem()
{
  if (m_solutionPdf       ) delete m_solutionPdf;

#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
  if (m_unifiedSolutionCdf) delete m_unifiedSolutionCdf;
  if (m_unifiedCdfValues  ) delete m_unifiedCdfValues;
  if (m_unifiedCdfGrids   ) delete m_unifiedCdfGrids;

  if (m_subSolutionCdf    ) delete m_subSolutionCdf;
  if (m_subCdfValues      ) delete m_subCdfValues;
  if (m_subCdfGrids       ) delete m_subCdfGrids;

  if (m_subSolutionMdf    ) delete m_subSolutionMdf;
#endif
#ifdef UQ_ALSO_COMPUTE_MDFS_WITHOUT_KDE
  if (m_subMdfValues      ) delete m_subMdfValues;
  if (m_subMdfGrids       ) delete m_subMdfGrids;
#endif
  if (m_solutionRealizer  ) delete m_solutionRealizer;

  if (m_mcSeqGenerator    ) delete m_mcSeqGenerator;

  if (m_qoiChain) {
    m_qoiChain->clear();
    delete m_qoiChain;
  }

  if (m_paramChain) {
    m_paramChain->clear();
    delete m_paramChain;
  }

  if (m_optionsObj        ) delete m_optionsObj;
}

// Statistical methods -----------------------------
template <class P_V,class P_M, class Q_V, class Q_M>
bool
  StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>::computeSolutionFlag() const
{
  return m_optionsObj->m_computeSolution;
}
//--------------------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
void
StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>::solveWithMonteCarlo(
  const McOptionsValues* alternativeOptionsValues)
{
  m_env.fullComm().Barrier();
  m_env.fullComm().syncPrintDebugMsg("Entering StatisticalForwardProblem<P_V,P_M>::solveWithMonteCarlo()",1,3000000);

  if (m_optionsObj->m_computeSolution == false) {
    if ((m_env.subDisplayFile())) {
      *m_env.subDisplayFile() << "In StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>::solveWithMonteCarlo()"
                              << ": avoiding solution, as requested by user"
                              << std::endl;
    }
    return;
  }
  if ((m_env.subDisplayFile())) {
    *m_env.subDisplayFile() << "In StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>::solveWithMonteCarlo()"
                            << ": computing solution, as requested by user"
                            << std::endl;
  }

  if (m_solutionPdf       ) delete m_solutionPdf;

#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
  if (m_unifiedSolutionCdf) delete m_unifiedSolutionCdf;
  if (m_unifiedCdfValues  ) delete m_unifiedCdfValues;
  if (m_unifiedCdfGrids   ) delete m_unifiedCdfGrids;

  if (m_subSolutionCdf    ) delete m_subSolutionCdf;
  if (m_subCdfValues      ) delete m_subCdfValues;
  if (m_subCdfGrids       ) delete m_subCdfGrids;

  if (m_subSolutionMdf    ) delete m_subSolutionMdf;
#endif
#ifdef UQ_ALSO_COMPUTE_MDFS_WITHOUT_KDE
  if (m_subMdfValues      ) delete m_subMdfValues;
  if (m_subMdfGrids       ) delete m_subMdfGrids;
#endif
  if (m_solutionRealizer  ) delete m_solutionRealizer;

  if (m_mcSeqGenerator    ) delete m_mcSeqGenerator;

  if (m_qoiChain) {
    m_qoiChain->clear();
    delete m_qoiChain;
  }

  if (m_paramChain) {
    m_paramChain->clear();
    delete m_paramChain;
  }

  Q_V numEvaluationPointsVec(m_qoiRv.imageSet().vectorSpace().zeroVector());
  numEvaluationPointsVec.cwSet(250.);

  // Compute output realizer: Monte Carlo approach
  m_paramChain = new SequenceOfVectors<P_V,P_M>(m_paramRv.imageSet().vectorSpace(),0,m_optionsObj->m_prefix+"paramChain");
  m_qoiChain   = new SequenceOfVectors<Q_V,Q_M>(m_qoiRv.imageSet().vectorSpace(),  0,m_optionsObj->m_prefix+"qoiChain"  );
  m_mcSeqGenerator = new MonteCarloSG<P_V,P_M,Q_V,Q_M>(m_optionsObj->m_prefix.c_str(),
                                                              alternativeOptionsValues,
                                                              m_paramRv,
                                                              m_qoiFunction);
  //m_qoiRv);
  m_mcSeqGenerator->generateSequence(*m_paramChain,
                                     *m_qoiChain);
  m_solutionRealizer = new SequentialVectorRealizer<Q_V,Q_M>((m_optionsObj->m_prefix+"Qoi").c_str(),
                                                                    *m_qoiChain);
  m_qoiRv.setRealizer(*m_solutionRealizer);

  // Compute output mdf: uniform sampling approach
#ifdef UQ_ALSO_COMPUTE_MDFS_WITHOUT_KDE
  m_subMdfGrids  = new ArrayOfOneDGrids <Q_V,Q_M>((m_optionsObj->m_prefix+"QoiMdf_").c_str(),m_qoiRv.imageSet().vectorSpace());
  m_subMdfValues = new ArrayOfOneDTables<Q_V,Q_M>((m_optionsObj->m_prefix+"QoiMdf_").c_str(),m_qoiRv.imageSet().vectorSpace());
  m_qoiChain->subUniformlySampledMdf(numEvaluationPointsVec, // input
                                     *m_subMdfGrids,         // output
                                     *m_subMdfValues);       // output

  m_subSolutionMdf = new SampledVectorMdf<Q_V,Q_M>((m_optionsObj->m_prefix+"Qoi").c_str(),
                                                          *m_subMdfGrids,
                                                          *m_subMdfValues);
  m_qoiRv.setMdf(*m_subSolutionMdf);
#endif

  // Compute output cdf: uniform sampling approach
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
  std::string subCoreName_qoiCdf(m_optionsObj->m_prefix+    "QoiCdf_");
  std::string uniCoreName_qoiCdf(m_optionsObj->m_prefix+"unifQoiCdf_");
  if (m_env.numSubEnvironments() == 1) subCoreName_qoiCdf = uniCoreName_qoiCdf;

  std::string subCoreName_solutionCdf(m_optionsObj->m_prefix+    "Qoi");
  std::string uniCoreName_solutionCdf(m_optionsObj->m_prefix+"unifQoi");
  if (m_env.numSubEnvironments() == 1) subCoreName_solutionCdf = uniCoreName_solutionCdf;

  m_subCdfGrids  = new ArrayOfOneDGrids <Q_V,Q_M>(subCoreName_qoiCdf.c_str(),m_qoiRv.imageSet().vectorSpace());
  m_subCdfValues = new ArrayOfOneDTables<Q_V,Q_M>(subCoreName_qoiCdf.c_str(),m_qoiRv.imageSet().vectorSpace());
  m_qoiChain->subUniformlySampledCdf(numEvaluationPointsVec, // input
                                     *m_subCdfGrids,         // output
                                     *m_subCdfValues);       // output

  m_subSolutionCdf = new SampledVectorCdf<Q_V,Q_M>(subCoreName_solutionCdf.c_str(),
                                                          *m_subCdfGrids,
                                                          *m_subCdfValues);
  m_qoiRv.setSubCdf(*m_subSolutionCdf);

  // Compute unified cdf if necessary
  if (m_env.numSubEnvironments() == 1) {
    m_qoiRv.setUnifiedCdf(*m_subSolutionCdf);
  }
  else {
    m_unifiedCdfGrids  = new ArrayOfOneDGrids <Q_V,Q_M>(uniCoreName_qoiCdf.c_str(),m_qoiRv.imageSet().vectorSpace());
    m_unifiedCdfValues = new ArrayOfOneDTables<Q_V,Q_M>(uniCoreName_qoiCdf.c_str(),m_qoiRv.imageSet().vectorSpace());
    m_qoiChain->unifiedUniformlySampledCdf(numEvaluationPointsVec, // input
                                           *m_unifiedCdfGrids,     // output
                                           *m_unifiedCdfValues);   // output

    m_unifiedSolutionCdf = new SampledVectorCdf<Q_V,Q_M>(uniCoreName_solutionCdf.c_str(),
                                                                *m_unifiedCdfGrids,
                                                                *m_unifiedCdfValues);
    m_qoiRv.setUnifiedCdf(*m_unifiedSolutionCdf);
  }
#endif
  // Compute (just unified one) covariance matrix, if requested
  // Compute (just unified one) correlation matrix, if requested
  P_M* pqCovarianceMatrix  = NULL;
  P_M* pqCorrelationMatrix = NULL;
  if (m_optionsObj->m_computeCovariances || m_optionsObj->m_computeCorrelations) {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>::solveWithMonteCarlo()"
                              << ", prefix = " << m_optionsObj->m_prefix
                              << ": instantiating cov and corr matrices"
                              << std::endl;
    }

    // Only compute correlations on the inter0Comm communicator
    if (m_env.subRank() == 0) {
      pqCovarianceMatrix = new P_M(m_env,
                                   m_paramRv.imageSet().vectorSpace().map(),       // number of rows
                                   m_qoiRv.imageSet().vectorSpace().dimGlobal());  // number of cols
      pqCorrelationMatrix = new P_M(m_env,
                                    m_paramRv.imageSet().vectorSpace().map(),      // number of rows
                                    m_qoiRv.imageSet().vectorSpace().dimGlobal()); // number of cols
      ComputeCovCorrMatricesBetweenVectorSequences(*m_paramChain,
                                                     *m_qoiChain,
                                                     std::min(m_paramRv.realizer().subPeriod(),m_qoiRv.realizer().subPeriod()), // FIX ME: might be INFINITY
                                                     *pqCovarianceMatrix,
                                                     *pqCorrelationMatrix);
    }
  }

  // Write data out
  if (m_env.subDisplayFile()) {
    if (pqCovarianceMatrix ) {
      *m_env.subDisplayFile() << "In StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>::solveWithMonteCarlo()"
                              << ", prefix = "                           << m_optionsObj->m_prefix
                              << ": contents of covariance matrix are\n" << *pqCovarianceMatrix
                              << std::endl;
    }
    if (pqCorrelationMatrix) {
      *m_env.subDisplayFile() << "In StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>::solveWithMonteCarlo()"
                              << ", prefix = "                            << m_optionsObj->m_prefix
                              << ": contents of correlation matrix are\n" << *pqCorrelationMatrix
                              << std::endl;
    }
  }

  // Open data output file
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "In StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>::solveWithMonteCarlo()"
                            << ", prefix = "                                        << m_optionsObj->m_prefix
                            << ": checking necessity of opening data output file '" << m_optionsObj->m_dataOutputFileName
                            << "'"
                            << std::endl;
  }
  FilePtrSetStruct filePtrSet;
  if (m_env.openOutputFile(m_optionsObj->m_dataOutputFileName,
                           UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT, // Yes, always ".m"
                           m_optionsObj->m_dataOutputAllowedSet,
                           false,
                           filePtrSet)) {
#ifdef UQ_ALSO_COMPUTE_MDFS_WITHOUT_KDE
    m_qoiRv.mdf().print(*filePtrSet.ofsVar);
#endif
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
    *filePtrSet.ofsVar << m_qoiRv.subCdf();
#endif

    //if (pqCovarianceMatrix ) *filePtrSet.ofsVar << *pqCovarianceMatrix;  // FIX ME: output matrix in matlab format
    //if (pqCorrelationMatrix) *filePtrSet.ofsVar << *pqCorrelationMatrix; // FIX ME: output matrix in matlab format

    // Write unified cdf if necessary
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
    if (m_env.numSubEnvironments() > 1) {
      if (m_qoiRv.imageSet().vectorSpace().numOfProcsForStorage() == 1) {
        if (m_env.inter0Rank() == 0) {
          *filePtrSet.ofsVar << m_qoiRv.unifiedCdf(); //*m_unifiedSolutionCdf;
        }
      }
      else {
        queso_error_msg("unified cdf writing, parallel vectors not supported yet");
      }
    }
#endif
    // Close data output file
    m_env.closeFile(filePtrSet,UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>::solveWithMonteCarlo()"
                              << ", prefix = "                 << m_optionsObj->m_prefix
                              << ": closed data output file '" << m_optionsObj->m_dataOutputFileName
                              << "'"
                              << std::endl;
    }
  }
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << std::endl;
  }

  if (pqCovarianceMatrix ) delete pqCovarianceMatrix;
  if (pqCorrelationMatrix) delete pqCorrelationMatrix;

  m_env.fullComm().syncPrintDebugMsg("Leaving StatisticalForwardProblem<P_V,P_M>::solveWithMonteCarlo()",1,3000000);
  m_env.fullComm().Barrier();

  return;
}
//--------------------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
const GenericVectorRV<Q_V,Q_M>&
StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>::qoiRv() const
{
  return m_qoiRv;
}
//--------------------------------------------------

#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
template <class P_V,class P_M,class Q_V,class Q_M>
const BaseVectorCdf<Q_V,Q_M>&
StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>::qoiRv_unifiedCdf() const
{
  if (m_env.numSubEnvironments() == 1) {
    return m_qoiRv.subCdf();
  }

  if (m_env.inter0Rank() < 0) {
    return m_qoiRv.subCdf();
  }


  //                    m_env.worldRank(),
  //                    "StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>::qoiRv_unifiedCdf()",
  //                    "variable is NULL");
  return m_qoiRv.unifiedCdf(); //*m_unifiedSolutionCdf;
}
#endif
//--------------------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
  const BaseVectorSequence<Q_V,Q_M>&
  StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>::getParamChain() const
{

  // Make sure this runs after the forward propagation
  // only then we obtain the actual realizations of the parameters
  queso_require_msg(m_paramChain, "m_paramChain is NULL");

  return *m_paramChain;

}
// I/O methods--------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
void
StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>::print(std::ostream& os) const
{
  return;
}

}  // End namespace QUESO

template class QUESO::StatisticalForwardProblem<QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix>;
