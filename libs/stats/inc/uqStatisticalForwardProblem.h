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

#ifndef __UQ_PROPAG_PROBLEM_H__
#define __UQ_PROPAG_PROBLEM_H__

#include <uqVectorFunction.h>
#include <uqMonteCarloSG.h>
#include <uqVectorRV.h>
#include <uqSequenceOfVectors.h>

#undef UQ_PROPAG_PROBLEM_READS_SOLVER_OPTION

#define UQ_PROPAG_PROBLEM_FILENAME_FOR_NO_FILE "."

// _ODV = option default value
#define UQ_PROPAG_PROBLEM_COMPUTE_SOLUTION_ODV      1
#define UQ_PROPAG_PROBLEM_COMPUTE_COVARIANCES_ODV   1
#define UQ_PROPAG_PROBLEM_COMPUTE_CORRELATIONS_ODV  1
#define UQ_PROPAG_PROBLEM_DATA_OUTPUT_FILE_NAME_ODV UQ_PROPAG_PROBLEM_FILENAME_FOR_NO_FILE
#define UQ_PROPAG_PROBLEM_DATA_OUTPUT_ALLOW_ODV     ""
#ifdef UQ_PROPAG_PROBLEM_READS_SOLVER_OPTION
#define UQ_PROPAG_PROBLEM_SOLVER_ODV                "mc" // Monte Carlo
#endif

/*! A templated class that represents statistical forward problems.
 */
template <class P_V,class P_M,class Q_V,class Q_M>
class uqStatisticalForwardProblemClass
{
public:
  uqStatisticalForwardProblemClass(const char*                                       prefix,      /*! The prefix.                  */
                                   const uqBaseVectorRVClass      <P_V,P_M>&         paramRv,     /*! The input random variable.   */
                                   const uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>& qoiFunction, /*! The function to compute qoi. */
                                   uqGenericVectorRVClass         <Q_V,Q_M>&         qoiRv);      /*! The output random variable.  */
 ~uqStatisticalForwardProblemClass();

        bool                             computeSolutionFlag() const;
        void                             solveWithMonteCarlo();
  const uqGenericVectorRVClass<Q_V,Q_M>& qoiRv              () const;
  const uqBaseVectorCdfClass  <Q_V,Q_M>& qoiRv_unifiedCdf   () const;

        void                             print              (std::ostream& os) const;

private:
        void                             commonConstructor  ();
        void                             defineMyOptions    (po::options_description& optionsDesc);
        void                             getMyOptionValues  (po::options_description& optionsDesc);

  const uqBaseEnvironmentClass&                     m_env;
        std::string                                 m_prefix;

        po::options_description*                    m_optionsDesc;
        std::string                                 m_option_help;
	std::string                                 m_option_computeSolution;
	std::string                                 m_option_computeCovariances;
	std::string                                 m_option_computeCorrelations;
        std::string                                 m_option_dataOutputFileName;
        std::string                                 m_option_dataOutputAllowedSet;
#ifdef UQ_PROPAG_PROBLEM_READS_SOLVER_OPTION
        std::string                                 m_option_solver;
#endif

        bool                                        m_computeSolution;
        bool                                        m_computeCovariances;
        bool                                        m_computeCorrelations;
        std::string                                 m_dataOutputFileName;
        std::set<unsigned int>                      m_dataOutputAllowedSet;
#ifdef UQ_PROPAG_PROBLEM_READS_SOLVER_OPTION
	std::string                                 m_solverString;
#endif

  const uqBaseVectorRVClass      <P_V,P_M>&         m_paramRv;
  const uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>& m_qoiFunction;
        uqGenericVectorRVClass   <Q_V,Q_M>&         m_qoiRv; // FIX ME: Maybe not always generic ?

        uqBaseVectorSequenceClass<Q_V,Q_M>*         m_paramChain;
        uqBaseVectorSequenceClass<Q_V,Q_M>*         m_qoiChain;
        uqMonteCarloSGClass      <P_V,P_M,Q_V,Q_M>* m_mcSeqGenerator;

        uqBaseVectorRealizerClass<Q_V,Q_M>*         m_solutionRealizer;
 
        uqArrayOfOneDGridsClass  <Q_V,Q_M>*         m_subMdfGrids;
        uqArrayOfOneDTablesClass <Q_V,Q_M>*         m_subMdfValues;
        uqBaseVectorMdfClass     <Q_V,Q_M>*         m_subSolutionMdf;

        uqArrayOfOneDGridsClass  <Q_V,Q_M>*         m_subCdfGrids;
        uqArrayOfOneDTablesClass <Q_V,Q_M>*         m_subCdfValues;
        uqBaseVectorCdfClass     <Q_V,Q_M>*         m_subSolutionCdf;

        uqArrayOfOneDGridsClass  <Q_V,Q_M>*         m_unifiedCdfGrids;
        uqArrayOfOneDTablesClass <Q_V,Q_M>*         m_unifiedCdfValues;
        uqBaseVectorCdfClass     <Q_V,Q_M>*         m_unifiedSolutionCdf;

        uqBaseJointPdfClass      <Q_V,Q_M>*         m_solutionPdf;
};

template<class P_V,class P_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M>& obj);

template <class P_V,class P_M,class Q_V,class Q_M>
uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M>::uqStatisticalForwardProblemClass(
  const char*                                       prefix,
  const uqBaseVectorRVClass      <P_V,P_M>&         paramRv,
  const uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>& qoiFunction,
        uqGenericVectorRVClass   <Q_V,Q_M>&         qoiRv)
  :
  m_env                        (paramRv.env()),
  m_prefix                     ((std::string)(prefix) + "fp_"),
  m_optionsDesc                (new po::options_description("UQ Propagation Problem")),
  m_option_help                (m_prefix + "help"                ),
  m_option_computeSolution     (m_prefix + "computeSolution"     ),
  m_option_computeCovariances  (m_prefix + "computeCovariances"  ),
  m_option_computeCorrelations (m_prefix + "computeCorrelations" ),
  m_option_dataOutputFileName  (m_prefix + "dataOutputFileName"  ),
  m_option_dataOutputAllowedSet(m_prefix + "dataOutputAllowedSet"),
#ifdef UQ_PROPAG_PROBLEM_READS_SOLVER_OPTION
  m_option_solver              (m_prefix + "solver"              ),
#endif
  m_computeSolution       (UQ_PROPAG_PROBLEM_COMPUTE_SOLUTION_ODV     ),
  m_computeCovariances    (UQ_PROPAG_PROBLEM_COMPUTE_COVARIANCES_ODV  ),
  m_computeCorrelations   (UQ_PROPAG_PROBLEM_COMPUTE_CORRELATIONS_ODV ),
  m_dataOutputFileName    (UQ_PROPAG_PROBLEM_DATA_OUTPUT_FILE_NAME_ODV),
//m_dataOutputAllowedSet  (),
#ifdef UQ_PROPAG_PROBLEM_READS_SOLVER_OPTION
  m_solverString          (UQ_PROPAG_PROBLEM_SOLVER_ODV               ),
#endif
  m_paramRv               (paramRv),
  m_qoiFunction           (qoiFunction),
  m_qoiRv                 (qoiRv),
  m_paramChain            (NULL),
  m_qoiChain              (NULL),
  m_mcSeqGenerator        (NULL),
  m_solutionRealizer      (NULL),
  m_subMdfGrids           (NULL),
  m_subMdfValues          (NULL),
  m_subSolutionMdf        (NULL),
  m_subCdfGrids           (NULL),
  m_subCdfValues          (NULL),
  m_subSolutionCdf        (NULL),
  m_unifiedCdfGrids       (NULL),
  m_unifiedCdfValues      (NULL),
  m_unifiedSolutionCdf    (NULL),
  m_solutionPdf           (NULL)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M>::constructor()"
                           << ": prefix = "              << m_prefix
                           << std::endl;
  }

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "In uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M>::constructor()"
                           << ": after getting values of options, state of object is:"
                           << "\n" << *this
                           << std::endl;
  }

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }
}

template <class P_V,class P_M,class Q_V,class Q_M>
uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M>::~uqStatisticalForwardProblemClass()
{
  if (m_solutionPdf       ) delete m_solutionPdf;

  if (m_unifiedSolutionCdf) delete m_unifiedSolutionCdf;
  if (m_unifiedCdfValues  ) delete m_unifiedCdfValues;
  if (m_unifiedCdfGrids   ) delete m_unifiedCdfGrids;

  if (m_subSolutionCdf    ) delete m_subSolutionCdf;
  if (m_subCdfValues      ) delete m_subCdfValues;
  if (m_subCdfGrids       ) delete m_subCdfGrids;

  if (m_subSolutionMdf    ) delete m_subSolutionMdf;
  if (m_subMdfValues      ) delete m_subMdfValues;
  if (m_subMdfGrids       ) delete m_subMdfGrids;

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

  if (m_optionsDesc       ) delete m_optionsDesc;
}

template<class P_V,class P_M,class Q_V,class Q_M>
void
uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M>::defineMyOptions(
  po::options_description& optionsDesc)
{
  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                                                       "produce help message for propagation problem")
    (m_option_computeSolution.c_str(),      po::value<bool       >()->default_value(UQ_PROPAG_PROBLEM_COMPUTE_SOLUTION_ODV     ), "compute solution process"                    )
    (m_option_computeCovariances.c_str(),   po::value<bool       >()->default_value(UQ_PROPAG_PROBLEM_COMPUTE_COVARIANCES_ODV  ), "compute pq covariances"                      )
    (m_option_computeCorrelations.c_str(),  po::value<bool       >()->default_value(UQ_PROPAG_PROBLEM_COMPUTE_CORRELATIONS_ODV ), "compute pq correlations"                     )
    (m_option_dataOutputFileName.c_str(),   po::value<std::string>()->default_value(UQ_PROPAG_PROBLEM_DATA_OUTPUT_FILE_NAME_ODV), "name of data output file"                    )
    (m_option_dataOutputAllowedSet.c_str(), po::value<std::string>()->default_value(UQ_PROPAG_PROBLEM_DATA_OUTPUT_ALLOW_ODV    ), "subEnvs that will write to data output file" )
#ifdef UQ_PROPAG_PROBLEM_READS_SOLVER_OPTION
    (m_option_solver.c_str(),               po::value<std::string>()->default_value(UQ_PROPAG_PROBLEM_SOLVER_ODV               ), "algorithm for propagation"                   )
#endif
  ;

  return;
}

template<class P_V,class P_M,class Q_V,class Q_M>
void
uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M>::getMyOptionValues(
  po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << optionsDesc
                              << std::endl;
    }
  }

  if (m_env.allOptionsMap().count(m_option_computeSolution.c_str())) {
    m_computeSolution = ((const po::variable_value&) m_env.allOptionsMap()[m_option_computeSolution.c_str()]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_computeCovariances.c_str())) {
    m_computeCovariances = ((const po::variable_value&) m_env.allOptionsMap()[m_option_computeCovariances.c_str()]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_computeCorrelations.c_str())) {
    m_computeCorrelations = ((const po::variable_value&) m_env.allOptionsMap()[m_option_computeCorrelations.c_str()]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_dataOutputFileName.c_str())) {
    m_dataOutputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_dataOutputFileName.c_str()]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_dataOutputAllowedSet.c_str())) {
    m_dataOutputAllowedSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_dataOutputAllowedSet.c_str()].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpAllow);

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_dataOutputAllowedSet.insert((unsigned int) tmpAllow[i]);
      }
    }
  }

#ifdef UQ_PROPAG_PROBLEM_READS_SOLVER_OPTION
  if (m_env.allOptionsMap().count(m_option_solver.c_str())) {
    m_solverString = ((const po::variable_value&) m_env.allOptionsMap()[m_option_solver.c_str()]).as<std::string>();
  }
#endif

  return;
}

template <class P_V,class P_M, class Q_V, class Q_M>
bool
  uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M>::computeSolutionFlag() const
{
  return m_computeSolution;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M>::solveWithMonteCarlo()
{
  m_env.fullComm().Barrier();
  m_env.syncPrintDebugMsg("Entering uqStatisticalForwardProblemClass<P_V,P_M>::solveWithMonteCarlo()",1,3000000,m_env.fullComm());

  if (m_computeSolution == false) {
    if ((m_env.subDisplayFile())) {
      *m_env.subDisplayFile() << "In uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M>::solveWithMonteCarlo()"
                              << ": avoiding solution, as requested by user"
                              << std::endl;
    }
    return;
  }
  if ((m_env.subDisplayFile())) {
    *m_env.subDisplayFile() << "In uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M>::solveWithMonteCarlo()"
                            << ": computing solution, as requested by user"
                            << std::endl;
  }

  if (m_solutionPdf       ) delete m_solutionPdf;

  if (m_unifiedSolutionCdf) delete m_unifiedSolutionCdf;
  if (m_unifiedCdfValues  ) delete m_unifiedCdfValues;
  if (m_unifiedCdfGrids   ) delete m_unifiedCdfGrids;

  if (m_subSolutionCdf    ) delete m_subSolutionCdf;
  if (m_subCdfValues      ) delete m_subCdfValues;
  if (m_subCdfGrids       ) delete m_subCdfGrids;

  if (m_subSolutionMdf    ) delete m_subSolutionMdf;
  if (m_subMdfValues      ) delete m_subMdfValues;
  if (m_subMdfGrids       ) delete m_subMdfGrids;

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
  m_paramChain = new uqSequenceOfVectorsClass<P_V,P_M>(m_paramRv.imageSet().vectorSpace(),0,m_prefix+"paramChain");
  m_qoiChain   = new uqSequenceOfVectorsClass<Q_V,Q_M>(m_qoiRv.imageSet().vectorSpace(),  0,m_prefix+"qoiChain"  );
  m_mcSeqGenerator = new uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>(m_prefix.c_str(),
                                                              m_paramRv,
                                                              m_qoiFunction,
                                                              m_qoiRv);
  m_mcSeqGenerator->generateSequence(*m_paramChain,
                                     *m_qoiChain);
  m_solutionRealizer = new uqSequentialVectorRealizerClass<Q_V,Q_M>((m_prefix+"Qoi").c_str(),
                                                                    *m_qoiChain);
  m_qoiRv.setRealizer(*m_solutionRealizer);

  // Compute output mdf: uniform sampling approach
  m_subMdfGrids  = new uqArrayOfOneDGridsClass <Q_V,Q_M>((m_prefix+"QoiMdf_").c_str(),m_qoiRv.imageSet().vectorSpace());
  m_subMdfValues = new uqArrayOfOneDTablesClass<Q_V,Q_M>((m_prefix+"QoiMdf_").c_str(),m_qoiRv.imageSet().vectorSpace());
  m_qoiChain->subUniformlySampledMdf(numEvaluationPointsVec, // input
                                     *m_subMdfGrids,         // output
                                     *m_subMdfValues);       // output

  m_subSolutionMdf = new uqSampledVectorMdfClass<Q_V,Q_M>((m_prefix+"Qoi").c_str(),
                                                          *m_subMdfGrids,
                                                          *m_subMdfValues);
  m_qoiRv.setMdf(*m_subSolutionMdf);

  // Compute output cdf: uniform sampling approach
  std::string subCoreName_qoiCdf(m_prefix+    "QoiCdf_");
  std::string uniCoreName_qoiCdf(m_prefix+"unifQoiCdf_");
  if (m_env.numSubEnvironments() == 1) subCoreName_qoiCdf = uniCoreName_qoiCdf;

  std::string subCoreName_solutionCdf(m_prefix+    "Qoi");
  std::string uniCoreName_solutionCdf(m_prefix+"unifQoi");
  if (m_env.numSubEnvironments() == 1) subCoreName_solutionCdf = uniCoreName_solutionCdf;

  m_subCdfGrids  = new uqArrayOfOneDGridsClass <Q_V,Q_M>(subCoreName_qoiCdf.c_str(),m_qoiRv.imageSet().vectorSpace());
  m_subCdfValues = new uqArrayOfOneDTablesClass<Q_V,Q_M>(subCoreName_qoiCdf.c_str(),m_qoiRv.imageSet().vectorSpace());
  m_qoiChain->subUniformlySampledCdf(numEvaluationPointsVec, // input
                                     *m_subCdfGrids,         // output
                                     *m_subCdfValues);       // output

  m_subSolutionCdf = new uqSampledVectorCdfClass<Q_V,Q_M>(subCoreName_solutionCdf.c_str(),
                                                          *m_subCdfGrids,
                                                          *m_subCdfValues);
  m_qoiRv.setSubCdf(*m_subSolutionCdf);

  // Compute unified cdf if necessary
  if (m_env.numSubEnvironments() == 1) {
    m_qoiRv.setUnifiedCdf(*m_subSolutionCdf);
  }
  else {
    m_unifiedCdfGrids  = new uqArrayOfOneDGridsClass <Q_V,Q_M>(uniCoreName_qoiCdf.c_str(),m_qoiRv.imageSet().vectorSpace());
    m_unifiedCdfValues = new uqArrayOfOneDTablesClass<Q_V,Q_M>(uniCoreName_qoiCdf.c_str(),m_qoiRv.imageSet().vectorSpace());
    m_qoiChain->unifiedUniformlySampledCdf(numEvaluationPointsVec, // input
                                           *m_unifiedCdfGrids,     // output
                                           *m_unifiedCdfValues);   // output

    m_unifiedSolutionCdf = new uqSampledVectorCdfClass<Q_V,Q_M>(uniCoreName_solutionCdf.c_str(),
                                                                *m_unifiedCdfGrids,
                                                                *m_unifiedCdfValues);
    m_qoiRv.setUnifiedCdf(*m_unifiedSolutionCdf);
  }

  // Compute (just unified one) covariance matrix, if requested
  // Compute (just unified one) correlation matrix, if requested
  P_M* pqCovarianceMatrix  = NULL;
  P_M* pqCorrelationMatrix = NULL;
  if (m_computeCovariances || m_computeCorrelations) {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M>::solveWithMonteCarlo()"
                             << ", prefix = " << m_prefix
                             << ": instantiating cov and corr matrices"
                             << std::endl;
    }
    pqCovarianceMatrix = new P_M(m_env,
                                 m_paramRv.imageSet().vectorSpace().map(),       // number of rows
                                 m_qoiRv.imageSet().vectorSpace().dimGlobal());  // number of cols
    pqCorrelationMatrix = new P_M(m_env,
                                  m_paramRv.imageSet().vectorSpace().map(),      // number of rows
                                  m_qoiRv.imageSet().vectorSpace().dimGlobal()); // number of cols
#if 0
    uqComputeCovCorrMatricesBetweenVectorRvs<P_V,P_M,Q_V,Q_M>(m_paramRv,
                                                              m_qoiRv,
                                                              std::min(m_paramRv.realizer().subPeriod(),m_qoiRv.realizer().subPeriod()), // FIX ME: might be INFINITY
                                                              *pqCovarianceMatrix,
                                                              *pqCorrelationMatrix);
#else
    uqComputeCovCorrMatricesBetweenVectorSequences(*m_paramChain,
                                                   *m_qoiChain,
                                                   std::min(m_paramRv.realizer().subPeriod(),m_qoiRv.realizer().subPeriod()), // FIX ME: might be INFINITY
                                                   *pqCovarianceMatrix,
                                                   *pqCorrelationMatrix);
#endif
  }

  // Open data output file
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "In uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M>::solveWithMonteCarlo()"
                           << ", prefix = "                                        << m_prefix
                           << ": checking necessity of opening data output file '" << m_dataOutputFileName
                           << "'"
                           << std::endl;
  }
  std::ofstream* ofsvar = NULL;
  m_env.openOutputFile(m_dataOutputFileName,
                       "m",
                       m_dataOutputAllowedSet,
                       false,
                       ofsvar);

  // Write data out
  if (m_env.subDisplayFile()) {
    if (pqCovarianceMatrix ) {
      *m_env.subDisplayFile() << "In uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M>::solveWithMonteCarlo()"
                             << ", prefix = "                           << m_prefix
                             << ": contents of covariance matrix are\n" << *pqCovarianceMatrix
                             << std::endl;
    }
    if (pqCorrelationMatrix) {
      *m_env.subDisplayFile() << "In uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M>::solveWithMonteCarlo()"
                             << ", prefix = "                            << m_prefix
                             << ": contents of correlation matrix are\n" << *pqCorrelationMatrix
                             << std::endl;
    }
  }

  if (ofsvar) {
    m_qoiRv.mdf().print(*ofsvar);
    *ofsvar << m_qoiRv.subCdf();

    //if (pqCovarianceMatrix ) *ofsvar << *pqCovarianceMatrix;  // FIX ME: output matrix in matlab format
    //if (pqCorrelationMatrix) *ofsvar << *pqCorrelationMatrix; // FIX ME: output matrix in matlab format

    // Write unified cdf if necessary
    if (m_env.numSubEnvironments() > 1) {
      if (m_qoiRv.imageSet().vectorSpace().zeroVector().numberOfProcessorsRequiredForStorage() == 1) {
        if (m_env.inter0Rank() == 0) {
          *ofsvar << m_qoiRv.unifiedCdf(); //*m_unifiedSolutionCdf;
        }
      }
      else {
        UQ_FATAL_TEST_MACRO(true,
                            m_env.fullRank(),
                            "uqStatisticalForwardProblem<P_V,P_M,Q_V,Q_M>::solveWithMonteCarlo()",
                            "unifed cdf writing, parallel vectors not supported yet");
      }
    }
  }

  // Close data output file
  if (ofsvar) {
    ofsvar->close();
    delete ofsvar;
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M>::solveWithMonteCarlo()"
                             << ", prefix = "                 << m_prefix
                             << ": closed data output file '" << m_dataOutputFileName
                             << "'"
                             << std::endl;
    }
  }
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << std::endl;
  }

  if (pqCovarianceMatrix ) delete pqCovarianceMatrix;
  if (pqCorrelationMatrix) delete pqCorrelationMatrix;

  m_env.syncPrintDebugMsg("Leaving uqStatisticalForwardProblemClass<P_V,P_M>::solveWithMonteCarlo()",1,3000000,m_env.fullComm());
  m_env.fullComm().Barrier();

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
const uqGenericVectorRVClass<Q_V,Q_M>& 
uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M>::qoiRv() const
{
  return m_qoiRv;
}

template <class P_V,class P_M,class Q_V,class Q_M>
const uqBaseVectorCdfClass<Q_V,Q_M>&
uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M>::qoiRv_unifiedCdf() const
{
  if (m_env.numSubEnvironments() == 1) {
    return m_qoiRv.subCdf();
  }

  if (m_env.inter0Rank() < 0) {
    return m_qoiRv.subCdf();
  }

  //UQ_FATAL_TEST_MACRO(m_unifiedSolutionCdf == NULL,
  //                    m_env.fullRank(),
  //                    "uqStatisticalForwardProblem<P_V,P_M,Q_V,Q_M>::qoiRv_unifiedCdf()",
  //                    "variable is NULL");
  return m_qoiRv.unifiedCdf(); //*m_unifiedSolutionCdf;
}


template <class P_V,class P_M,class Q_V,class Q_M>
void
uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M>::print(std::ostream& os) const
{
  os <<         m_option_computeSolution      << " = " << m_computeSolution
     << "\n" << m_option_computeCovariances   << " = " << m_computeCovariances
     << "\n" << m_option_computeCorrelations  << " = " << m_computeCorrelations
     << "\n" << m_option_dataOutputFileName   << " = " << m_dataOutputFileName;
  os << "\n" << m_option_dataOutputAllowedSet << " = ";
  for (std::set<unsigned int>::iterator setIt = m_dataOutputAllowedSet.begin(); setIt != m_dataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
#ifdef UQ_PROPAG_PROBLEM_READS_SOLVER_OPTION
       << "\n" << m_option_solver << " = " << m_solverString
#endif
  os << std::endl;
}

template<class P_V,class P_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_PROPAG_PROBLEM_H__
