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

#ifndef __UQ_SFP_H__
#define __UQ_SFP_H__

#include <queso/StatisticalForwardProblemOptions.h>
#include <queso/VectorFunction.h>
#include <queso/MonteCarloSG.h>
#include <queso/VectorRV.h>
#include <queso/SequenceOfVectors.h>

namespace QUESO {

/*! \file uqStatisticalForwardProblem.h
    \brief Class to solve a Statistical Forward Problem
*/

/*! \class StatisticalForwardProblem
 *  \brief This templated class represents a Statistical Forward Problem.
 * 
 * This templated class represents a statistical forward problem. It is templated on the types
 * 'P_V' and 'Q_V' of vectors and types 'P_M' and 'Q_M' of matrices, where 'P_' stands for 
 * 'parameter' and 'Q_' stands for 'quantities of interest'. Conceptually, a statistical forward 
 * problem has two input entities and one output entity.\n
 *
 * The input entities of a statistical forward problem are:
<list type=number>
<item> the input (parameter) RV, an instance of class 'BaseVectorRV<P_V,P_M>', and
<item> the QoI function, an instance of class 'BaseVectorFunction<P_V,P_M,Q_V,Q_M>'.
</list>
  * Let \f$ q() \f$  denote the mathematical QoI function and  \f$ x \f$  denote a vector of 
  * parameters. The QoI function object stores the routine that computes  \f$ q(x) \f$  and 
  * whatever data necessary by such routine. See file 'libs/basic/inc/VectorFunction.h' 
  * for more details.\n
  *
  * The output entity of a statistical forward problem is:
<list type=number>
<item> the QoI RV, another instance of class 'BaseVectorRV<P_V,P_M>'.
</list>   
  * The QoI RV stores the solution according to the Bayesian approach. The solution of a SPF 
  * is computed by calling 'solveWithMonteCarlo()'. More operations, with different methods, 
  * will be available in the future.\n
  *
  * The solution process might demand extra objects to be passed through the chosen solution 
  * operation interface. This distinction is important: this class separates 'what the problem 
  * is' from 'how the problem is solved'. Upon return from a solution operation, the QoI RV is 
  * available through the operation 'qoiRv()'. Such QoI RV is able to provide:
<list type=number>
<item> a vector realizer through the operation 'qoiRv().realizer()', which returns an
       instance of the class 'BaseVectorRealizer<Q_V,Q_M>'.
</list>*/

template <class P_V,class P_M,class Q_V,class Q_M>
class StatisticalForwardProblem
{
public:
 //! @name Constructor/Destructor methods
 //@{
 //! Constructor.
 /*! Requirements: 1) the image set of the vector random variable 'paramRv' and the domain set of 
  * the QoI function 'qoiFunction' should belong to vector spaces of equal dimensions and 2) the 
  * image set of the QoI function 'qoiFunction' and the image set of the vector random variable 
  * 'qoiRv' should belong to vector spaces of equal dimensions. If the requirements are satisfied, 
  * the constructor then reads input options that begin with the string '\<prefix\>fp_'. Options 
  * reading is handled by class 'StatisticalForwardProblemOptions'. If no options input file 
  * is provided, the construction assigns \c alternativeOptionsValues to the options of the SFP.*/
  StatisticalForwardProblem(const char*                                       prefix,
                                   const SfpOptionsValues*                    alternativeOptionsValues, // dakota
                                   const BaseVectorRV      <P_V,P_M>&         paramRv,
                                   const BaseVectorFunction<P_V,P_M,Q_V,Q_M>& qoiFunction,
                                   GenericVectorRV         <Q_V,Q_M>&         qoiRv);
 
  //! Destructor
  ~StatisticalForwardProblem();
  //@}

  //! @name Statistical methods
  //@{
  //! Whether or not compute the solution.  
  bool                                   computeSolutionFlag() const;
  
  //! Solves the problem through Monte Carlo algorithm.
  /*! Requirements: none at this moment. If the requirements are satisfied, this operation checks
   * the member flag 'm_computeSolution' (one of the options read from the input file during 
   * construction). If the flag is 'false', the operation returns immediately, computing nothing.
   * If the flag is 'true', the operation sets the member variable 'm_qoiRv' accordingly. 
   * The operation:
<list type=number>
<item> instantiates 'SequenceOfVectors<P_V,P_M>' (the input sequence of vectors),
<item> instantiates 'SequenceOfVectors<Q_V,Q_M>' (the output sequence of vectors),
<item> instantiates 'MonteCarloSG<P_V,P_M,Q_V,Q_M>' (the Monte Carlo algorithm),
<item> populates the output sequence with the Monte Carlo algorithm,
<item> sets the realizer of 'm_qoiRv' with the contents of the output sequence, and
</list>*/
//<item> computes the CDFs of the components of 'm_qoiRv' as instances of 'SampledVectorCdf<Q_V,Q_M>'
  void                                   solveWithMonteCarlo(const McOptionsValues* alternativeOptionsValues); // dakota
  
  //! Returns the QoI RV; access to private attribute m_qoiRv.
  const GenericVectorRV<Q_V,Q_M>& qoiRv              () const;

  //! Returns the parameter chain; access to private attribute m_paramChain.
  const BaseVectorSequence<Q_V,Q_M>& getParamChain   () const;
  
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
  
  //<item> CDFs of QoI components through the operation 'qoiRv().unifiedCdf()',
  //     which returns an instance of the class 'BaseVectorCdf<Q_V,Q_M>'
  const BaseVectorCdf  <Q_V,Q_M>& qoiRv_unifiedCdf   () const;
  
#endif
  //@}
  
  //! @name I/O methods
  //@{
  //! TODO: Prints the sequence.
  /*! \todo: implement me!*/
  void print(std::ostream& os) const;
  //@}

private:
  //! TODO: Common constructor.
  /*! \todo: implement me!*/
  void  commonConstructor();

  const BaseEnvironment&                     m_env;

  const BaseVectorRV      <P_V,P_M>&         m_paramRv;
  const BaseVectorFunction<P_V,P_M,Q_V,Q_M>& m_qoiFunction;
        GenericVectorRV   <Q_V,Q_M>&         m_qoiRv; // FIX ME: Maybe not always generic ?

        BaseVectorSequence<Q_V,Q_M>*         m_paramChain;
        BaseVectorSequence<Q_V,Q_M>*         m_qoiChain;
        MonteCarloSG      <P_V,P_M,Q_V,Q_M>* m_mcSeqGenerator;

        BaseVectorRealizer<Q_V,Q_M>*         m_solutionRealizer;
 
        BaseJointPdf      <Q_V,Q_M>*         m_solutionPdf;

        SfpOptionsValues                     m_alternativeOptionsValues;
        StatisticalForwardProblemOptions*    m_optionsObj;

#ifdef UQ_ALSO_COMPUTE_MDFS_WITHOUT_KDE
        ArrayOfOneDGrids  <Q_V,Q_M>*         m_subMdfGrids;
        ArrayOfOneDTables <Q_V,Q_M>*         m_subMdfValues;
#endif
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
        BaseVectorMdf     <Q_V,Q_M>*         m_subSolutionMdf;
        ArrayOfOneDGrids  <Q_V,Q_M>*         m_subCdfGrids;
        ArrayOfOneDTables <Q_V,Q_M>*         m_subCdfValues;
        BaseVectorCdf     <Q_V,Q_M>*         m_subSolutionCdf;

        ArrayOfOneDGrids  <Q_V,Q_M>*         m_unifiedCdfGrids;
        ArrayOfOneDTables <Q_V,Q_M>*         m_unifiedCdfValues;
        BaseVectorCdf     <Q_V,Q_M>*         m_unifiedSolutionCdf;
#endif
};
//! Prints the object \c obj, overloading an operator.
template<class P_V,class P_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>& obj);

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
  m_alternativeOptionsValues(),
  m_optionsObj              (NULL)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = " << prefix
                            << ", alternativeOptionsValues = " << alternativeOptionsValues
                            << ", m_env.optionsInputFileName() = " << m_env.optionsInputFileName()
                            << std::endl;
  }

  if (alternativeOptionsValues) m_alternativeOptionsValues = *alternativeOptionsValues;
  if (m_env.optionsInputFileName() == "") {
    m_optionsObj = new StatisticalForwardProblemOptions(m_env,prefix,m_alternativeOptionsValues);
  }
  else {
    m_optionsObj = new StatisticalForwardProblemOptions(m_env,prefix);
    m_optionsObj->scanOptionsValues();
  }

  UQ_FATAL_TEST_MACRO(paramRv.imageSet().vectorSpace().dimLocal() != qoiFunction.domainSet().vectorSpace().dimLocal(),
                      m_env.worldRank(),
                      "StatisticalForwardProblem<P_V,P_M>::constructor()",
                      "'paramRv' and 'qoiFunction' are related to vector spaces of different dimensions");

  UQ_FATAL_TEST_MACRO(qoiFunction.imageSet().vectorSpace().dimLocal() != qoiRv.imageSet().vectorSpace().dimLocal(),
                      m_env.worldRank(),
                      "StatisticalForwardProblem<P_V,P_M>::constructor()",
                      "'qoiFunction' and 'qoiRv' are related to vector spaces of different dimensions");

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
  return m_optionsObj->m_ov.m_computeSolution;
}
//--------------------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
void
StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>::solveWithMonteCarlo(
  const McOptionsValues* alternativeOptionsValues)
{
  m_env.fullComm().Barrier();
  m_env.fullComm().syncPrintDebugMsg("Entering StatisticalForwardProblem<P_V,P_M>::solveWithMonteCarlo()",1,3000000);

  if (m_optionsObj->m_ov.m_computeSolution == false) {
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
  if (m_optionsObj->m_ov.m_computeCovariances || m_optionsObj->m_ov.m_computeCorrelations) {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>::solveWithMonteCarlo()"
                              << ", prefix = " << m_optionsObj->m_prefix
                              << ": instantiating cov and corr matrices"
                              << std::endl;
    }
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
                            << ": checking necessity of opening data output file '" << m_optionsObj->m_ov.m_dataOutputFileName
                            << "'"
                            << std::endl;
  }
  FilePtrSetStruct filePtrSet;
  if (m_env.openOutputFile(m_optionsObj->m_ov.m_dataOutputFileName,
                           UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT, // Yes, always ".m"
                           m_optionsObj->m_ov.m_dataOutputAllowedSet,
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
        UQ_FATAL_TEST_MACRO(true,
                            m_env.worldRank(),
                            "StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>::solveWithMonteCarlo()",
                            "unified cdf writing, parallel vectors not supported yet");
      }
    }
#endif
    // Close data output file
    m_env.closeFile(filePtrSet,UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>::solveWithMonteCarlo()"
                              << ", prefix = "                 << m_optionsObj->m_prefix
                              << ": closed data output file '" << m_optionsObj->m_ov.m_dataOutputFileName
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

  //UQ_FATAL_TEST_MACRO(m_unifiedSolutionCdf == NULL,
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
  UQ_FATAL_TEST_MACRO(m_paramChain == NULL,
                      m_env.worldRank(),
                      (std::string)("StatisticalForwardProblem<V,M,V,M>::getParamChain()"),
                      "m_paramChain is NULL");

  return *m_paramChain;
 
}
// I/O methods--------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
void
StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>::print(std::ostream& os) const
{
  return;
}
//--------------------------------------------------
// Operator declared outside class definition ------
//--------------------------------------------------
template<class P_V,class P_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>& obj)
{
  obj.print(os);

  return os;
}

}  // End namespace QUESO

#endif // __UQ_SFP_H__
