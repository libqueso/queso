//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012 The PECOS Development Team
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

#ifndef __UQ_SIP_H__
#define __UQ_SIP_H__

#include <uqStatisticalInverseProblemOptions.h>
#include <uqMetropolisHastingsSG1.h>
#include <uqMLSampling1.h>
#include <uqInstantiateIntersection.h>
#include <uqVectorRV.h>
#include <uqScalarFunction.h>
//#include <grvy.h>

/*! This templated class represents a statistical inverse problem.
    It is templated on the type 'P_V' of vector and type 'P_M' of matrix,
    where 'P_' stands for 'parameter'.
    Some templated classes might also use 'Q_V' and 'Q_M' when referring to a vector type and a matrix type,
    where 'Q_' stands for 'quantities of interest'.
    See e.g. the class 'uqStatisticalForwardProblemClass'.

    Conceptually, a statistical inverse problem has two input entities and one output entity.
*/
/*! -------------------------------------------------------------
*/
/*! The input entities of a statistical inverse problem are:
<list type=number>
<item> the prior rv, an instance of class 'uqBaseVectorRVClass<P_V,P_M>', and
<item> the likelihood function, an instance of class 'uqBaseScalarFunctionClass<P_V,P_M>'.
</list>
    Let '\pi(.)' denote the mathematical likelihood function and 'x' denote a vector of parameters.
    The likelihood function object stores the routine that computes \pi(x) and whatever data necessary by such routine.
    The routine in the likelihood function object can compute either the actual value \pi(x) or the value ln[\pi(x)].
    The constructor of the 'scalar function' asks for the user to specify which value the routine is actually computing,
    so that the 'scalar function' class can properly implement both class operations 'actualValue()' and 'minus2LnValue()'.
    See files 'libs/basic/inc/uqScalarFunction.h' and 'libs/stats/inc/uqJointPdf.h' for more details.
*/
/*! -------------------------------------------------------------
*/
/*! The output entity of a statistical inverse problem is:
<list type=number>
<item> the posterior rv, another instance of class 'uqBaseVectorRVClass<P_V,P_M>'.
</list>   
    The posterior rv stores the solution according to the Bayesian approach.
    A similar situation occurs e.g. in the case of a system Ax=b of linear equations,
    where 'A' and 'b' are inputs, and 'x' stores the solution of the inverse problem.
*/
/*! -------------------------------------------------------------
*/
/*! The solution of a statistical inverse problem is computed by calling one of the following operations:
<list type=number>
<item> 'solveWithBayesMetropolisHastings(...)'.
<item> 'solveWithBayesMLSampling()'.
</list> 
    More operations, with different methods, will be available in the future.
*/
/*! The solution process might demand extra objects to be passed through the chosen solution operation interface.
    This distinction is important: this class separates 'what the problem is' from 'how the problem is solved'.
*/
/*! -------------------------------------------------------------
*/
/*! Upon return from a solution operation, the posterior rv is available through
    the operation 'postRv()'. Such posterior rv is able to provide:
<list type=number>
<item> a joint pdf (up to a multiplicative constant) through the operation 'postRv().pdf()',
       which returns an instance of the class 'uqBaseJointPdfClass<P_V,P_M>', and
<item> a vector realizer through the operation 'postRv().realizer()', which returns an
       instance of the class 'uqBaseVectorRealizerClass<P_V,P_M>'.
</list>
*/
/*! -------------------------------------------------------------
*/
/*! If the options request data to be written in the output file (MATLAB .m format only, for now),
    the user can run 'grep zeros \<OUTPUT FILE NAME\>' after the solution procedure ends
    in order to check which MATLAB variables are defined and set.
    The names of the varibles are self explanatory.
*/
template <class P_V,class P_M>
class uqStatisticalInverseProblemClass
{
public:
  uqStatisticalInverseProblemClass(const char*                               prefix,
                                   const uqSipOptionsValuesClass*            alternativeOptionsValues, // dakota
                                   const uqBaseVectorRVClass      <P_V,P_M>& priorRv,            
                                   const uqBaseScalarFunctionClass<P_V,P_M>& likelihoodFunction, 
                                         uqGenericVectorRVClass   <P_V,P_M>& postRv);
 ~uqStatisticalInverseProblemClass();

        bool                             computeSolutionFlag             () const;
        void                             solveWithBayesMetropolisHastings(const uqMhOptionsValuesClass* alternativeOptionsValues, // dakota
                                                                          const P_V&                    initialValues,
                                                                          const P_M*                    initialProposalCovMatrix);
        void                             solveWithBayesMLSampling        ();
  const uqBaseVectorRVClass   <P_V,P_M>& priorRv                         () const;
  const uqGenericVectorRVClass<P_V,P_M>& postRv                          () const;
	double                           logEvidence                     () const;
	double                           meanLogLikelihood               () const;
	double                           eig                             () const;

        void                             print                           (std::ostream& os) const;

private:
  const uqBaseEnvironmentClass&                 m_env;

  const uqBaseVectorRVClass        <P_V,P_M>&   m_priorRv;
  const uqBaseScalarFunctionClass  <P_V,P_M>&   m_likelihoodFunction;
        uqGenericVectorRVClass     <P_V,P_M>&   m_postRv;

        uqVectorSetClass           <P_V,P_M>*   m_solutionDomain;
        uqBaseJointPdfClass        <P_V,P_M>*   m_solutionPdf;
        uqBaseVectorMdfClass       <P_V,P_M>*   m_subSolutionMdf;
        uqBaseVectorCdfClass       <P_V,P_M>*   m_subSolutionCdf;
        uqBaseVectorRealizerClass  <P_V,P_M>*   m_solutionRealizer;

        uqMetropolisHastingsSGClass<P_V,P_M>*   m_mhSeqGenerator;
        uqMLSamplingClass          <P_V,P_M>*   m_mlSampler;
        uqBaseVectorSequenceClass  <P_V,P_M>*   m_chain;
#ifdef UQ_ALSO_COMPUTE_MDFS_WITHOUT_KDE
        uqArrayOfOneDGridsClass    <P_V,P_M>*   m_subMdfGrids;
        uqArrayOfOneDTablesClass   <P_V,P_M>*   m_subMdfValues;
#endif

        uqSipOptionsValuesClass                  m_alternativeOptionsValues;
        uqStatisticalInverseProblemOptionsClass* m_optionsObj;
};

template<class P_V,class P_M>
std::ostream& operator<<(std::ostream& os, const uqStatisticalInverseProblemClass<P_V,P_M>& obj);

/*! Constructor. */
/*! Requirements:
<list type=number>
<item> the image set of the vector random variable 'priorRv',
       the domain set of the likelihood function 'likelihoodFunction' and
       the image set of the vector random variable 'postRv'
       should belong to vector spaces of equal dimensions.
</list>
*/
/*! If the requirements are satisfied, the constructor then reads input options that begin with the string '\<prefix\>ip_'.
    For instance, if 'prefix' is 'pROblem_775_', then the constructor will read all options that begin with 'pROblem_775_ip_'.
    Options reading is handled by class 'uqStatisticalInverseProblemOptionsClass'.
*/
/*! Input options are read from the QUESO input file, whose name is required by the constructor of the QUESO environment class.
    The QUESO environment class is instantiated at the application level, right after 'MPI_Init(&argc,&argv)'. 
    The QUESO environment is required by reference by many constructors in the QUESO library, and is available by reference from many classes as well.
*/
template <class P_V,class P_M>
uqStatisticalInverseProblemClass<P_V,P_M>::uqStatisticalInverseProblemClass(
  /*! The prefix                 */ const char*                               prefix,
  /*! Options (if no input file) */ const uqSipOptionsValuesClass*            alternativeOptionsValues, // dakota
  /*! The prior rv               */ const uqBaseVectorRVClass      <P_V,P_M>& priorRv,
  /*! The likelihood function    */ const uqBaseScalarFunctionClass<P_V,P_M>& likelihoodFunction,
  /*! The posterior rv           */       uqGenericVectorRVClass   <P_V,P_M>& postRv)
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
  m_alternativeOptionsValues(),
  m_optionsObj              (NULL)
{
#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "Entering uqSipClass" << std::endl;
#endif
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering uqStatisticalInverseProblemClass<P_V,P_M>::constructor()"
                            << ": prefix = " << prefix
                            << ", alternativeOptionsValues = " << alternativeOptionsValues
                            << ", m_env.optionsInputFileName() = " << m_env.optionsInputFileName()
                            << std::endl;
  }

  if (alternativeOptionsValues) m_alternativeOptionsValues = *alternativeOptionsValues;
  if (m_env.optionsInputFileName() == "") {
    m_optionsObj = new uqStatisticalInverseProblemOptionsClass(m_env,prefix,m_alternativeOptionsValues);
  }
  else {
    m_optionsObj = new uqStatisticalInverseProblemOptionsClass(m_env,prefix);
    m_optionsObj->scanOptionsValues();
  }
#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "In uqSipClass, finished scanning options" << std::endl;
#endif

  UQ_FATAL_TEST_MACRO(priorRv.imageSet().vectorSpace().dimLocal() != likelihoodFunction.domainSet().vectorSpace().dimLocal(),
                      m_env.worldRank(),
                      "uqStatisticalInverseProblemClass<P_V,P_M>::constructor()",
                      "'priorRv' and 'likelihoodFunction' are related to vector spaces of different dimensions");

  UQ_FATAL_TEST_MACRO(priorRv.imageSet().vectorSpace().dimLocal() != postRv.imageSet().vectorSpace().dimLocal(),
                      m_env.worldRank(),
                      "uqStatisticalInverseProblemClass<P_V,P_M>::constructor()",
                      "'priorRv' and 'postRv' are related to vector spaces of different dimensions");

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving uqStatisticalInverseProblemClass<P_V,P_M>::constructor()"
                            << ": prefix = " << m_optionsObj->m_prefix
                            << std::endl;
  }

  return;
}

/*! Destructor. */
template <class P_V,class P_M>
uqStatisticalInverseProblemClass<P_V,P_M>::~uqStatisticalInverseProblemClass()
{
  if (m_chain) {
    m_chain->clear();
    delete m_chain;
  }
  if (m_mlSampler       ) delete m_mlSampler;
  if (m_mhSeqGenerator  ) delete m_mhSeqGenerator;
  if (m_solutionRealizer) delete m_solutionRealizer;
  if (m_subSolutionCdf  ) delete m_subSolutionCdf;
  if (m_subSolutionMdf  ) delete m_subSolutionMdf;
  if (m_solutionPdf     ) delete m_solutionPdf;
  if (m_solutionDomain  ) delete m_solutionDomain;
  if (m_optionsObj      ) delete m_optionsObj;
}

/*! Operation to solve the problem through Bayes formula and a Metropolis-Hastings algorithm. */
/*! Requirements:
<list type=number>
<item> 'initialValues' should have the same number of components as member variable 'm_priorRv'
<item> if 'initialProposalCovMatrix' is not NULL, it should be square and its size should be equal to the size of 'initialValues'
</list>
*/
/*! If the requirements are satisfied, this operation checks the member flag 'm_computeSolution' (one of the options read from the input file during construction).
*/
/*! If the flag is 'false', the operation returns immediately, computing nothing.
 */
/*! If the flag is 'true', the operation sets the member variable 'm_postRv' accordingly. The operation:
<list type=number>
<item> sets the pdf of 'm_postRv' equal to an instance of 'uqBayesianJointPdfClass<P_V,P_M>',
<item> instantiates 'uqSequenceOfVectorsClass<P_V,P_M>' (the chain),
<item> instantiates 'uqMetropolisHastingsSGClass<P_V,P_M>' (tha Metropolis-Hastings algorithm),
<item> populates the chain with the Metropolis-Hastings algorithm, and
<item> sets the realizer of 'm_postRv' with the contents of the chain.
</list>
*/
template <class P_V,class P_M>
void
uqStatisticalInverseProblemClass<P_V,P_M>::solveWithBayesMetropolisHastings(
  const uqMhOptionsValuesClass* alternativeOptionsValues, // dakota
  const P_V&                    initialValues,
  const P_M*                    initialProposalCovMatrix)
{
  //grvy_timer_begin("BayesMetropolisHastings"); TODO: revisit timing output
  //std::cout << "proc " << m_env.fullRank() << ", HERE sip 000" << std::endl;
  m_env.fullComm().Barrier();
  //std::cout << "proc " << m_env.fullRank() << ", HERE sip 001" << std::endl;
  m_env.fullComm().syncPrintDebugMsg("Entering uqStatisticalInverseProblemClass<P_V,P_M>::solveWithBayesMetropolisHastings()",1,3000000);

  if (m_optionsObj->m_ov.m_computeSolution == false) {
    if ((m_env.subDisplayFile())) {
      *m_env.subDisplayFile() << "In uqStatisticalInverseProblemClass<P_V,P_M>::solveWithBayesMetropolisHastings()"
                              << ": avoiding solution, as requested by user"
                              << std::endl;
    }
    return;
  }
  if ((m_env.subDisplayFile())) {
    *m_env.subDisplayFile() << "In uqStatisticalInverseProblemClass<P_V,P_M>::solveWithBayesMetropolisHastings()"
                            << ": computing solution, as requested by user"
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(m_priorRv.imageSet().vectorSpace().dimLocal() != initialValues.sizeLocal(),
                      m_env.worldRank(),
                      "uqStatisticalInverseProblemClass<P_V,P_M>::solveWithBayesMetropolisHastings()",
                      "'m_priorRv' and 'initialValues' should have equal dimensions");

  if (initialProposalCovMatrix) {
    UQ_FATAL_TEST_MACRO(m_priorRv.imageSet().vectorSpace().dimLocal() != initialProposalCovMatrix->numRowsLocal(),
                        m_env.worldRank(),
                        "uqStatisticalInverseProblemClass<P_V,P_M>::solveWithBayesMetropolisHastings()",
                        "'m_priorRv' and 'initialProposalCovMatrix' should have equal dimensions");
    UQ_FATAL_TEST_MACRO(initialProposalCovMatrix->numCols() != initialProposalCovMatrix->numRowsGlobal(),
                        m_env.worldRank(),
                        "uqStatisticalInverseProblemClass<P_V,P_M>::solveWithBayesMetropolisHastings()",
                        "'initialProposalCovMatrix' should be a square matrix");
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
  m_solutionDomain = uqInstantiateIntersection(m_priorRv.pdf().domainSet(),m_likelihoodFunction.domainSet());

  m_solutionPdf = new uqBayesianJointPdfClass<P_V,P_M>(m_optionsObj->m_prefix.c_str(),
                                                       m_priorRv.pdf(),
                                                       m_likelihoodFunction,
                                                       1.,
                                                       *m_solutionDomain);

  m_postRv.setPdf(*m_solutionPdf);

  // Compute output realizer: Metropolis-Hastings approach
  m_chain = new uqSequenceOfVectorsClass<P_V,P_M>(m_postRv.imageSet().vectorSpace(),0,m_optionsObj->m_prefix+"chain");
  m_mhSeqGenerator = new uqMetropolisHastingsSGClass<P_V,P_M>(m_optionsObj->m_prefix.c_str(), // dakota
                                                              alternativeOptionsValues,
                                                              m_postRv,
                                                              initialValues,
                                                              initialProposalCovMatrix);

  m_mhSeqGenerator->generateSequence(*m_chain,NULL,NULL);

  m_solutionRealizer = new uqSequentialVectorRealizerClass<P_V,P_M>(m_optionsObj->m_prefix.c_str(),
                                                                   *m_chain);

  m_postRv.setRealizer(*m_solutionRealizer);

  m_env.fullComm().syncPrintDebugMsg("In uqStatisticalInverseProblemClass<P_V,P_M>::solveWithBayesMetropolisHastings(), code place 1",3,3000000);
  //m_env.fullComm().Barrier();

  // Compute output mdf: uniform sampling approach
#ifdef UQ_ALSO_COMPUTE_MDFS_WITHOUT_KDE
  m_subMdfGrids  = new uqArrayOfOneDGridsClass <P_V,P_M>((m_optionsObj->m_prefix+"Mdf_").c_str(),m_postRv.imageSet().vectorSpace());
  m_subMdfValues = new uqArrayOfOneDTablesClass<P_V,P_M>((m_optionsObj->m_prefix+"Mdf_").c_str(),m_postRv.imageSet().vectorSpace());
  m_chain->subUniformlySampledMdf(numEvaluationPointsVec, // input
                                  *m_subMdfGrids,         // output
                                  *m_subMdfValues);       // output
  m_subSolutionMdf = new uqSampledVectorMdfClass<P_V,P_M>(m_optionsObj->m_prefix.c_str(),
                                                          *m_subMdfGrids,
                                                          *m_subMdfValues);
  m_postRv.setMdf(*m_subSolutionMdf);

  if ((m_optionsObj->m_ov.m_dataOutputFileName                       != UQ_SIP_FILENAME_FOR_NO_FILE                    ) &&
      (m_optionsObj->m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_ov.m_dataOutputAllowedSet.end())) {
    if (m_env.subRank() == 0) {
      // Write data output file
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "Opening data output file '" << m_optionsObj->m_ov.m_dataOutputFileName
                                << "' for calibration problem with problem with prefix = " << m_optionsObj->m_prefix
                                << std::endl;
      }

      // Open file
      // Always write at the end of an eventual pre-existing file
      std::ofstream* ofsvar = new std::ofstream((m_optionsObj->m_ov.m_dataOutputFileName+"_sub"+m_env.subIdString()+".m").c_str(), std::ofstream::out | std::ofstream::in | std::ofstream::ate);
      if ((ofsvar            == NULL ) ||
          (ofsvar->is_open() == false)) {
        delete ofsvar;
        ofsvar = new std::ofstream((m_optionsObj->m_ov.m_dataOutputFileName+"_sub"+m_env.subIdString()+".m").c_str(), std::ofstream::out | std::ofstream::trunc);
      }
      UQ_FATAL_TEST_MACRO((ofsvar && ofsvar->is_open()) == false,
                          m_env.worldRank(),
                          "uqStatisticalInverseProblem<P_V,P_M>::solveWithBayesMetropolisHastings()",
                          "failed to open file");

      m_postRv.mdf().print(*ofsvar);

      // Close file
      //ofsvar->close();
      delete ofsvar;
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "Closed data output file '" << m_optionsObj->m_ov.m_dataOutputFileName
                                << "' for calibration problem with problem with prefix = " << m_optionsObj->m_prefix
                                << std::endl;
      }
    }
  }
#endif
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << std::endl;
  }

  m_env.fullComm().syncPrintDebugMsg("Leaving uqStatisticalInverseProblemClass<P_V,P_M>::solveWithBayesMetropolisHastings()",1,3000000);
  m_env.fullComm().Barrier();
  // grvy_timer_end("BayesMetropolisHastings"); TODO: revist timers
  return;
}

template <class P_V,class P_M>
void
uqStatisticalInverseProblemClass<P_V,P_M>::solveWithBayesMLSampling()
{
  m_env.fullComm().Barrier();
  m_env.fullComm().syncPrintDebugMsg("Entering uqStatisticalInverseProblemClass<P_V,P_M>::solveWithBayesMLSampling()",1,3000000);

  if (m_optionsObj->m_ov.m_computeSolution == false) {
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
  if (m_mhSeqGenerator  ) delete m_mhSeqGenerator;
  if (m_solutionRealizer) delete m_solutionRealizer;
  if (m_subSolutionCdf  ) delete m_subSolutionCdf;
  if (m_subSolutionMdf  ) delete m_subSolutionMdf;
  if (m_solutionPdf     ) delete m_solutionPdf;
  if (m_solutionDomain  ) delete m_solutionDomain;

  P_V numEvaluationPointsVec(m_priorRv.imageSet().vectorSpace().zeroVector());
  numEvaluationPointsVec.cwSet(250.);

  // Compute output pdf up to a multiplicative constant: Bayesian approach
  m_solutionDomain = uqInstantiateIntersection(m_priorRv.pdf().domainSet(),m_likelihoodFunction.domainSet());

  m_solutionPdf = new uqBayesianJointPdfClass<P_V,P_M>(m_optionsObj->m_prefix.c_str(),
                                                       m_priorRv.pdf(),
                                                       m_likelihoodFunction,
                                                       1.,
                                                       *m_solutionDomain);

  m_postRv.setPdf(*m_solutionPdf);

  // Compute output realizer: ML approach
  m_chain = new uqSequenceOfVectorsClass<P_V,P_M>(m_postRv.imageSet().vectorSpace(),0,m_optionsObj->m_prefix+"chain");
  m_mlSampler = new uqMLSamplingClass<P_V,P_M>(m_optionsObj->m_prefix.c_str(),
                                             //m_postRv,
                                               m_priorRv,
                                               m_likelihoodFunction);
  //                                           initialValues,
  //                                           initialProposalCovMatrix);

  m_mlSampler->generateSequence(*m_chain,NULL,NULL);

  m_solutionRealizer = new uqSequentialVectorRealizerClass<P_V,P_M>(m_optionsObj->m_prefix.c_str(),
                                                                    *m_chain);

  m_postRv.setRealizer(*m_solutionRealizer);

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << std::endl;
  }

  m_env.fullComm().syncPrintDebugMsg("Leaving uqStatisticalInverseProblemClass<P_V,P_M>::solveWithBayesMLSampling()",1,3000000);
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
double uqStatisticalInverseProblemClass<P_V,P_M>::logEvidence() const
{
  UQ_FATAL_TEST_MACRO(m_mlSampler == NULL,
                      m_env.worldRank(),
                      "uqStatisticalInverseProblemClass<P_V,P_M>::logEvidence()",
                      "m_mlSampler is NULL");
  return m_mlSampler->logEvidence();
}

template <class P_V,class P_M>
double uqStatisticalInverseProblemClass<P_V,P_M>::meanLogLikelihood() const
{
  UQ_FATAL_TEST_MACRO(m_mlSampler == NULL,
                      m_env.worldRank(),
                      "uqStatisticalInverseProblemClass<P_V,P_M>::meanLogLikelihood()",
                      "m_mlSampler is NULL");
  return m_mlSampler->meanLogLikelihood();
}

template <class P_V,class P_M>
double uqStatisticalInverseProblemClass<P_V,P_M>::eig() const
{
  UQ_FATAL_TEST_MACRO(m_mlSampler == NULL,
                      m_env.worldRank(),
                      "uqStatisticalInverseProblemClass<P_V,P_M>::eig()",
                      "m_mlSampler is NULL");
  return m_mlSampler->eig();
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
