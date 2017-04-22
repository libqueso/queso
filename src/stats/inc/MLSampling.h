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

#ifndef UQ_MULTI_LEVEL_SAMPLING_H
#define UQ_MULTI_LEVEL_SAMPLING_H

#define ML_NEW_CODE_2009_12_29

#include <queso/MLSamplingOptions.h>
#include <queso/MetropolisHastingsSG.h>
#include <queso/FiniteDistribution.h>
#include <queso/VectorRV.h>
#include <queso/GenericVectorRV.h>
#include <queso/VectorSpace.h>
#include <queso/MarkovChainPositionData.h>
#include <queso/ScalarFunctionSynchronizer.h>
#include <queso/SequenceOfVectors.h>
#include <queso/ArrayOfSequences.h>
#ifdef QUESO_HAS_GLPK
#include <glpk.h>
#endif
#include <sys/time.h>
#include <fstream>

#define ML_CHECKPOINT_FIXED_AMOUNT_OF_DATA 6

//---------------------------------------------------------

namespace QUESO {

class GslVector;
class GslMatrix;

#ifdef QUESO_HAS_GLPK
struct BIP_routine_struct {
  const BaseEnvironment* env;
  unsigned int                  currLevel;
};

void BIP_routine(glp_tree *tree, void *info);
#endif

//---------------------------------------------------------

struct ExchangeInfoStruct
{
  int          originalNodeOfInitialPosition;
  unsigned int originalIndexOfInitialPosition;
  int          finalNodeOfInitialPosition;
  unsigned int numberOfPositions;
};

//---------------------------------------------------------

template <class P_V = GslVector>
struct BalancedLinkedChainControlStruct
{
  P_V*         initialPosition;
  double       initialLogPrior;
  double       initialLogLikelihood;
  unsigned int numberOfPositions;
};

template <class P_V = GslVector>
struct BalancedLinkedChainsPerNodeStruct
{
  std::vector<BalancedLinkedChainControlStruct<P_V> > balLinkedChains;
};

//---------------------------------------------------------

struct UnbalancedLinkedChainControlStruct
{
  unsigned int initialPositionIndexInPreviousChain; // KEY new: add logPrior and logLikelihood
  unsigned int numberOfPositions;
};

struct UnbalancedLinkedChainsPerNodeStruct
{
  std::vector<UnbalancedLinkedChainControlStruct> unbLinkedChains;
};

//---------------------------------------------------------

//---------------------------------------------------------
 /*! \file MLSampling.h
 * \class MLSampling
 * \brief A templated class that represents a Multilevel generator of samples.
 *
 * A templated class that represents a Multilevel sampler.  Options reading is handled by class
 * 'MLSamplingOptions'.
 * It implements the method: S. H. Cheung and E. E. Prudencio. Parallel adaptive multilevel sampling
 * algorithms for the Bayesian analysis of mathematical models. International Journal for Uncertainty
 * Quantification, 2(3):215237, 2012.  */
 /* If options request data to be written in the output file (MATLAB .m format
 * only, for now), the user can check which MATLAB variables are defined and set by running
 * 'grep zeros <OUTPUT FILE NAME>' after the solution procedures ends.
 * The names of the variables are chosen to be self explanatory. */

template <class P_V = GslVector, class P_M = GslMatrix>
class MLSampling
{
public:
 //! @name Constructor/Destructor methods
 //@{
 //! Constructor
 /*! Constructor: instantiates an object of the class given a prefix, the prior RV and the likelihood function.
  * Internally, this method set the domain of the target PDF to be the intersection of the domains of the prior
  * PDF and of the likelihood function. Also the current level and exponent of the Multilevel method are set to
  * zero. */
  MLSampling(/*! Prefix                  */ const char*                               prefix,
                    /*! The prior rv            */ const BaseVectorRV      <P_V,P_M>& priorRv,
                    /*! The likelihood function */ const BaseScalarFunction<P_V,P_M>& likelihoodFunction);
 //! Destructor
 ~MLSampling();
 //@}
   //! @name Statistical methods
 //@{
 //! Method to generate the chain.
 /*! Requirement: the vector space 'm_vectorSpace' should have dimension equal to the size of a
 * vector in 'workingChain'. If the requirement is satisfied, this operation sets the size and
 * the contents of 'workingChain' using the algorithm options set in the constructor. If not NULL,
 * 'workingLogLikelihoodValues' and 'workingLogTargetValues' are set accordingly. This operation
 * currently implements the  The Parallel Adaptive Multilevel Stochastic Simulation Algorithm
 * (PAMSSA).\n
 * It consists of 11 steps:
 <list type=number>
 <item> read options;
 <item> save chain and corresponding target pdf values from previous level
 <item> compute currExponent and sequence of weights for current level and update 'm_logEvidenceFactors'
 <item> create covariance matrix for current level
 <item> create *unified* finite distribution for current level
 <item> decide on using balanced chains or not
 <item> plan for number of linked chains for each node so that all nodes generate the closest possible to the same number of positions
 <item> create vector RV for current level
 <item> scale unified covariance matrix until min <= rejection rate <= max
 <item> sample vector RV of current level
 <item> filter chain if requested.
 </list>
 * At the end, the method computes information gain:
 * \f$ ln( \pi(D|M) ) = E[ln( \pi(D|\theta,M) )] - E[ln( \pi(\theta|D,M) / \pi(\theta|M) )] \f$
 *
 * \see S. H. Cheung and E. E. Prudencio. Parallel adaptive multilevel sampling algorithms for the Bayesian analysis of mathematical models. International Journal for Uncertainty Quantification, 2(3):215237, 2012.  */
  void   generateSequence (BaseVectorSequence<P_V,P_M>& workingChain,
                           ScalarSequence<double>*      workingLogLikelihoodValues,
                           ScalarSequence<double>*      workingLogTargetValues);

  //! Method to calculate the logarithm of the evidence.
  /*! Access to the private member: \c m_logEvidence.*/
  double logEvidence      () const;

  //! Method to calculate the mean of the logarithm of the likelihood.
  /*! Access to the private member: \c m_meanLogLikelihood. */
  double meanLogLikelihood() const;

  //! Calculates the expected information gain value, EIG.
  /*! The expected information gain is defined as the expected log ratio between the posterior and
   prior distribution for the parameters in the statistical model. Recalling the Bayes formula:
   \f[
   \pi (\theta|D, M_j) = \frac{f(D|\theta, M_j) \cdot \pi (\theta | M_j)}{\pi(D, M_j)}
   \f]
   where \f$ \pi (\theta|D, M_j) \f$ is the posterior PDF, \f$ \pi (\theta | M_j) \f$ is the prior,
   \f$ f(D|\theta, M_j)\f$ is the likelihood function and \f$ \pi(D, M_j) \f$ is the evidence for a
   given set of parameters \f$ \theta \f$, data \f$ D \f$  and model class \f$ M_j \f$. Then the EIG
   can be calculated as:
   \f[
   EIG = log \left(\frac{ \pi (\theta|D, M_j)}{\pi (\theta | M_j) } \right)
       = log \left(\frac{f(D|\theta, M_j) }{\pi(D, M_j)} \right)
       = log \left(\frac{f(D|\theta, M_j) }{\pi(D, M_j)} \right)
       = log \left(f(D|\theta, M_j)\right) - log\left(f(D|\theta, M_j) \right)
   \f]

   \see Long, Quan and Scavino, Marco and Tempone, Raul and Wang, Suojin, Fast estimation of expected information gains for Bayesian experimental designs based on Laplace approximations. Computer Methods In Applied Mechanics And Engineering, 259:24-39,2013. DOI = 10.1016/j.cma.2013.02.017.   */
  double eig              () const;
  //@}

  //! @name I/O methods
  //@{
  //! TODO: Prints the sequence.
  /*! \todo: implement me!*/
  void   print           (std::ostream& os) const;
  //@}


  friend std::ostream& operator<<(std::ostream& os,
      const MLSampling<P_V,P_M>& obj) {
    obj.print(os);
    return os;
  }

private:
 //! Writes checkpoint data for the ML method.
 /*!*/
 /*!@param[in]  currExponent, currEta, currChain, currLogLikelihoodValues, currLogTargetValues.*/
  void   checkpointML                  (double                                          currExponent,                       // input
                                        double                                          currEta,                            // input
                                        const SequenceOfVectors<P_V,P_M>&        currChain,                          // input
                                        const ScalarSequence<double>&            currLogLikelihoodValues,            // input
                                        const ScalarSequence<double>&            currLogTargetValues);               // input

 //! Restarts ML algorithm.
 /*! This method reads the control file and determines the number of lines on it; then it reads the stored
  * values, calls MPI_Bcast and processes the read data in all available MPI nodes. */
 /*!@param[out]  currExponent, currEta, currChain, currLogLikelihoodValues, currLogTargetValues.*/
 void   restartML                     (double&                                         currExponent,                       // output
                                        double&                                         currEta,                            // output
                                        SequenceOfVectors<P_V,P_M>&              currChain,                          // output
                                        ScalarSequence<double>&                  currLogLikelihoodValues,            // output
                                        ScalarSequence<double>&                  currLogTargetValues);               // output

 //! Generates the sequence at the level 0.
 /*! @param[in]  currOptions
  @param[out] unifiedRequestedNumSamples, currChain, currLogLikelihoodValues, currLogTargetValues*/
  void   generateSequence_Level0_all   (const MLSamplingLevelOptions&            currOptions,                        // input
                                        unsigned int&                                   unifiedRequestedNumSamples,         // output
                                        SequenceOfVectors<P_V,P_M>&              currChain,                          // output
                                        ScalarSequence<double>&                  currLogLikelihoodValues,            // output
                                        ScalarSequence<double>&                  currLogTargetValues);               // output

  //! Reads options for the ML algorithm (Step 01 from ML algorithm).
 /*! This method is responsible for the Step 01 in the ML algorithm implemented/described in the method MLSampling<P_V,P_M>::generateSequence.*/
 /*! @param[in]  currOptions */
 /*! @param[out] unifiedRequestedNumSamples*/
  void   generateSequence_Step01_inter0(const MLSamplingLevelOptions*            currOptions,                        // input
                                        unsigned int&                                   unifiedRequestedNumSamples);        // output

 //! Saves chain and corresponding target pdf values from previous level (Step 02 from ML algorithm).
 /*! This method is responsible for the Step 02 in the ML algorithm implemented/described in the method MLSampling<P_V,P_M>::generateSequence.*/
 /*! @param[in]  currOptions
  * @param[in,out] currChain,currLogLikelihoodValues,currLogTargetValues
  * @param[out] prevChain,prevLogLikelihoodValues, prevLogTargetValues,indexOfFirstWeight,indexOfLastWeight*/
  void   generateSequence_Step02_inter0(const MLSamplingLevelOptions*            currOptions,                        // input
                                        SequenceOfVectors<P_V,P_M>&              currChain,                          // input/output
                                        ScalarSequence<double>&                  currLogLikelihoodValues,            // input/output
                                        ScalarSequence<double>&                  currLogTargetValues,                // input/output
                                        SequenceOfVectors<P_V,P_M>&              prevChain,                          // output
                                        ScalarSequence<double>&                  prevLogLikelihoodValues,            // output
                                        ScalarSequence<double>&                  prevLogTargetValues,                // output
                                        unsigned int&                                   indexOfFirstWeight,                 // output
                                        unsigned int&                                   indexOfLastWeight);                 // output

 //! Computes \c currExponent and sequence of weights for current level and update 'm_logEvidenceFactors' (Step 03 from ML algorithm).
 /*! This method is responsible for the Step 03 in the ML algorithm implemented/described in the method MLSampling<P_V,P_M>::generateSequence.*/
 /*! @param[in] currOptions, prevLogLikelihoodValues, prevExponent,failedExponent
  * @param[out] currExponent, weightSequence */
  void   generateSequence_Step03_inter0(const MLSamplingLevelOptions*            currOptions,                        // input
                                        const ScalarSequence<double>&            prevLogLikelihoodValues,            // input
                                        double                                          prevExponent,                       // input
                                        double                                          failedExponent,                     // input // gpmsa1
                                        double&                                         currExponent,                       // output
                                        ScalarSequence<double>&                  weightSequence);                    // output

  //! Creates covariance matrix for current level (Step 04 from ML algorithm).
  /*! This method is responsible for the Step 04 in the ML algorithm implemented/described in the method MLSampling<P_V,P_M>::generateSequence.*/
  /*! @param[in] prevChain, weightSequence
  * @param[out] unifiedCovMatrix */
  void   generateSequence_Step04_inter0(const SequenceOfVectors<P_V,P_M>&        prevChain,                          // input
                                        const ScalarSequence<double>&            weightSequence,                     // input
                                        P_M&                                            unifiedCovMatrix);                  // output

  //! Creates \b unified finite distribution for current level (Step 05 from ML algorithm).
  /*! This method is responsible for the Step 05 in the ML algorithm implemented/described in the method MLSampling<P_V,P_M>::generateSequence.*/
  /*! @param[in] unifiedRequestedNumSamples, weightSequence
      @param[out] unifiedIndexCountersAtProc0Only, unifiedWeightStdVectorAtProc0Only */
  void   generateSequence_Step05_inter0(unsigned int                                    unifiedRequestedNumSamples,         // input
                                        const ScalarSequence<double>&            weightSequence,                     // input
                                        std::vector<unsigned int>&                      unifiedIndexCountersAtProc0Only,    // output
                                        std::vector<double>&                            unifiedWeightStdVectorAtProc0Only); // output

  //! Decides on wheter or not to use balanced chains (Step 06 from ML algorithm).
  /*! This method is responsible for the Step 06 in the ML algorithm implemented/described in the method MLSampling<P_V,P_M>::generateSequence.*/
  /*! @param[in] currOptions, indexOfFirstWeight, indexOfLastWeight, unifiedIndexCountersAtProc0Only
      @param[out] useBalancedChains, exchangeStdVec*/
  void   generateSequence_Step06_all   (const MLSamplingLevelOptions*            currOptions,                        // input
                                        unsigned int                                    indexOfFirstWeight,                 // input
                                        unsigned int                                    indexOfLastWeight,                  // input
                                        const std::vector<unsigned int>&                unifiedIndexCountersAtProc0Only,    // input
                                        bool&                                           useBalancedChains,                  // output
                                        std::vector<ExchangeInfoStruct>&              exchangeStdVec);                    // output

  //! Plans for number of linked chains for each node so that all nodes generate the closest possible to the same number of positions (Step 07 from ML algorithm).
  /*! This method is responsible for the Step 07 in the ML algorithm implemented/described in the method MLSampling<P_V,P_M>::generateSequence.*/
  /*! @param[in] useBalancedChains,indexOfFirstWeight,indexOfLastWeight, unifiedIndexCountersAtProc0Only,currOptions,prevChain,prevExponent,currExponent,prevLogLikelihoodValues,prevLogTargetValues
      @param[in,out] exchangeStdVec
      @param[out] unbalancedLinkControl, balancedLinkControl*/
  void   generateSequence_Step07_inter0(bool                                            useBalancedChains,                  // input
                                        unsigned int                                    indexOfFirstWeight,                 // input
                                        unsigned int                                    indexOfLastWeight,                  // input
                                        const std::vector<unsigned int>&                unifiedIndexCountersAtProc0Only,    // input
                                        UnbalancedLinkedChainsPerNodeStruct&          unbalancedLinkControl,              // (possible) output
                                        const MLSamplingLevelOptions*            currOptions,                        // input
                                        const SequenceOfVectors<P_V,P_M>&        prevChain,                          // input
                                        double                                   prevExponent,                       // input
                                        double                                   currExponent,                       // input
                                        const ScalarSequence<double>&            prevLogLikelihoodValues,            // input
                                        const ScalarSequence<double>&            prevLogTargetValues,                // input
                                        std::vector<ExchangeInfoStruct>&              exchangeStdVec,                     // (possible) input/output
                                        BalancedLinkedChainsPerNodeStruct<P_V>&       balancedLinkControl);               // (possible) output

  //! Creates a vector RV for current level (Step 08 from ML algorithm).
  /*! This method is responsible for the Step 08 in the ML algorithm implemented/described in the method MLSampling<P_V,P_M>::generateSequence.*/
  /*! @param[in,out] currPdf
      @param[out] currRv */
  void   generateSequence_Step08_all   (BayesianJointPdf<P_V,P_M>&               currPdf,                            // input/output
                                        GenericVectorRV<P_V,P_M>&                currRv);                            // output

  //! Scales the unified covariance matrix until min <= rejection rate <= max (Step 09 from ML algorithm).
  /*! This method is responsible for the Step 09 in the ML algorithm implemented/described in the method MLSampling<P_V,P_M>::generateSequence.*/
  /*! @param[in] prevChain, indexOfFirstWeight, indexOfLastWeight, unifiedWeightStdVectorAtProc0Only, weightSequence, prevEta,currRv, currOptions,
     @param[in,out]  unifiedCovMatrix
     @param[out] currEta */
  void   generateSequence_Step09_all   (const SequenceOfVectors<P_V,P_M>&        prevChain,                          // input
                                        double                                   prevExponent,                       // input
                                        double                                   currExponent,                       // input
                                        const ScalarSequence<double>&            prevLogLikelihoodValues,            // input
                                        const ScalarSequence<double>&            prevLogTargetValues,                // input
                                        unsigned int                                    indexOfFirstWeight,                 // input
                                        unsigned int                                    indexOfLastWeight,                  // input
                                        const std::vector<double>&                      unifiedWeightStdVectorAtProc0Only,  // input
                                        const ScalarSequence<double>&            weightSequence,                     // input
                                        double                                          prevEta,                            // input
                                        const GenericVectorRV<P_V,P_M>&          currRv,                             // input
                                        MLSamplingLevelOptions*                  currOptions,                        // input (changed temporarily internally)
                                        P_M&                                            unifiedCovMatrix,                   // input/output
                                        double&                                         currEta);                           // output

  //! Samples the vector RV of current level (Step 10 from ML algorithm).
  /*! This method is responsible for the Step 10 in the ML algorithm implemented/described in the method MLSampling<P_V,P_M>::generateSequence.*/
  /*! @param[in] currOptions, unifiedCovMatrix, currRv, useBalancedChains, unbalancedLinkControl, indexOfFirstWeight, prevChain, balancedLinkControl
     @param[out] currChain, cumulativeRawChainRunTime, cumulativeRawChainRejections, currLogLikelihoodValues,currLogTargetValues*/
  void   generateSequence_Step10_all   (MLSamplingLevelOptions&                  currOptions,                        // input (changed temporarily internally)
                                        const P_M&                                      unifiedCovMatrix,                   // input
                                        const GenericVectorRV  <P_V,P_M>&        currRv,                             // input
                                        bool                                            useBalancedChains,                  // input
                                        const UnbalancedLinkedChainsPerNodeStruct&    unbalancedLinkControl,              // input // Round Rock
                                        unsigned int                                    indexOfFirstWeight,                 // input // Round Rock
                                        const SequenceOfVectors<P_V,P_M>&        prevChain,                          // input // Round Rock
                                        double                                   prevExponent,                       // input
                                        double                                   currExponent,                       // input
                                        const ScalarSequence<double>&            prevLogLikelihoodValues,            // input
                                        const ScalarSequence<double>&            prevLogTargetValues,                // input
                                        const BalancedLinkedChainsPerNodeStruct<P_V>& balancedLinkControl,                // input // Round Rock
                                        SequenceOfVectors      <P_V,P_M>&        currChain,                          // output
                                        double&                                         cumulativeRawChainRunTime,          // output
                                        unsigned int&                                   cumulativeRawChainRejections,       // output
                                        ScalarSequence         <double>*         currLogLikelihoodValues,            // output
                                        ScalarSequence         <double>*         currLogTargetValues);               // output

  //! Filters chain (Step 11 from ML algorithm).
  /*! This method is responsible for the Step 11 in the ML algorithm implemented/described in the method MLSampling<P_V,P_M>::generateSequence.*/
  /*! @param[in] currOptions, unifiedRequestedNumSamples, cumulativeRawChainRejections,
      @param[in,out] currChain, currLogLikelihoodValues, currLogTargetValues
      @param[out]  unifiedNumberOfRejections */
  void   generateSequence_Step11_inter0(const MLSamplingLevelOptions*            currOptions,                        // input
                                        unsigned int                                    unifiedRequestedNumSamples,         // input
                                        unsigned int                                    cumulativeRawChainRejections,       // input
                                        SequenceOfVectors<P_V,P_M>&              currChain,                          // input/output
                                        ScalarSequence<double>&                  currLogLikelihoodValues,            // input/output
                                        ScalarSequence<double>&                  currLogTargetValues,                // input/output
                                        unsigned int&                                   unifiedNumberOfRejections);         // output

  /*! @param[in] unifiedRequestedNumSamples, unifiedWeightStdVectorAtProc0Only
      @param[out] unifiedIndexCountersAtProc0Only*/
  void   sampleIndexes_proc0           (unsigned int                                    unifiedRequestedNumSamples,         // input
                                        const std::vector<double>&                      unifiedWeightStdVectorAtProc0Only,  // input
                                        std::vector<unsigned int>&                      unifiedIndexCountersAtProc0Only);   // output

  /*! @param[in] currOptions, indexOfFirstWeight, indexOfLastWeight, unifiedIndexCountersAtProc0Only
      @param[out] exchangeStdVec*/
  bool   decideOnBalancedChains_all    (const MLSamplingLevelOptions*            currOptions,                        // input
                                        unsigned int                                    indexOfFirstWeight,                 // input
                                        unsigned int                                    indexOfLastWeight,                  // input
                                        const std::vector<unsigned int>&                unifiedIndexCountersAtProc0Only,    // input
                                        std::vector<ExchangeInfoStruct>&              exchangeStdVec);                    // output

   /*! @param[in] currOptions, prevChain, exchangeStdVec
    *  @param[out] exchangeStdVec, balancedLinkControl*/
   void   prepareBalLinkedChains_inter0 (const MLSamplingLevelOptions*            currOptions,                        // input
                                         const SequenceOfVectors<P_V,P_M>&        prevChain,                          // input
                                         double                                   prevExponent,                       // input
                                         double                                   currExponent,                       // input
                                         const ScalarSequence<double>&            prevLogLikelihoodValues,            // input
                                         const ScalarSequence<double>&            prevLogTargetValues,                // input
                                         std::vector<ExchangeInfoStruct>&              exchangeStdVec,                     // input/output
                                         BalancedLinkedChainsPerNodeStruct<P_V>&       balancedLinkControl);               // output

  /*! @param[in] indexOfFirstWeight, indexOfLastWeight, unifiedIndexCountersAtProc0Only
   *  @param[out] unbalancedLinkControl*/
   void   prepareUnbLinkedChains_inter0 (unsigned int                                    indexOfFirstWeight,                 // input
                                        unsigned int                                    indexOfLastWeight,                  // input
                                        const std::vector<unsigned int>&                unifiedIndexCountersAtProc0Only,    // input
                                        UnbalancedLinkedChainsPerNodeStruct&          unbalancedLinkControl);             // output

   /*! @param[in] inputOptions, unifiedCovMatrix, rv, balancedLinkControl,
    *  @param[out] workingChain, cumulativeRunTime, cumulativeRejections, currLogLikelihoodValues, currLogTargetValues*/
   void   generateBalLinkedChains_all   (MLSamplingLevelOptions&                  inputOptions,                       // input, only m_rawChainSize changes
                                        const P_M&                                      unifiedCovMatrix,                   // input
                                        const GenericVectorRV  <P_V,P_M>&        rv,                                 // input
                                        const BalancedLinkedChainsPerNodeStruct<P_V>& balancedLinkControl,                // input // Round Rock
                                        SequenceOfVectors      <P_V,P_M>&        workingChain,                       // output
                                        double&                                         cumulativeRunTime,                  // output
                                        unsigned int&                                   cumulativeRejections,               // output
                                        ScalarSequence         <double>*         currLogLikelihoodValues,            // output
                                        ScalarSequence         <double>*         currLogTargetValues);               // output

    /*! @param[in] inputOptions, unifiedCovMatrix, rv, unbalancedLinkControl, indexOfFirstWeight, prevChain
     *  @param[out] workingChain, cumulativeRunTime, cumulativeRejections, currLogLikelihoodValues, currLogTargetValues*/
    void   generateUnbLinkedChains_all   (MLSamplingLevelOptions&                  inputOptions,                       // input, only m_rawChainSize changes
                                          const P_M&                                      unifiedCovMatrix,                   // input
                                          const GenericVectorRV  <P_V,P_M>&        rv,                                 // input
                                          const UnbalancedLinkedChainsPerNodeStruct&    unbalancedLinkControl,              // input // Round Rock
                                          unsigned int                                    indexOfFirstWeight,                 // input // Round Rock
                                          const SequenceOfVectors<P_V,P_M>&        prevChain,                          // input // Round Rock
                                          double                                   prevExponent,                       // input
                                          double                                   currExponent,                       // input
                                          const ScalarSequence<double>&            prevLogLikelihoodValues,            // input
                                          const ScalarSequence<double>&            prevLogTargetValues,                // input
                                          SequenceOfVectors      <P_V,P_M>&        workingChain,                       // output
                                          double&                                         cumulativeRunTime,                  // output
                                          unsigned int&                                   cumulativeRejections,               // output
                                          ScalarSequence         <double>*         currLogLikelihoodValues,            // output
                                          ScalarSequence         <double>*         currLogTargetValues);               // output

#ifdef QUESO_HAS_GLPK
  /*! @param[in] exchangeStdVec
   *  @param[out] exchangeStdVec*/
  void   solveBIP_proc0                (std::vector<ExchangeInfoStruct>&              exchangeStdVec);                    // input/output
#endif

  /*! @param[in] currOptions, exchangeStdVec
   *  @param[out] exchangeStdVec*/
  void   justBalance_proc0             (const MLSamplingLevelOptions*            currOptions,                        // input
                                        std::vector<ExchangeInfoStruct>&              exchangeStdVec);                    // input/output

  /*! @param[in] prevChain, exchangeStdVec, finalNumChainsPerNode, finalNumPositionsPerNode
   *  @param[out] balancedLinkControl*/
  void   mpiExchangePositions_inter0   (const SequenceOfVectors<P_V,P_M>&        prevChain,                          // input
                                        double                                   prevExponent,                       // input
                                        double                                   currExponent,                       // input
                                        const ScalarSequence<double>&            prevLogLikelihoodValues,            // input
                                        const ScalarSequence<double>&            prevLogTargetValues,                // input
                                        const std::vector<ExchangeInfoStruct>&        exchangeStdVec,                     // input
                                        const std::vector<unsigned int>&                finalNumChainsPerNode,              // input
                                        const std::vector<unsigned int>&                finalNumPositionsPerNode,           // input
                                        BalancedLinkedChainsPerNodeStruct<P_V>&       balancedLinkControl);               // output

  // Private variables
  //! Queso enviroment.
  const BaseEnvironment&             m_env;

  //!  Prior RV.
  const BaseVectorRV      <P_V,P_M>& m_priorRv;

  //! Likelihood function.
  const BaseScalarFunction<P_V,P_M>& m_likelihoodFunction;

  //! Vector space.
  const VectorSpace       <P_V,P_M>& m_vectorSpace;

   //! Domain of the target PDF: intersection of the domains of the prior PDf and likelihood function.
  VectorSet         <P_V,P_M>* m_targetDomain;

  unsigned int                        m_numDisabledParameters; // gpmsa2

  std::vector<bool>                   m_parameterEnabledStatus; // gpmsa2

   //! Options for the ML algorithm.
   MLSamplingOptions            m_options;

   //! Current level.
   unsigned int                        m_currLevel;          // restart

   //! Curret step.
   unsigned int                        m_currStep;

   //! Exponent for debugging.
   double                              m_debugExponent;
	std::vector<double>                 m_logEvidenceFactors; // restart
        double                              m_logEvidence;
        double                              m_meanLogLikelihood;
        double                              m_eig;
};

}  // End namespace QUESO

#endif // UQ_MULTI_LEVEL_SAMPLING_H
