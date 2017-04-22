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

#ifndef UQ_MH_SG_H
#define UQ_MH_SG_H

#include <queso/MetropolisHastingsSGOptions.h>
#include <queso/TKGroup.h>
#include <queso/VectorRV.h>
#include <queso/VectorSpace.h>
#include <queso/MarkovChainPositionData.h>
#include <queso/ScalarFunctionSynchronizer.h>
#include <queso/SequenceOfVectors.h>
#include <queso/ArrayOfSequences.h>
#include <sys/time.h>
#include <fstream>
#include <queso/SharedPtr.h>

namespace QUESO {

class GslVector;
class GslMatrix;

template <class P_V, class P_M>
class Algorithm;

//--------------------------------------------------
// MHRawChainInfoStruct --------------------------
//--------------------------------------------------
 /*! \file MetropolisHastingsSG.h
 * \brief A templated class that represents a Metropolis-Hastings generator of samples and a struct which stores its info.
 *
 * \struct  MHRawChainInfoStruct
 * \brief A struct that represents a Metropolis-Hastings sample.
 *
 * Some of the information about the  Metropolis-Hastings sample generator includes the allowed number
 * of delayed rejections, number of rejections, number of positions in or out of target support, and
 * so on. This struct is responsible for the storage of such info. */

struct MHRawChainInfoStruct
{
 //! @name Constructor/Destructor methods
 //@{
 //! Constructor.
  MHRawChainInfoStruct();

  //! Copy constructor.
  MHRawChainInfoStruct(const MHRawChainInfoStruct& rhs);

  //! Destructor
  ~MHRawChainInfoStruct();
  //@}

  //! @name Set methods
  //@{
  //! Assignment operator.
  MHRawChainInfoStruct& operator= (const MHRawChainInfoStruct& rhs);

  //! Addition assignment operator.
  MHRawChainInfoStruct& operator+=(const MHRawChainInfoStruct& rhs);
  //@}

   //! @name Misc methods
  //@{
  //! Copies Metropolis-Hastings chain info from \c src to \c this.
  void copy  (const MHRawChainInfoStruct& src);

  //! Resets Metropolis-Hastings chain info.
  void reset ();

  //! Calculates the MPI sum of \c this.
  void mpiSum(const MpiComm& comm, MHRawChainInfoStruct& sumInfo);
  //@}

  double       runTime;
  double       candidateRunTime;
  double       targetRunTime;
  double       mhAlphaRunTime;
  double       drAlphaRunTime;
  double       drRunTime;
  double       amRunTime;

  unsigned int numTargetCalls;
  unsigned int numDRs;
  unsigned int numOutOfTargetSupport;
  unsigned int numOutOfTargetSupportInDR;
  unsigned int numRejections;

};

//--------------------------------------------------
// MetropolisHastingsSG----------------------
//--------------------------------------------------

/*!\class MetropolisHastingsSG
 * \brief A templated class that represents a Metropolis-Hastings generator of samples.
 *
 * This class implements a Metropolis-Hastings generator of samples. 'SG' stands for 'Sequence Generator'.
 * Options reading is handled by class 'MetropolisHastingsOptions'. If options request data to be
 * written in the output file (MATLAB .m format only, for now), the user can check which MATLAB variables
 * are defined and set by running 'grep zeros <OUTPUT FILE NAME>' after the solution procedures ends.
 * The names of the variables are self explanatory. */

template <class P_V = GslVector, class P_M = GslMatrix>
class MetropolisHastingsSG
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor.
  /*!
   * This method reads the options from the options input file. It calls
   * \c commonConstructor().
   *
   * Requirements:
   *   -# the image set of the vector random variable 'sourceRv' should belong
   *      to a vector space of dimension equal to the size of the vector
   *      'initialPosition' and
   *   -# if 'inputProposalCovMatrix' is not NULL, is should be square and its
   *      size should be equal to the size of 'initialPosition'.
   *
   * If \c alternativeOptionsValues is NULL and an input file is specified, the
   * constructor reads input options that begin with the string
   * '\<prefix\>mh_'.  For instance, if 'prefix' is 'pROblem_775_ip_', then
   * the constructor will read all options that begin with 'pROblem_775_ip_mh_'.
   *
   * If \c alternativeOptionsValues is not NULL, the input file is ignored and
   * construction copies the object pointed to by \c alternativeOptionsValues
   * to and stores the copy internally.  Users may delete the object poined to
   * by \c alternativeOptionsValues.  Users cannot change the options object
   * after MetropolisHastingsSG has been constructed.
   *
   * Options reading is handled by class 'MetropolisHastingsOptions'.
   */
  MetropolisHastingsSG(const char*                         prefix,
		       const MhOptionsValues*       alternativeOptionsValues, // dakota
		       const BaseVectorRV<P_V,P_M>& sourceRv,
		       const P_V&                          initialPosition,
		       const P_M*                          inputProposalCovMatrix);

  //! Constructor.
  MetropolisHastingsSG(const char*                         prefix,
		       const MhOptionsValues*       alternativeOptionsValues, // dakota
		       const BaseVectorRV<P_V,P_M>& sourceRv,
		       const P_V&                          initialPosition,
		       double                              initialLogPrior,
		       double                              initialLogLikelihood,
		       const P_M*                          inputProposalCovMatrix);

  //! Constructor.
  MetropolisHastingsSG(const MLSamplingLevelOptions& mlOptions,
		       const BaseVectorRV<P_V,P_M>&  sourceRv,
		       const P_V&                           initialPosition,
		       const P_M*                           inputProposalCovMatrix);

  //! Constructor.
  MetropolisHastingsSG(const MLSamplingLevelOptions& mlOptions,
		       const BaseVectorRV<P_V,P_M>&  sourceRv,
		       const P_V&                           initialPosition,
		       double                               initialLogPrior,
		       double                               initialLogLikelihood,
		       const P_M*                           inputProposalCovMatrix);

  //! Destructor
  ~MetropolisHastingsSG();
  //@}

  //! @name Statistical methods
  //@{
 //! Method to generate the chain.
 /*! Requirement: the vector space 'm_vectorSpace' should have dimension equal to the size of a
  * vector in 'workingChain'. If the requirement is satisfied, this operation sets the size and
  * the contents of 'workingChain' using the algorithm options set in the constructor. If not NULL,
  * 'workingLogLikelihoodValues' and 'workingLogTargetValues' are set accordingly. This operation
  * currently implements the DRAM algorithm (Heikki Haario, Marko Laine, Antonietta Mira and
  * Eero Saksman, "DRAM: Efficient Adaptive MCMC", Statistics and Computing (2006), 16:339-354),
  * as a translation of the core routine at the MCMC toolbox for MATLAB, available at
  *  \htmlonly helios.fmi.fi/~lainema/mcmc/‎ \endhtmlonly (accessed in July 3rd, 2013). Indeed,
  * the example available in examples/statisticalInverseProblem is closely related to the
  * 'normal example' in the toolbox. Over time, though:
  <list type=number>
  <item> the whole set of QUESO classes took shape, focusing not only on Markov Chains, but on
  statistical forward problems and model validation as well;
  <item> the interfaces to this Metropolis-Hastings class changed;
  <item> QUESO has parallel capabilities;
  <item> TK (transition kernel) class has been added in order to have both DRAM with Stochastic
  Newton capabilities.
  </list>
  */
  void         generateSequence   (BaseVectorSequence<P_V,P_M>& workingChain,
                                   ScalarSequence<double>*      workingLogLikelihoodValues,
                                   ScalarSequence<double>*      workingLogTargetValues);

  //! Gets information from the raw chain.
  void         getRawChainInfo    (MHRawChainInfoStruct& info) const;

   //@}

  //! Returns the underlying transition kernel for this sequence generator
  const BaseTKGroup<P_V, P_M> & transitionKernel() const;

  //! @name I/O methods
  //@{
  //! TODO: Prints the sequence.
  /*! \todo: implement me!*/
  void   print                    (std::ostream& os) const;
  friend std::ostream& operator<<(std::ostream& os,
      const MetropolisHastingsSG<P_V,P_M>& obj)
  {
    obj.print(os);

    return os;
  }
  //@}

private:
  //! Reads the options values from the options input file.
  /*!  This method \b actually reads the options input file, such as the value for the Delayed Rejection
   * scales, the presence of Hessian covariance matrices and reads the user-provided initial proposal
   * covariance matrix (alternative to local Hessians).*/
  void   commonConstructor        ();

  //!  This method generates the chain.
  /*! Given the size of the chain, the values of the first position in the chain. It checks if the position
   * is out of target pdf support, otherwise it calculates the value of both its likelihood and prior and
   * adds the position to the chain. For the next positions, once they are generated, some tests are
   * performed (such as unicity and the value of alpha) and the steps for the first position are repeated,
   * including the optional Delayed Rejection and the adaptive Metropolis (adaptation of covariance matrix)
   * steps.*/
  void   generateFullChain        (const P_V&                          valuesOf1stPosition,
                                   unsigned int                        chainSize,
                                   BaseVectorSequence<P_V,P_M>& workingChain,
                                   ScalarSequence<double>*      workingLogLikelihoodValues,
                                   ScalarSequence<double>*      workingLogTargetValues);

  //! Adaptive Metropolis method that deals with adapting the proposal covariance matrix
  void adapt(unsigned int positionId,
      BaseVectorSequence<P_V, P_M> & workingChain);

  //! Does delayed rejection
  /*!
   * When faced with an imminent rejection, this method computes a series of
   * new candidates in the same proposal direction but with smaller proposal
   * step sizes.  For each of these new candidates, we check if it will be
   * accepted, if not we repeat the process until all the candidates have been
   * tested.
   *
   * If there is a candidate that will be accepted, this method will return \c
   * true, otherwise it returns \c false.
   *
   * \c currentCandidateData is updated whenever a new proposal is generated
   * throughout the delayed rejection procedure.
   *
   * \c delayedRejection promises not to change \c currentPositionData, because
   * changing the current position of the Markov chain would be grossly
   * inappropriate.
   */
  bool delayedRejection(unsigned int positionId,
      const MarkovChainPositionData<P_V> & currentPositionData,
      MarkovChainPositionData<P_V> & currentCandidateData);

  //! This method reads the chain contents.
  void   readFullChain            (const std::string&                  inputFileName,
                                   const std::string&                  inputFileType,
                                   unsigned int                        chainSize,
                                   BaseVectorSequence<P_V,P_M>& workingChain);

  //! This method updates the adapted covariance matrix
  /*! This function is called is the option to used adaptive Metropolis was chosen by the user
   * (via options input file). It performs an adaptation of covariance matrix. */
  void   updateAdaptedCovMatrix   (const BaseVectorSequence<P_V,P_M>&  subChain,
                                   unsigned int                               idOfFirstPositionInSubChain,
                                   double&                                    lastChainSize,
                                   P_V&                                       lastMean,
                                   P_M&                                       lastAdaptedCovMatrix);

  //! Calculates acceptance ratio.
  /*! The acceptance ratio is used to decide whether to accept or reject a candidate. */
  double alpha                    (const std::vector<MarkovChainPositionData<P_V>*>& inputPositions,
                                   const std::vector<unsigned int                        >& inputTKStageIds);

  //! Decides whether or not to accept alpha.
  /*! If either alpha is negative or greater than one, its value will not be accepted.*/
  bool   acceptAlpha              (double                                     alpha);

  //! Writes information about the Markov chain in a file.
  /*! It writes down the alpha quotients, the number of rejected positions, number of positions out of
   * target support, the name of the components and the chain runtime.*/
  int    writeInfo                (const BaseVectorSequence<P_V,P_M>&  workingChain,
                                   std::ofstream&                             ofsvar) const;

  const BaseEnvironment & m_env;
  const VectorSpace <P_V,P_M> & m_vectorSpace;
  const BaseJointPdf<P_V,P_M> & m_targetPdf;
  P_V m_initialPosition;
  P_M m_initialProposalCovMatrix;
  bool m_nullInputProposalCovMatrix;
  unsigned int m_numDisabledParameters; // gpmsa2
  std::vector<bool> m_parameterEnabledStatus; // gpmsa2
  typename ScopedPtr<const ScalarFunctionSynchronizer<P_V,P_M> >::Type m_targetPdfSynchronizer;

  typename SharedPtr<BaseTKGroup<P_V,P_M> >::Type m_tk;
  typename SharedPtr<Algorithm<P_V, P_M> >::Type m_algorithm;
  unsigned int m_positionIdForDebugging;
  unsigned int m_stageIdForDebugging;
  std::vector<unsigned int> m_idsOfUniquePositions;
  std::vector<double> m_logTargets;
  std::vector<double> m_alphaQuotients;
  double m_lastChainSize;
  typename ScopedPtr<P_V>::Type m_lastMean;
  typename ScopedPtr<P_M>::Type m_lastAdaptedCovMatrix;
  unsigned int m_numPositionsNotSubWritten;

  MHRawChainInfoStruct m_rawChainInfo;

  ScopedPtr<const MhOptionsValues>::Type m_optionsObj;
  ScopedPtr<MetropolisHastingsSGOptions>::Type m_oldOptions;

	bool m_computeInitialPriorAndLikelihoodValues;
	double m_initialLogPriorValue;
	double m_initialLogLikelihoodValue;

  void transformInitialCovMatrixToGaussianSpace(const BoxSubset<P_V, P_M> &
      boxSubset);

  bool m_userDidNotProvideOptions;

  unsigned int m_latestDirtyCovMatrixIteration;
};

}  // End namespace QUESO

#endif // UQ_MH_SG_H
