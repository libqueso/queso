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

#ifndef UQ_SIP_H
#define UQ_SIP_H

#include <queso/StatisticalInverseProblemOptions.h>
#include <queso/MetropolisHastingsSG.h>
#include <queso/MLSampling.h>
#include <queso/InstantiateIntersection.h>
#include <queso/VectorRealizer.h>
#include <queso/SequentialVectorRealizer.h>
#include <queso/VectorRV.h>
#include <queso/ScalarFunction.h>
#include <queso/GPMSA.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*! \file StatisticalInverseProblem.h
    \brief Class to solve a Statistical Inverse Problem
*/

/*! \class StatisticalInverseProblem
 *  \brief This templated class represents a Statistical Inverse Problem.
 *
 * This class is templated on the type 'P_V' of vector and type 'P_M' of matrix, where
 * 'P_' stands for 'parameter'. Some templated classes might also use 'Q_V' and 'Q_M'
 * when referring to a vector type and a matrix type, where 'Q_' stands for 'quantities
 * of interest'. Conceptually, a statistical inverse problem has two input entities and
 * one output entity.\n
 *
 * The input entities of a statistical inverse problem are:
<list type=number>
<item> the prior RV, an instance of class 'BaseVectorRV<P_V,P_M>', and
<item> the likelihood function, an instance of class 'BaseScalarFunction<P_V,P_M>'.
</list>
 * Let \f$ \pi()\f$ denote the mathematical likelihood function and \f$ x \f$ denote a
 * vector of parameters. The likelihood function object stores the routine that computes
 * \f$ \pi(x) \f$ and whatever data necessary by such routine. The routine in the
 * likelihood function object can compute either the actual value \f$ \pi(x) \f$ or the
 * value \f$ \ln[\pi(x)] \f$. See files
 * 'basic/inc/ScalarFunction.h' and 'stats/inc/JointPdf.h' for more details.\n
 *
 * The output entity of a statistical inverse problem is the posterior RV, another instance
 * of class 'BaseVectorRV<P_V,P_M>', which stores the solution according to the Bayesian
 * approach. Upon return from a solution operation, the posterior RV is available through the
 * operation 'postRv()'. Such posterior RV is able to provide:
<list type=number>
<item> a joint pdf (up to a multiplicative constant) through the operation 'postRv().pdf()',
       which returns an instance of the class 'BaseJointPdf<P_V,P_M>', and
<item> a vector realizer through the operation 'postRv().realizer()', which returns an
       instance of the class 'BaseVectorRealizer<P_V,P_M>'.
</list>*/

/* OLD STUFF: The constructor of the 'scalar function' asks for the user to specify which
 * value the routine is actually computing, so that the 'scalar function' class can properly
 * implements both class operations 'actualValue()' and 'minus2LnValue()'*/

template <class P_V = GslVector, class P_M = GslMatrix>
class StatisticalInverseProblem
{
public:
 //! @name Constructor/Destructor methods
 //@{
 //! Constructor.
 /*! Requirements: 1) the image set of the vector random variable 'priorRv', 2) the domain set of
 * the likelihood function 'likelihoodFunction' and 3) the image set of the vector random variable
 * 'postRv' should belong to vector spaces of equal dimensions. If the requirements are satisfied,
 * the constructor then reads input options that begin with the string '\<prefix\>ip_'. If no
 * options input file is provided, the construction assigns \c alternativeOptionsValues to the
 * options of the SIP.*/
  StatisticalInverseProblem(const char*                               prefix,
                                   const SipOptionsValues*            alternativeOptionsValues, // dakota
                                   const BaseVectorRV      <P_V,P_M>& priorRv,
                                   const BaseScalarFunction<P_V,P_M>& likelihoodFunction,
                                         GenericVectorRV   <P_V,P_M>& postRv);

  //! Constructor for statistical inverse problems to be sovled using GPMSA
  /*!
   * Requirements:
   *   1) the factory for the GPMSA object
   *   2) the image set of the vector random variable \c postRv (obtainable via
   *      GaussianProcessFactory::prior method) should be equal to the full
   *      prior image set (including all the hyperparameters)
   *
   * If the requirements are satisfied, the constructor then reads input
   * options that begin with the string '\<prefix\>ip_'.  If no options input
   * file is provided, the construction assigns \c alternativeOptionsValues to
   * the options of the statistical inverse problem.
   */
  StatisticalInverseProblem(const char * prefix,
                            const SipOptionsValues * alternativeOptionsValues,
                            const GPMSAFactory<P_V,P_M> & gpmsaFactory,
                            GenericVectorRV<P_V,P_M> & postRv);

  //! Destructor
  ~StatisticalInverseProblem();
  //@}

  //! @name Statistical methods
  //@{
  //! Whether or not compute the solution.
  bool                             computeSolutionFlag             () const;

  //! Solves the problem via Bayes formula and a Metropolis-Hastings algorithm.
  /*!
   * Requirements:
   * <list type=number>
   * <item> 'initialValues' should have the same number of components as member
   *        variable 'm_priorRv'
   * <item> if 'initialProposalCovMatrix' is not NULL, it should be square and
   *        its size should be equal to the size of 'initialValues'.
   * </list>
   *
   * If the requirements are satisfied, this methods checks the member flag
   * 'm_computeSolution' (one of the options read from the input file during
   * construction). If the flag is 'false', the operation returns immediately,
   * computing nothing; otherwise, the operation sets the member variable
   * 'm_postRv' accordingly. The operation:
   * <list type=number>
   * <item> sets the pdf of 'm_postRv' equal to an instance of
   *        BayesianJointPdf<P_V,P_M>,
   * <item> instantiates SequenceOfVectors<P_V,P_M> (the chain),
   * <item> instantiates MetropolisHastingsSG<P_V,P_M> (the Metropolis-Hastings
   *        algorithm),
   * <item> populates the chain with the Metropolis-Hastings algorithm, and
   * <item> sets the realizer of 'm_postRv' with the contents of the chain.
   * </list>
   *
   * If \c seedWithMapEstimator() is called before this method, then
   * \c initialValues is used as the initial condition for an optimization
   * procedure, the result of which will be used as the seed for the chain.
   */
  void solveWithBayesMetropolisHastings(const MhOptionsValues* alternativeOptionsValues, // dakota
					const P_V&                    initialValues,
					const P_M*                    initialProposalCovMatrix);

  //! Seeds the chain with the result of a deterministic optimisation
  /*!
   * This only works for Metropolis-Hastings right now.  Multi-level is not
   * currently supported.
   */
  void seedWithMAPEstimator();

  //! Solves with Bayes Multi-Level (ML) sampling.
  void                             solveWithBayesMLSampling        ();

  //! Return the underlying MetropolisHastingSG object
  const MetropolisHastingsSG<P_V, P_M> & sequenceGenerator() const;

  //! Returns the Prior RV; access to private attribute m_priorRv.
  const BaseVectorRV   <P_V,P_M>& priorRv                   () const;

  //! Returns the Posterior RV; access to private attribute m_postrRv.
  /*! The Posterior RV contains the solution of the Bayes problem.*/
  const GenericVectorRV<P_V,P_M>& postRv                    () const;

  //! Returns the MCMC chain; access to private attribute m_chain.
  /*! Only valid after solve has been called.*/
  const BaseVectorSequence<P_V,P_M>& chain() const;

  //! Returns log likelihood values; access to private attribute m_logLikelihoodValues.
  /*! Only valid for MH and only after solve has been called.*/
  const ScalarSequence<double>& logLikelihoodValues() const;

  //! Returns log target values; access to private attribute m_logTargetValues.
  /*! Only valid for MH and only after solve has been called.*/
  const ScalarSequence<double>& logTargetValues() const;

  //! Returns the logarithm value of the evidence. Related to ML.
  double                           logEvidence                     () const;

   //! Returns the mean of the logarithm value of the likelihood. Related to ML.
  double                           meanLogLikelihood               () const;

  //\TODO Related to ML.
  double                           eig                             () const;
  //@}

  //! @name I/O methods
  //@{
  //! TODO: Prints the sequence.
  /*! \todo: implement me!*/
  void print(std::ostream& os) const;
  friend std::ostream& operator<<(std::ostream& os,
      const StatisticalInverseProblem<P_V,P_M>& obj) {
    obj.print(os);
    return os;
  }
  //@}
private:
  const BaseEnvironment&                 m_env;

  const BaseVectorRV        <P_V,P_M>&   m_priorRv;
  const BaseScalarFunction  <P_V,P_M>&   m_likelihoodFunction;
        GenericVectorRV     <P_V,P_M>&   m_postRv;

        VectorSet           <P_V,P_M>*   m_solutionDomain;
        BaseJointPdf        <P_V,P_M>*   m_solutionPdf;
        BaseVectorMdf       <P_V,P_M>*   m_subSolutionMdf;
        BaseVectorCdf       <P_V,P_M>*   m_subSolutionCdf;
        BaseVectorRealizer  <P_V,P_M>*   m_solutionRealizer;

        MetropolisHastingsSG<P_V,P_M>*   m_mhSeqGenerator;
        MLSampling          <P_V,P_M>*   m_mlSampler;
        BaseVectorSequence  <P_V,P_M>*   m_chain;
        ScalarSequence      <double>*    m_logLikelihoodValues;
        ScalarSequence      <double>*    m_logTargetValues;

        const SipOptionsValues * m_optionsObj;

        bool                              m_seedWithMAPEstimator;

#ifdef UQ_ALSO_COMPUTE_MDFS_WITHOUT_KDE
        ArrayOfOneDGrids    <P_V,P_M>*   m_subMdfGrids;
        ArrayOfOneDTables   <P_V,P_M>*   m_subMdfValues;
#endif

  bool m_userDidNotProvideOptions;
};

}  // End namespace QUESO

#endif // UQ_SIP_H
