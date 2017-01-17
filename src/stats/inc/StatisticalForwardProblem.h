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

#ifndef UQ_SFP_H
#define UQ_SFP_H

#include <queso/StatisticalForwardProblemOptions.h>
#include <queso/VectorFunction.h>
#include <queso/MonteCarloSG.h>
#include <queso/VectorRV.h>
#include <queso/GenericVectorRV.h>
#include <queso/SequenceOfVectors.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*! \file StatisticalForwardProblem.h
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

template <class P_V = GslVector, class P_M = GslMatrix, class Q_V = GslVector, class Q_M = GslMatrix>
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
  friend std::ostream& operator<<(std::ostream& os,
      const StatisticalForwardProblem<P_V,P_M,Q_V,Q_M>& obj) {
    obj.print(os);
    return os;
  }
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

        const SfpOptionsValues * m_optionsObj;

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

  bool m_userDidNotProvideOptions;
};

}  // End namespace QUESO

#endif // UQ_SFP_H
