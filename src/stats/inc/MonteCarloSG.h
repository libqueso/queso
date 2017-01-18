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

#ifndef UQ_MOC_SG_H
#define UQ_MOC_SG_H

#include <queso/VectorRV.h>
#include <queso/VectorFunction.h>
#include <queso/VectorFunctionSynchronizer.h>
#include <queso/MonteCarloSGOptions.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*!
 * \file MonteCarloSG.h
 * \brief A templated class that implements a Monte Carlo generator of samples.
 *
 * \class MonteCarloSG
 * \brief A templated class that implements a Monte Carlo generator of samples.
 *
 * This class implements a Monte Carlo generator of samples. 'SG' stands for 'Sequence Generator'.
 * Options reading is handled by class 'MonteCarloOptions'. If options request data to be
 * written in the output file (MATLAB .m format only, for now), the user can check which MATLAB
 * variables are defined and set by running 'grep zeros <OUTPUT FILE NAME>' after the solution
 * procedures ends. The names of the variables are self explanatory. */

template <class P_V = GslVector, class P_M = GslMatrix, class Q_V = GslVector, class Q_M = GslMatrix>
class MonteCarloSG
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor.
 /*! Requirements: 1) the image set of the vector random variable 'paramRv' and 2) the domain set of the
 * QoI function 'qoiFunction' should belong to vector spaces of equal dimensions. If the requirements
 * are satisfied, the constructor then reads input options that begin with the string '\<prefix\>_mc_'.
 * For instance, if 'prefix' is 'pROblem_775_fp_', then the constructor will read all options that begin
 * with 'pROblem_775_fp_mc_'. Options reading is handled by class 'MonteCarloOptions'.*/
  MonteCarloSG(const char*                                       prefix,
                      const McOptionsValues*                     alternativeOptionsValues, // dakota
                      const BaseVectorRV      <P_V,P_M>&         paramRv,
                      const BaseVectorFunction<P_V,P_M,Q_V,Q_M>& qoiFunction);

  //! Destructor.
  ~MonteCarloSG();
  //@}

  //! @name Statistical methods
  //@{
  //! Generates the QoI (output) sequence, it calls internGenerateSequence().
  /*! This method checks for a parallel environment (and uses it if available) and calls the private
   * member function internGenerateSequence() to generate a sequence of values of the Quantity of
   * interest (QoI).*/
  void generateSequence(BaseVectorSequence<P_V,P_M>& workingPSeq,
                        BaseVectorSequence<Q_V,Q_M>& workingQSeq);
  //@}

  //! @name I/O methods
  //@{
  //! Prints the sequence.
  void print                      (std::ostream& os) const;
  friend std::ostream& operator<<(std::ostream& os,
      const MonteCarloSG<P_V,P_M,Q_V,Q_M>& obj) {
    obj.print(os);
    return os;
  }
  //@}
private:
  //! Generates the QoI (output) sequence; it calls actualGenerateSequence().
  /*!  This method calls the private member actualGenerateSequence() to \b actually generate a
   * sequence of values of the Quantity of interest (QoI). It writes both the QoI and the parameter
   * sequences to output files. If a parallel environment is available, the user has the choice of
   * (via the input options file) which sub-environment will write data to files. The unified
   * sequence of data (QoI and/or Parameter sequences) of all the sub-environments (or of the single
   * environment, in case of a serial run) is also written to files.*/
  void internGenerateSequence(const BaseVectorRV      <P_V,P_M>& paramRv,
                                    BaseVectorSequence<P_V,P_M>& workingPSeq,
                                    BaseVectorSequence<Q_V,Q_M>& workingQSeq);


  //!  This method \b actually generates the QoI sequence.
  /*! Requirements: 1) the vector space containing the domain set of the QoI function 'm_qoiFunction'
   * should have dimension equal to the size of a vector in 'workingPSeq' and 2) the vector space
   * containing the image set of the qoi function 'm_qoiFunction' should have dimension equal to the
   * size of a vector in 'workingQSeq'. If the requirements are satisfied, this operation sets the
   * size and the contents of 'workingPSeq' and 'workingQSeq' using the algorithm options set in the
   * constructor. */
  void actualGenerateSequence(const BaseVectorRV      <P_V,P_M>& paramRv,
                                    BaseVectorSequence<P_V,P_M>& workingPSeq,
                                    BaseVectorSequence<Q_V,Q_M>& workingQSeq,
                                    unsigned int                        seqSize);
  //! Reads the sequence.
  void actualReadSequence    (const BaseVectorRV      <P_V,P_M>& paramRv,
                              const std::string&                        dataInputFileName,
                              const std::string&                        dataInputFileType,
                                    BaseVectorSequence<P_V,P_M>& workingPSeq,
                                    BaseVectorSequence<Q_V,Q_M>& workingQSeq,
                                    unsigned int                        seqSize);

  const BaseEnvironment&                             m_env;
  const BaseVectorRV              <P_V,P_M>&         m_paramRv;
  const BaseVectorFunction        <P_V,P_M,Q_V,Q_M>& m_qoiFunction;
  const VectorSpace               <P_V,P_M>&         m_paramSpace;
  const VectorSpace               <Q_V,Q_M>&         m_qoiSpace;
  const VectorFunctionSynchronizer<P_V,P_M,Q_V,Q_M>* m_qoiFunctionSynchronizer;
  unsigned int                                              m_numPsNotSubWritten;
  unsigned int                                              m_numQsNotSubWritten;

  const McOptionsValues *                            m_optionsObj;

  bool m_userDidNotProvideOptions;
};

}  // End namespace QUESO

#endif // UQ_MOC_SG_H
