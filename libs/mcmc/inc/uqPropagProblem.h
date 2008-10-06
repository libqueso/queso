/* uq/libs/mcmc/inc/uqPropagProblem.h
 *
 * Copyright (C) 2008 The PECOS Team, http://queso.ices.utexas.edu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_PROPAG_PROBLEM_H__
#define __UQ_PROPAG_PROBLEM_H__

#include <uqVectorFunction.h> // For substep 4 (or 1) in appls. with propagation

#include <uqMonteCarloSG.h>
#include <uqVectorRV.h>

#undef UQ_PROPAG_PROBLEM_READS_SOLVER_OPTION

#define UQ_PROPAG_PROBLEM_FILENAME_FOR_NO_OUTPUT_FILE "."

// _ODV = option default value
#define UQ_PROPAG_PROBLEM_COMPUTE_SOLUTION_ODV 1
#define UQ_PROPAG_PROBLEM_OUTPUT_FILE_NAME_ODV UQ_PROPAG_PROBLEM_FILENAME_FOR_NO_OUTPUT_FILE
#ifdef UQ_PROPAG_PROBLEM_READS_SOLVER_OPTION
#define UQ_PROPAG_PROBLEM_SOLVER_ODV           "mc" // Monte Carlo
#endif

template <class P_V,class P_M,class Q_V,class Q_M>
class uqPropagProblemClass
{
public:
  uqPropagProblemClass(const char*                                       prefix,
                       const uqBaseVectorRVClass      <P_V,P_M>&         paramRv,
                       const uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>& qoiFunction,
                             uqBaseVectorRVClass      <Q_V,Q_M>&         qoiRV);
 ~uqPropagProblemClass();

        void solveWithMonteCarlo();

        void print              (std::ostream& os) const;

private:
        void commonConstructor  ();
        void defineMyOptions    (po::options_description& optionsDesc);
        void getMyOptionValues  (po::options_description& optionsDesc);

  const uqEnvironmentClass&                         m_env;
        std::string                                 m_prefix;

        po::options_description*                    m_optionsDesc;
        std::string                                 m_option_help;
	std::string                                 m_option_computeSolution;
        std::string                                 m_option_outputFileName;
#ifdef UQ_PROPAG_PROBLEM_READS_SOLVER_OPTION
        std::string                                 m_option_solver;
#endif

        bool                                        m_computeSolution;
        std::string                                 m_outputFileName;
#ifdef UQ_PROPAG_PROBLEM_READS_SOLVER_OPTION
	std::string                                 m_solverString;
#endif

  const uqBaseVectorRVClass      <P_V,P_M>&         m_paramRv;
  const uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>& m_qoiFunction;
        uqBaseVectorRVClass      <Q_V,Q_M>&         m_qoiRv;

        uqBaseVectorPdfClass     <Q_V,Q_M>*         m_solutionPdf;
        uqBaseVectorMdfClass     <Q_V,Q_M>*         m_solutionMdf;
        uqBaseVectorCdfClass     <Q_V,Q_M>*         m_solutionCdf;
        uqBaseVectorRealizerClass<Q_V,Q_M>*         m_solutionRealizer;

        uqMonteCarloSGClass      <P_V,P_M,Q_V,Q_M>* m_mcSeqGenerator;
        uqBaseVectorSequenceClass<Q_V,Q_M>*         m_chain;
        uqArrayOfOneDGridsClass  <Q_V,Q_M>*         m_mdfGrids;
        uqArrayOfOneDTablesClass <Q_V,Q_M>*         m_mdfValues;
        uqArrayOfOneDGridsClass  <Q_V,Q_M>*         m_cdfGrids;
        uqArrayOfOneDTablesClass <Q_V,Q_M>*         m_cdfValues;
};

template<class P_V,class P_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const uqPropagProblemClass<P_V,P_M,Q_V,Q_M>& obj);

template <class P_V,class P_M,class Q_V,class Q_M>
uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::uqPropagProblemClass(
  const char*                                       prefix,
  const uqBaseVectorRVClass      <P_V,P_M>&         paramRv,
  const uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>& qoiFunction,
        uqBaseVectorRVClass      <Q_V,Q_M>&         qoiRv)
  :
  m_env                   (paramRv.env()),
  m_prefix                ((std::string)(prefix) + "pro_"),
  m_optionsDesc           (new po::options_description("UQ Propagation Problem")),
  m_option_help           (m_prefix + "help"           ),
  m_option_computeSolution(m_prefix + "computeSolution"),
  m_option_outputFileName (m_prefix + "outputFileName" ),
#ifdef UQ_PROPAG_PROBLEM_READS_SOLVER_OPTION
  m_option_solver         (m_prefix + "solver"         ),
#endif
  m_computeSolution       (UQ_PROPAG_PROBLEM_COMPUTE_SOLUTION_ODV),
  m_outputFileName        (UQ_PROPAG_PROBLEM_OUTPUT_FILE_NAME_ODV),
#ifdef UQ_PROPAG_PROBLEM_READS_SOLVER_OPTION
  m_solverString          (UQ_PROPAG_PROBLEM_SOLVER_ODV),
#endif
  m_paramRv               (paramRv),
  m_qoiFunction           (qoiFunction),
  m_qoiRv                 (qoiRv),
  m_solutionPdf           (NULL),
  m_solutionMdf           (NULL),
  m_solutionCdf           (NULL),
  m_solutionRealizer      (NULL),
  m_mcSeqGenerator        (NULL),
  m_chain                 (NULL)
{
  if (m_env.rank() == 0) std::cout << "Entering uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::constructor()"
                                   << ": prefix = "              << m_prefix
                                   << std::endl;

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.rank() == 0) std::cout << "In uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::constructor()"
                                   << ": after getting values of options, state of object is:"
                                   << "\n" << *this
                                   << std::endl;

  if (m_env.rank() == 0) std::cout << "Leaving uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::constructor()"
                                   << ": prefix = "              << m_prefix
                                   << std::endl;
}

template <class P_V,class P_M,class Q_V,class Q_M>
uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::~uqPropagProblemClass()
{
  if (m_chain) {
    m_chain->clear();
    delete m_chain;
  }
  if (m_mcSeqGenerator  ) delete m_mcSeqGenerator;
  if (m_solutionRealizer) delete m_solutionRealizer;
  if (m_solutionCdf     ) delete m_solutionCdf;
  if (m_solutionMdf     ) delete m_solutionMdf;
  if (m_solutionPdf     ) delete m_solutionPdf;
  if (m_optionsDesc     ) delete m_optionsDesc;
}

template<class P_V,class P_M,class Q_V,class Q_M>
void
uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::defineMyOptions(
  po::options_description& optionsDesc)
{
  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                                             "produce help message for propagation problem")
    (m_option_computeSolution.c_str(), po::value<bool       >()->default_value(UQ_PROPAG_PROBLEM_COMPUTE_SOLUTION_ODV), "compute solution process"                       )
    (m_option_outputFileName.c_str(),  po::value<std::string>()->default_value(UQ_PROPAG_PROBLEM_OUTPUT_FILE_NAME_ODV), "name of output file"                         )
#ifdef UQ_PROPAG_PROBLEM_READS_SOLVER_OPTION
    (m_option_solver.c_str(),          po::value<std::string>()->default_value(UQ_PROPAG_PROBLEM_SOLVER_ODV          ), "algorithm for propagation"                   )
#endif
  ;

  return;
}

template<class P_V,class P_M,class Q_V,class Q_M>
void
  uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::getMyOptionValues(
  po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    std::cout << optionsDesc
              << std::endl;
  }

  if (m_env.allOptionsMap().count(m_option_computeSolution.c_str())) {
    m_computeSolution = m_env.allOptionsMap()[m_option_computeSolution.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_outputFileName.c_str())) {
    m_outputFileName = m_env.allOptionsMap()[m_option_outputFileName.c_str()].as<std::string>();
  }

#ifdef UQ_PROPAG_PROBLEM_READS_SOLVER_OPTION
  if (m_env.allOptionsMap().count(m_option_solver.c_str())) {
    m_solverString = m_env.allOptionsMap()[m_option_solver.c_str()].as<std::string>();
  }
#endif

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::solveWithMonteCarlo()
{
  if (m_computeSolution == false) {
    if ((m_env.rank() == 0)) {
      std::cout << "In uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::solveWithMonteCarlo()"
                << ": computeping solution, as requested by user"
                << std::endl;
    }
    return;
  }

  if (m_solutionPdf     ) delete m_solutionPdf;
  if (m_solutionMdf     ) delete m_solutionMdf;
  if (m_solutionCdf     ) delete m_solutionCdf;
  if (m_solutionRealizer) delete m_solutionRealizer;
  if (m_mcSeqGenerator  ) delete m_mcSeqGenerator;

  // Compute output realizer: Monte Carlo approach
  m_chain = new uqSequenceOfVectorsClass<Q_V,Q_M>(m_qoiRv.imageSpace(),0,m_prefix+"chain");
  m_mcSeqGenerator = new uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>(m_prefix.c_str(),
                                                              m_paramRv,
                                                              m_qoiFunction,
                                                              m_qoiRv);
  m_mcSeqGenerator->generateSequence(*m_chain);
  m_solutionRealizer = new uqSequentialVectorRealizerClass<Q_V,Q_M>(m_prefix.c_str(),
                                                                   *m_chain);
  m_qoiRv.setRealizer(*m_solutionRealizer);

  // Compute output mdf: uniform sampling approach
  m_mdfGrids  = new uqArrayOfOneDGridsClass <Q_V,Q_M>((m_prefix+"mdf_").c_str(),m_qoiRv.imageSpace());
  m_mdfValues = new uqArrayOfOneDTablesClass<Q_V,Q_M>((m_prefix+"mdf_").c_str(),m_qoiRv.imageSpace());
  Q_V* numIntervalsVec = m_qoiRv.imageSpace().newVector(250.);
  m_chain->uniformlySampledMdf(*numIntervalsVec, // input
                               *m_mdfGrids,      // output
                               *m_mdfValues);    // output
  delete numIntervalsVec;
  m_solutionMdf = new uqSampledVectorMdfClass<Q_V,Q_M>(m_prefix.c_str(),
                                                      *m_mdfGrids,
                                                      *m_mdfValues);
  m_qoiRv.setMdf(*m_solutionMdf);

  // Compute output cdf: uniform sampling approach
  m_cdfGrids  = new uqArrayOfOneDGridsClass <Q_V,Q_M>((m_prefix+"cdf_").c_str(),m_qoiRv.imageSpace());
  m_cdfValues = new uqArrayOfOneDTablesClass<Q_V,Q_M>((m_prefix+"cdf_").c_str(),m_qoiRv.imageSpace());
  numIntervalsVec = m_qoiRv.imageSpace().newVector(250.);
  m_chain->uniformlySampledCdf(*numIntervalsVec, // input
                               *m_cdfGrids,      // output
                               *m_cdfValues);    // output
  delete numIntervalsVec;
  m_solutionCdf = new uqSampledVectorCdfClass<Q_V,Q_M>(m_prefix.c_str(),
                                                      *m_cdfGrids,
                                                      *m_cdfValues);
  m_qoiRv.setCdf(*m_solutionCdf);

  if (m_outputFileName != UQ_PROPAG_PROBLEM_FILENAME_FOR_NO_OUTPUT_FILE) {
    // Write output file
    if (m_env.rank() == 0) {
      std::cout << "Opening output file '" << m_outputFileName
                << "' for propagation problem with problem with prefix = " << m_prefix
                << std::endl;
    }

    // Open file
    std::ofstream* ofs = new std::ofstream(m_outputFileName.c_str(), std::ofstream::out | std::ofstream::in | std::ofstream::ate);
    if ((ofs            == NULL ) ||
        (ofs->is_open() == false)) {
      delete ofs;
      ofs = new std::ofstream(m_outputFileName.c_str(), std::ofstream::out | std::ofstream::trunc);
    }
    UQ_FATAL_TEST_MACRO((ofs && ofs->is_open()) == false,
                        m_env.rank(),
                        "uqPropagProblem<P_V,P_M,Q_V,Q_M>::solveWithBayesMarkovChain()",
                        "failed to open file");

    m_qoiRv.mdf().printContents(*ofs);
    m_qoiRv.cdf().printContents(*ofs);

    // Close file
    ofs->close();
    delete ofs;
    if (m_env.rank() == 0) {
      std::cout << "Closed output file '" << m_outputFileName
                << "' for propagation problem with problem with prefix = " << m_prefix
                << std::endl;
    }
  }
  if (m_env.rank() == 0) {
    std::cout << std::endl;
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqPropagProblemClass<P_V,P_M,Q_V,Q_M>::print(std::ostream& os) const
{
  os << "\n" << m_option_computeSolution << " = " << m_computeSolution
     << "\n" << m_option_outputFileName  << " = " << m_outputFileName
#ifdef UQ_PROPAG_PROBLEM_READS_SOLVER_OPTION
     << "\n" << m_option_solver          << " = " << m_solverString
#endif
     << std::endl;
}

template<class P_V,class P_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const uqPropagProblemClass<P_V,P_M,Q_V,Q_M>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_PROPAG_PROBLEM_H__
