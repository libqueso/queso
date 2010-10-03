//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010 The PECOS Development Team
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

#ifndef __UQ_SFP_OPTIONS_H__
#define __UQ_SFP_OPTIONS_H__

#include <uqEnvironment.h>
//#include <uqMonteCarloSGOptions.h>

#undef UQ_SFP_READS_SOLVER_OPTION

#define UQ_SFP_FILENAME_FOR_NO_FILE "."

// _ODV = option default value
#define UQ_SFP_COMPUTE_SOLUTION_ODV        1
#define UQ_SFP_COMPUTE_COVARIANCES_ODV     1
#define UQ_SFP_COMPUTE_CORRELATIONS_ODV    1
#define UQ_SFP_DATA_OUTPUT_FILE_NAME_ODV   UQ_SFP_FILENAME_FOR_NO_FILE
#define UQ_SFP_DATA_OUTPUT_ALLOWED_SET_ODV ""
#ifdef UQ_SFP_READS_SOLVER_OPTION
#define UQ_SFP_SOLVER_ODV                  "mc" // Monte Carlo
#endif

class uqSfpOptionsValuesClass
{
public:
  uqSfpOptionsValuesClass            ();
  uqSfpOptionsValuesClass            (const uqSfpOptionsValuesClass& src);
  uqSfpOptionsValuesClass& operator= (const uqSfpOptionsValuesClass& rhs);
 ~uqSfpOptionsValuesClass            ();

  bool                   m_computeSolution;
  bool                   m_computeCovariances;
  bool                   m_computeCorrelations;
  std::string            m_dataOutputFileName;
  std::set<unsigned int> m_dataOutputAllowedSet;
#ifdef UQ_SFP_READS_SOLVER_OPTION
  std::string            m_solverString;
#endif

  //uqMcOptionsValuesClass m_mcOptionsValues;

private:
  void copy(const uqSfpOptionsValuesClass& src);
};

class uqStatisticalForwardProblemOptionsClass
{
public:
  uqStatisticalForwardProblemOptionsClass(const uqBaseEnvironmentClass& env, const char* prefix);
  uqStatisticalForwardProblemOptionsClass(const uqBaseEnvironmentClass& env, const char* prefix, const uqSfpOptionsValuesClass& alternativeOptionsValues);
 ~uqStatisticalForwardProblemOptionsClass();

  void scanOptionsValues();
  void print            (std::ostream& os) const;

  uqSfpOptionsValuesClass       m_ov;
  std::string                   m_prefix;

private:
  void   defineMyOptions  (po::options_description& optionsDesc) const;
  void   getMyOptionValues(po::options_description& optionsDesc);

  const uqBaseEnvironmentClass& m_env;

  po::options_description*      m_optionsDesc;
  std::string                   m_option_help;
  std::string                   m_option_computeSolution;
  std::string                   m_option_computeCovariances;
  std::string                   m_option_computeCorrelations;
  std::string                   m_option_dataOutputFileName;
  std::string                   m_option_dataOutputAllowedSet;
#ifdef UQ_SFP_READS_SOLVER_OPTION
  std::string                   m_option_solver;
#endif
};

std::ostream& operator<<(std::ostream& os, const uqStatisticalForwardProblemOptionsClass& obj);
#endif // __UQ_SFP_OPTIONS_H__
