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

#ifndef __UQ_SIP_OPTIONS_H__
#define __UQ_SIP_OPTIONS_H__

#include <uqEnvironment.h>
//#include <uqMetropolisHastingsSGOptions.h>

#undef UQ_SIP_READS_SOLVER_OPTION

#define UQ_SIP_FILENAME_FOR_NO_FILE "."

// _ODV = option default value
#define UQ_SIP_COMPUTE_SOLUTION_ODV        1
#define UQ_SIP_DATA_OUTPUT_FILE_NAME_ODV   UQ_SIP_FILENAME_FOR_NO_FILE
#define UQ_SIP_DATA_OUTPUT_ALLOWED_SET_ODV ""
#ifdef UQ_SIP_READS_SOLVER_OPTION
#define UQ_SIP_SOLVER_ODV                  "bayes_mc" // Bayesian formula + Metropolis-Hastings
#endif

class uqSipOptionsValuesClass
{
public:
  uqSipOptionsValuesClass            ();
  uqSipOptionsValuesClass            (const uqSipOptionsValuesClass& src);
  uqSipOptionsValuesClass& operator= (const uqSipOptionsValuesClass& rhs);
 ~uqSipOptionsValuesClass            ();

  bool                   m_computeSolution;
  std::string            m_dataOutputFileName;
  std::set<unsigned int> m_dataOutputAllowedSet;
#ifdef UQ_SIP_READS_SOLVER_OPTION
  std::string            m_solverString;
#endif

  //uqMhOptionsValuesClass m_mhOptionsValues;

private:
  void copy(const uqSipOptionsValuesClass& src);
};

class uqStatisticalInverseProblemOptionsClass
{
public:
  uqStatisticalInverseProblemOptionsClass(const uqBaseEnvironmentClass& env, const char* prefix);
  uqStatisticalInverseProblemOptionsClass(const uqBaseEnvironmentClass& env, const char* prefix, const uqSipOptionsValuesClass& alternativeOptionsValues);
 ~uqStatisticalInverseProblemOptionsClass();

  void scanOptionsValues();
  void print            (std::ostream& os) const;

  uqSipOptionsValuesClass       m_ov;
  std::string                   m_prefix;

private:
  void   defineMyOptions  (po::options_description& optionsDesc) const;
  void   getMyOptionValues(po::options_description& optionsDesc);

  const uqBaseEnvironmentClass& m_env;

  po::options_description*      m_optionsDesc;
  std::string                   m_option_help;
  std::string                   m_option_computeSolution;
  std::string                   m_option_dataOutputFileName;
  std::string                   m_option_dataOutputAllowedSet;
#ifdef UQ_SIP_READS_SOLVER_OPTION
  std::string                   m_option_solver;
#endif
};

std::ostream& operator<<(std::ostream& os, const uqStatisticalInverseProblemOptionsClass& obj);
#endif // __UQ_SIP_OPTIONS_H__
