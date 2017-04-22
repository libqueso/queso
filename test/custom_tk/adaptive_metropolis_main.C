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

#include <queso/EnvironmentOptions.h>
#include <queso/GenericScalarFunction.h>
#include <queso/GslMatrix.h>
#include <queso/UniformVectorRV.h>
#include <queso/StatisticalInverseProblem.h>

template <class V = QUESO::GslVector, class M = QUESO::GslMatrix>
class Likelihood : public QUESO::BaseScalarFunction<V, M>
{
public:
  Likelihood(const char * prefix, const QUESO::VectorSet<V, M> & domainSet)
    : QUESO::BaseScalarFunction<V, M>(prefix, domainSet)
  {
  }

  ~Likelihood()
  {
  }

  virtual double lnValue(const V & domainVector, const V * domainDirection,
      V * gradVector, M * hessianMatrix, V * hessianEffect) const
  {
    return -0.5 * domainVector[0] * domainVector[0];
  }

  virtual double actualValue(const V & domainVector, const V * domainDirection,
      V * gradVector, M * hessianMatrix, V * hessianEffect) const
  {
    return std::exp(this->lnValue(domainVector, domainDirection, gradVector,
          hessianMatrix, hessianEffect));
  }

  double m_obsStdDev;
};

int main(int argc, char ** argv)
{
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc,&argv);
#endif

  QUESO::EnvOptionsValues options;
  options.m_numSubEnvironments = 1;
  options.m_subDisplayFileName = ".";
  options.m_subDisplayAllowAll = 0;
  options.m_subDisplayAllowedSet.insert(0);
  options.m_seed = 1.0;
  options.m_checkingLevel = 1;
  options.m_displayVerbosity = 0;

#ifdef QUESO_HAS_MPI
  QUESO::FullEnvironment env(MPI_COMM_WORLD, "", "", &options);
#else
  QUESO::FullEnvironment env("", "", &options);
#endif

  QUESO::VectorSpace<> paramSpace(env, "param_", 1, NULL);

  QUESO::GslVector paramMins(paramSpace.zeroVector());
  QUESO::GslVector paramMaxs(paramSpace.zeroVector());

  paramMins.cwSet(-INFINITY);
  paramMaxs.cwSet( INFINITY);

  QUESO::BoxSubset<> paramDomain("param_", paramSpace, paramMins, paramMaxs);

  Likelihood<> likelihood("", paramDomain);

  // Step 4 of 5: Instantiate the inverse problem
  QUESO::UniformVectorRV<> priorRv("prior_", paramDomain);
  QUESO::GenericVectorRV<> postRv("post_", paramSpace);
  QUESO::StatisticalInverseProblem<> ip("", NULL, priorRv, likelihood, postRv);

  // Step 5 of 5: Solve the inverse problem
  QUESO::GslVector paramInitials(paramSpace.zeroVector());

  QUESO::GslMatrix proposalCovMatrix(paramSpace.zeroVector());
  proposalCovMatrix(0,0) = 1e-6;

  QUESO::MhOptionsValues mh_options;
  mh_options.m_tk = "my_tk";
  mh_options.m_rawChainSize = 50;
  mh_options.m_amInitialNonAdaptInterval = 10;
  mh_options.m_amAdaptInterval = 10;
  mh_options.m_amAdaptedMatricesDataOutputFileName = "output_test_custom_tk_am/matrices";
  mh_options.m_amAdaptedMatricesDataOutputPeriod = 1;

  ip.solveWithBayesMetropolisHastings(&mh_options, paramInitials, &proposalCovMatrix);

  // Now read in test output and regression output and check they're the same
  QUESO::ScalarSequence<> computed_mat10(env, 1, "");
  QUESO::ScalarSequence<> computed_mat20(env, 1, "");
  QUESO::ScalarSequence<> computed_mat30(env, 1, "");
  QUESO::ScalarSequence<> computed_mat40(env, 1, "");

  QUESO::ScalarSequence<> expected_mat10(env, 1, "");
  QUESO::ScalarSequence<> expected_mat20(env, 1, "");
  QUESO::ScalarSequence<> expected_mat30(env, 1, "");
  QUESO::ScalarSequence<> expected_mat40(env, 1, "");

  computed_mat10.unifiedReadContents("output_test_custom_tk_am/matrices_am10_sub0", "m", 1);
  computed_mat20.unifiedReadContents("output_test_custom_tk_am/matrices_am20_sub0", "m", 1);
  computed_mat30.unifiedReadContents("output_test_custom_tk_am/matrices_am30_sub0", "m", 1);
  computed_mat40.unifiedReadContents("output_test_custom_tk_am/matrices_am40_sub0", "m", 1);

  std::string regressionsFileName = "custom_tk/output_test_custom_tk_am_regression";
  const char * regressions_srcdir = std::getenv("srcdir");
  if (regressions_srcdir) {
    regressionsFileName = regressions_srcdir + ('/' + regressionsFileName);
  }

  expected_mat10.unifiedReadContents(regressionsFileName + ('/' + std::string("matrices_am10_sub0")), "m", 1);
  expected_mat20.unifiedReadContents(regressionsFileName + ('/' + std::string("matrices_am20_sub0")), "m", 1);
  expected_mat30.unifiedReadContents(regressionsFileName + ('/' + std::string("matrices_am30_sub0")), "m", 1);
  expected_mat40.unifiedReadContents(regressionsFileName + ('/' + std::string("matrices_am40_sub0")), "m", 1);

  queso_require_equal_to_msg(computed_mat10[0], expected_mat10[0], "adapted matrix at iter 10 not equal");
  queso_require_equal_to_msg(computed_mat20[0], expected_mat20[0], "adapted matrix at iter 20 not equal");
  queso_require_equal_to_msg(computed_mat30[0], expected_mat30[0], "adapted matrix at iter 30 not equal");
  queso_require_equal_to_msg(computed_mat40[0], expected_mat40[0], "adapted matrix at iter 40 not equal");

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
