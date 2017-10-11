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

#include <queso/GPMSA.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/SimulationOutputMesh.h>

namespace QUESO {

template <class V, class M>
GPMSAEmulator<V, M>::GPMSAEmulator(
    const VectorSet<V, M> & /* domain */,
    const VectorSpace<V, M> & m_scenarioSpace,
    const VectorSpace<V, M> & m_parameterSpace,
    const VectorSpace<V, M> & m_simulationOutputSpace,
    const unsigned int m_numSimulations,
    const unsigned int m_numExperiments,
    const std::vector<typename SharedPtr<V>::Type> & m_simulationScenarios,
    const std::vector<typename SharedPtr<V>::Type> & m_simulationParameters,
    const std::vector<typename SharedPtr<V>::Type> & m_simulationOutputs,
    const std::vector<typename SharedPtr<V>::Type> & m_experimentScenarios,
    const std::vector<typename SharedPtr<V>::Type> & m_experimentOutputs,
    const std::vector<typename SharedPtr<V>::Type> & m_discrepancyBases,
    const std::vector<typename SharedPtr<M>::Type> & m_observationErrorMatrices,
    const typename SharedPtr<M>::Type & m_observationErrorMatrix,
    const ConcatenatedVectorRV<V, M> & m_totalPrior,
    const V & residual_in,
    const M & BT_Wy_B_inv_in,
    const M & KT_K_inv_in,
    const GPMSAOptions & opts,
    const std::vector<typename SharedPtr<SimulationOutputMesh<V> >::Type> & m_simulationMeshes_in,
    const std::vector<std::vector<SimulationOutputPoint> > & m_experimentPoints_in,
    const std::vector<std::vector<unsigned int> > & m_experimentVariables_in)
  :
  BaseScalarFunction<V, M>("", m_totalPrior.imageSet()),
  m_scenarioSpace(m_scenarioSpace),
  m_parameterSpace(m_parameterSpace),
  m_simulationOutputSpace(m_simulationOutputSpace),
  m_numSimulations(m_numSimulations),
  m_numExperiments(m_numExperiments),
  m_simulationScenarios(m_simulationScenarios),
  m_simulationParameters(m_simulationParameters),
  m_simulationOutputs(m_simulationOutputs),
  m_experimentScenarios(m_experimentScenarios),
  m_experimentOutputs(m_experimentOutputs),
  m_discrepancyBases(m_discrepancyBases),
  m_observationErrorMatrices(m_observationErrorMatrices),
  m_observationErrorMatrix(m_observationErrorMatrix),
  m_totalPrior(m_totalPrior),
  residual(residual_in),
  BT_Wy_B_inv(BT_Wy_B_inv_in),
  KT_K_inv(KT_K_inv_in),
  m_opts(opts),
  m_simulationMeshes(m_simulationMeshes_in),
  m_experimentPoints(m_experimentPoints_in),
  m_experimentVariables(m_experimentVariables_in),
  m_numExperimentOutputs(0)
{
  queso_assert_greater(m_numSimulations, 0);

  queso_assert_equal_to
    (m_simulationOutputs[0]->map().Comm().NumProc(), 1);
}

template <class V, class M>
GPMSAEmulator<V, M>::~GPMSAEmulator()
{
  // heap allocations get deleted via ScopedPtr
}

template <class V, class M>
double
GPMSAEmulator<V, M>::lnValue(const V & domainVector,
                                       const V * /* domainDirection */,
                                       V * /* gradVector */,
                                       M * /* hessianMatrix */,
                                       V * /* hessianEffect */) const
{
  // Components of domainVector:
  // theta(1)                     // = "theta", "t" in Higdon et. al. 2008
  // theta(2)
  // ...
  // theta(dimParameter)
  // truncation_error_precision   // = "lambda_eta", only exists in
  //                                   vector case
  // emulator_precision(1)        // = "lambda_{w i}" in vector case
  // ...                          // = "lambda_eta" in scalar case,
  //                                   which is unrelated to the
  //                                   vector case lambda_eta
  // emulator_precision(num_svd_terms)
  // emulator_corr_strength(1)    // = "rho_{eta k}"
  // ...                          // dimScenario = "p_x", dimParameter = "p_t"
  // emulator_corr_strength(dimScenario + dimParameter)
  // discrepancy_precision(1)     // = "lambda_delta" in scalar case,
  //                              //   "lambda_{v i}" in vector
  // ...
  // discrepancy_precision(F)     // "F" == num_discrepancy_groups
  // discrepancy_corr_strength(1) // = "rho_{delta k}" in scalar case,
  // ...                          //   "rho_{v i}" in vector
  // discrepancy_corr_strength(F*dimScenario)
  // emulator_data_precision(1)   // = "small white noise",
                                  // "small ridge", "nugget"
  // observation_error_precision  // = "lambda_y", iff user-requested

  // Other variables:
  // m_numSimulations             // = "m"
  // m_numExperiments             // = "n"
  // m_simulationScenarios        // = "eta"
  // m_experimentScenarios        // = "y"
  // m_experimentErrors           // = "Sigma_y"; now obsoleted by "W_y"
  // numOutputs                   // = "n_eta"
  // m_numExperimentOutputs       // = "n_y" := sum(n_y_i)
  // dimScenario                  // = "p_x"
  // num_svd_terms                // = "p_eta"
  // num_discrepancy_bases        // = "p_delta"
  // m_TruncatedSVD_simulationOutputs  // = "K_eta"
  // covMatrix                    // = "Sigma_D" in scalar case,
  //                                   "Sigma_zhat" in vector
  // m_discrepancyMatrices        // = "D_i"
  // m_observationErrorMatrices   // = "W_i^{-1}"
  // m_observationErrorMatrix     // = "W_y^{-1}"
  // m_BMatrix                    // = "B"
  // num_discrepancy_groups       // = "F"
  // m_emulatorPrecisionShapeVec          // = "a_eta"
  // 1.0/m_emulatorPrecisionScaleVec      // = "b_eta"
  // m_observationalPrecisionShapeVec     // = "a_y"
  // 1.0/m_observationalPrecisionScaleVec // = "b_y"

  // Construct covariance matrix
  const unsigned int totalRuns = this->m_numExperiments + this->m_numSimulations;
  const unsigned int numSimulationOutputs = this->m_simulationOutputSpace.dimLocal();
  const unsigned int num_discrepancy_bases = m_discrepancyBases.size();
  const unsigned int residualSize = (numSimulationOutputs == 1) ?
    (this->m_numSimulations + this->m_numExperiments) :
    totalRuns * num_svd_terms + m_numExperiments * num_discrepancy_bases;

  const unsigned int first_multivariate_index = m_simulationMeshes.empty() ?
    0 : (m_simulationMeshes.back()->first_solution_index() +
         m_simulationMeshes.back()->n_outputs());
  const unsigned int n_multivariate_indices =
    numSimulationOutputs - first_multivariate_index;
  const unsigned int n_variables = m_simulationMeshes.size() + n_multivariate_indices;
  const unsigned int num_discrepancy_groups = n_variables;

  double prodScenario = 1.0;
  double prodParameter = 1.0;
  double prodDiscrepancy = 1.0;
  unsigned int dimScenario = (this->m_scenarioSpace).dimLocal();
  unsigned int dimParameter = (this->m_parameterSpace).dimLocal();

  // Length of prior+hyperprior inputs
  unsigned int dimSum = 1 +
                        (this->num_svd_terms < num_nonzero_eigenvalues) +
                        m_opts.m_calibrateObservationalPrecision +
                        num_svd_terms +
                        dimParameter +
                        dimParameter +
                        dimScenario +
                        num_discrepancy_groups +
                        (num_discrepancy_groups * dimScenario);  // yum

  // Offset for Sigma_eta equivalent in vector case
  const unsigned int offset1 = (numSimulationOutputs == 1) ?
    0 : m_numExperiments * num_discrepancy_bases;

  // Offset for Sigma_w in vector case
  const unsigned int offset1b = offset1 +
    m_numExperiments * num_svd_terms;

  // Offset for lambda_eta term in zhat covariance in vector case
  const unsigned int offset2 = (numSimulationOutputs == 1) ?
    0 : m_numExperiments * (num_discrepancy_bases + num_svd_terms);

  // This is cumbersome.  All I want is a matrix.
  const MpiComm & comm = domainVector.map().Comm();
  Map z_map(residualSize, 0, comm);
  M covMatrix(this->m_env, z_map, residualSize);

  typename SharedPtr<V>::Type domainVectorParameter
    (new V(*(this->m_simulationParameters[0])));
  for (unsigned int k = 0; k < dimParameter; k++) {
    queso_assert (!queso_isnan(domainVector[k]));
    (*domainVectorParameter)[k] = domainVector[k];
  }

  // This for loop is a disaster and could do with a *lot* of optimisation
  for (unsigned int i = 0; i < totalRuns; i++) {

     // Scenario and uncertain input variables, *not* normalized
     typename SharedPtr<V>::Type scenario1;
     typename SharedPtr<V>::Type parameter1;

    // Decide whether to do experiment part of the covariance matrix
    // Get i-th simulation-parameter pair
    if (i < this->m_numExperiments) {
      // Experiment scenario (known)
      scenario1 = (this->m_experimentScenarios)[i];

      // Experiment parameter (unknown)
      parameter1 = domainVectorParameter;
    }
    else {
      scenario1 =
        (this->m_simulationScenarios)[i-this->m_numExperiments];
      parameter1 =
        (this->m_simulationParameters)[i-this->m_numExperiments];
    }

    for (unsigned int j = 0; j < totalRuns; j++) {

       // Scenario and uncertain input variables, *not* normalized
       typename SharedPtr<V>::Type scenario2;
       typename SharedPtr<V>::Type parameter2;

      if (j < this->m_numExperiments) {
        scenario2 = (this->m_experimentScenarios)[j];
        parameter2 = domainVectorParameter;
      }
      else {
        scenario2 =
          (this->m_simulationScenarios)[j-this->m_numExperiments];
        parameter2 =
          (this->m_simulationParameters)[j-this->m_numExperiments];
      }

      // Emulator component       // = first term in (1)
      prodScenario = 1.0;
      unsigned int emulatorCorrStrStart =
        dimParameter + (this->num_svd_terms < num_nonzero_eigenvalues) + num_svd_terms;
      for (unsigned int k = 0; k < dimScenario; k++) {
        const double & emulator_corr_strength =
          domainVector[emulatorCorrStrStart+k];
        double scenario_param1 =
          m_opts.normalized_scenario_parameter(k, (*scenario1)[k]);
        double scenario_param2 =
          m_opts.normalized_scenario_parameter(k, (*scenario2)[k]);
        prodScenario *= std::pow(emulator_corr_strength,
                                 4.0 * (scenario_param1 - scenario_param2) *
                                       (scenario_param1 - scenario_param2));
      }

      queso_assert (!queso_isnan(prodScenario));

      // = second term in (1)
      prodParameter = 1.0;
      for (unsigned int k = 0; k < dimParameter; k++) {
        queso_assert (!queso_isnan(domainVector[emulatorCorrStrStart+dimScenario+k]));
        queso_assert (!queso_isnan((*parameter1)[k]));
        queso_assert (!queso_isnan((*parameter2)[k]));
        const double & emulator_corr_strength =
          domainVector[emulatorCorrStrStart+dimScenario+k];
        double uncertain_param1 =
          m_opts.normalized_uncertain_parameter(k, (*parameter1)[k]);
        double uncertain_param2 =
          m_opts.normalized_uncertain_parameter(k, (*parameter2)[k]);
        prodParameter *= std::pow(
            emulator_corr_strength,
            4.0 * (uncertain_param1 - uncertain_param2) *
                  (uncertain_param1 - uncertain_param2));
      }

      queso_assert (!queso_isnan(prodParameter));

      // Sigma_eta in scalar case,
      // [Sigma_u, Sigma_uw; Sigma_uw^T, Sigma_w] in vector case
      for (unsigned int basis = 0; basis != num_svd_terms; ++basis)
        {
          // coefficient in (1)
          // The relevant precision for Sigma_eta is lambda_eta; for
          // Sigma_uw etc. it's lambda_wi and we skip lambda_eta
          const double relevant_precision =
            domainVector[dimParameter + basis +
                         (this->num_svd_terms<num_nonzero_eigenvalues)];
          queso_assert_greater(relevant_precision, 0.0);

          const unsigned int stridei =
            (i < this->m_numExperiments) ? m_numExperiments : m_numSimulations;
          const unsigned int offseti =
            (i < this->m_numExperiments) ? offset1 : offset1b - m_numExperiments;
          const unsigned int stridej =
            (j < this->m_numExperiments) ? m_numExperiments : m_numSimulations;
          const unsigned int offsetj =
            (j < this->m_numExperiments) ? offset1 : offset1b - m_numExperiments;

          covMatrix(offseti+basis*stridei+i,
                    offsetj+basis*stridej+j) =
            prodScenario * prodParameter / relevant_precision;
        }

      // If we're in the experiment cross correlation part, need extra
      // foo: Sigma_delta/Sigma_v and Sigma_y
      if (i < this->m_numExperiments && j < this->m_numExperiments) {
        typename SharedPtr<V>::Type cross_scenario1 = (this->m_experimentScenarios)[i];
        typename SharedPtr<V>::Type cross_scenario2 = (this->m_experimentScenarios)[j];
        unsigned int discrepancyCorrStrStart =
          dimParameter + num_svd_terms + dimParameter + dimScenario + num_discrepancy_groups +
          (this->num_svd_terms<num_nonzero_eigenvalues);

        // Loop over discrepancy groups.  Keep track of which
        // submatrix we're on.
        unsigned int cov_matrix_offset = 0;
        for (unsigned int disc_grp = 0; disc_grp < num_discrepancy_groups; disc_grp++) {
          // If this is a functional case then we may have many
          // indices within this one discrepancy group
          const unsigned int disc_grp_size = (disc_grp < m_simulationMeshes.size()) ?
            m_simulationMeshes[disc_grp]->n_outputs() : 1;

          prodDiscrepancy = 1.0;
          for (unsigned int k = 0; k < dimScenario; k++) {
            const double & discrepancy_corr_strength =
              domainVector[discrepancyCorrStrStart+(disc_grp*dimScenario)+k];
            double cross_scenario_param1 =
              m_opts.normalized_scenario_parameter(k, (*cross_scenario1)[k]);
            double cross_scenario_param2 =
              m_opts.normalized_scenario_parameter(k, (*cross_scenario2)[k]);
            prodDiscrepancy *=
              std::pow(discrepancy_corr_strength, 4.0 *
                       (cross_scenario_param1 - cross_scenario_param2) *
                       (cross_scenario_param1 - cross_scenario_param2));
          }

          queso_assert (!queso_isnan(prodDiscrepancy));

          unsigned int discrepancyPrecisionStart = dimParameter +
                                                   (num_svd_terms<num_nonzero_eigenvalues) +
                                                   num_svd_terms +
                                                   dimScenario +
                                                   dimParameter;

          const double discrepancy_precision =
            domainVector[discrepancyPrecisionStart+disc_grp];

          queso_assert_greater(discrepancy_precision, 0);

          // Sigma_delta term from below (3) in univariate case
          // Sigma_v term from p. 576 in multivariate case
          const double R_v = prodDiscrepancy / discrepancy_precision;

          // In the functional case, we have many solution indices
          // corresponding to the same discrepancy group.
          for (unsigned int disc_grp_entry = 0; disc_grp_entry !=
               disc_grp_size; ++disc_grp_entry)
            {
              covMatrix(cov_matrix_offset+i,
                        cov_matrix_offset+j) += R_v;
              cov_matrix_offset += m_numExperiments;
            }
        }

        if (numSimulationOutputs == 1)
          {
            // Experimental error comes in via W_y now.
/*
            // Sigma_y term from below (3)
            const double experimentalError =
              this->m_experimentErrors(i,j);
*/
            const double experimentalError =
              (*this->m_observationErrorMatrix)(i,j);

            queso_assert_greater_equal (experimentalError, 0);

            const double lambda_y =
              m_opts.m_calibrateObservationalPrecision ?
              domainVector[dimSum-1] : 1.0;

            covMatrix(i,j) += experimentalError / lambda_y;
          }
      }
    }

    // Add small white noise component to diagonal to make stuff +ve def
    // = "small ridge"
    // Barely alluded to, never described, in Higdon et. al.
    const double emulator_data_precision =
      domainVector[dimSum-1-
                   m_opts.m_calibrateObservationalPrecision];

    queso_assert_greater(emulator_data_precision, 0);
    double nugget = 1.0 / emulator_data_precision;

    // Compute the offset occupied by the \Sigma_v matrix
    unsigned int discrepancy_offset =
      numSimulationOutputs == 1 ? 0 : num_discrepancy_bases;

    discrepancy_offset *= m_numExperiments;

    for (unsigned int basis = 0; basis < num_svd_terms; basis++) {
      covMatrix(discrepancy_offset+basis*totalRuns+i,
                discrepancy_offset+basis*totalRuns+i) += nugget;
    }
  }

  // If we're in the multivariate case, we've built the full Sigma_z
  // matrix; now add the remaining Sigma_zhat terms
  if (numSimulationOutputs > 1)
    {
      const double lambda_y =
        m_opts.m_calibrateObservationalPrecision ?
        domainVector[dimSum-1] : 1.0;
      const double inv_lambda_y = 1.0/lambda_y;

      unsigned int BT_Wy_B_size = BT_Wy_B_inv.numCols();
      for (unsigned int i=0; i != BT_Wy_B_size; ++i)
        for (unsigned int j=0; j != BT_Wy_B_size; ++j)
          covMatrix(i,j) += BT_Wy_B_inv(i,j) * inv_lambda_y;

      // Only add KT_K_inv if we're actually calibrating the truncation error
      // precision term
      if (num_svd_terms < num_nonzero_eigenvalues) {
        const double trunc_err_precision = domainVector[dimParameter];
        const double inv_trunc_err_precision = 1.0 / trunc_err_precision;

        unsigned int KT_K_size = KT_K_inv.numCols();
        for (unsigned int i=0; i != KT_K_size; ++i)
          for (unsigned int j=0; j != KT_K_size; ++j)
            covMatrix(i+offset2,j+offset2) +=
              KT_K_inv(i,j) * inv_trunc_err_precision;
      }
    }


  // Solve covMatrix * sol = residual
  // = Sigma_D^-1 * (D - mu 1) from (3)
  V sol(covMatrix.invertMultiply(residual));

  // Premultiply by residual^T as in (3)
  double minus_2_log_lhd = 0.0;
  // There's no dot product function in GslVector.
  for (unsigned int i = 0; i < residualSize; i++) {
    minus_2_log_lhd += sol[i] * residual[i];
  }

if (queso_isnan(minus_2_log_lhd))
  for (unsigned int i = 0; i < residualSize; i++) {
    if (queso_isnan(sol[i]))
      std::cout << "NaN sol[" << i << ']' << std::endl;
    if (queso_isnan(residual[i]))
      std::cout << "NaN residual[" << i << ']' << std::endl;

    std::cout << "Covariance Matrix:" << std::endl;
    covMatrix.print(std::cout);
  }

// std::cout << "minus_2_log_lhd = " << minus_2_log_lhd << std::endl;

  double cov_det = covMatrix.determinant();

  if (cov_det <= 0)
    {
      std::cout << "Non-positive determinant for covMatrix = " << std::endl;
      covMatrix.print(std::cout);
      queso_error();
    }

  queso_assert_greater(minus_2_log_lhd, 0);

  minus_2_log_lhd += std::log(covMatrix.determinant());

  // Multiply by -1/2 coefficient from (3)
  return -0.5 * minus_2_log_lhd;
}

template <class V, class M>
double
GPMSAEmulator<V, M>::actualValue(const V & /* domainVector */,
                                           const V * /* domainDirection */,
                                           V * /* gradVector */,
                                           M * /* hessianMatrix */,
                                           V * /* hessianEffect */) const
{
  // Do nothing?
  return 1.0;
}

template <class V, class M>
GPMSAFactory<V, M>::GPMSAFactory(
    const BaseEnvironment & env,
    GPMSAOptions * opts,
    const BaseVectorRV<V, M> & parameterPrior,
    const VectorSpace<V, M> & scenarioSpace,
    const VectorSpace<V, M> & parameterSpace,
    const VectorSpace<V, M> & simulationOutputSpace,
    unsigned int numSimulations,
    unsigned int numExperiments)
  :
    m_env(env),
    m_parameterPrior(parameterPrior),
    m_scenarioSpace(scenarioSpace),
    m_parameterSpace(parameterSpace),
    m_simulationOutputSpace(simulationOutputSpace),
    m_numSimulations(numSimulations),
    m_numExperiments(numExperiments),
    m_simulationScenarios(numSimulations),
    m_simulationParameters(numSimulations),
    m_simulationOutputs(numSimulations),
    m_experimentScenarios(numExperiments),
    m_experimentOutputs(numExperiments),
    m_numSimulationAdds(0),
    m_numExperimentAdds(0),
    num_svd_terms(0),
    num_nonzero_eigenvalues(0),
    priors(),
    m_constructedGP(false)
{
  // Set up a nonsense, "placeholder" discrepancy basis; we'll test
  // for this later to determine whether or not the user has
  // requested a non-default discrepancy basis.
  m_discrepancyBases.push_back(NULL);

  // DM: Not sure if the logic in these 3 if-blocks is correct
  if ((opts == NULL) && (this->m_env.optionsInputFileName() == "")) {
    queso_error_msg("Must options object or an input file");
  }

  if (opts != NULL) {
    allocated_m_opts = false;
    this->m_opts = opts;
  }
  else {
    // Create a default one
    allocated_m_opts = true;
    this->m_opts = new GPMSAOptions(this->m_env, "");
  }

  // FIXME: WTF? - RHS
  // this->m_constructedGP = false;
}

template <class V, class M>
GPMSAFactory<V, M>::~GPMSAFactory()
{
  if (this->allocated_m_opts)
    delete this->m_opts;
}

template <class V, class M>
const GPMSAOptions &
GPMSAFactory<V, M>::options() const
{
  return *this->m_opts;
}


template <class V, class M>
GPMSAOptions &
GPMSAFactory<V, M>::options()
{
  return *this->m_opts;
}


template <class V, class M>
unsigned int
GPMSAFactory<V, M>::numSimulations() const
{
  return this->m_numSimulations;
}

template <class V, class M>
unsigned int
GPMSAFactory<V, M>::numExperiments() const
{
  return this->m_numExperiments;
}

template <class V, class M>
const VectorSpace<V, M> &
GPMSAFactory<V, M>::scenarioSpace() const
{
  return this->m_scenarioSpace;
}

template <class V, class M>
const VectorSpace<V, M> &
GPMSAFactory<V, M>::parameterSpace() const
{
  return this->m_parameterSpace;
}

template <class V, class M>
const VectorSpace<V, M> &
GPMSAFactory<V, M>::simulationOutputSpace() const
{
  return this->m_simulationOutputSpace;
}

template <class V, class M>
const V &
GPMSAFactory<V, M>::simulationScenario(
    unsigned int simulationId) const
{
  queso_require_less_msg(simulationId, m_simulationScenarios.size(), "simulationId is too large");

  queso_require_msg(m_simulationScenarios[simulationId], "vector is NULL");

  return *(this->m_simulationScenarios[simulationId]);
}

template <class V, class M>
const std::vector<typename SharedPtr<V>::Type> &
GPMSAFactory<V, M>::simulationScenarios() const
{
  return this->m_simulationScenarios;
}

template <class V, class M>
const V &
GPMSAFactory<V, M>::simulationParameter(
    unsigned int simulationId) const
{
  queso_require_less_msg(simulationId, m_simulationParameters.size(), "simulationId is too large");

  queso_require_msg(m_simulationParameters[simulationId], "vector is NULL");

  return *(this->m_simulationParameters[simulationId]);
}

template <class V, class M>
const std::vector<typename SharedPtr<V>::Type> &
GPMSAFactory<V, M>::simulationParameters() const
{
  return this->m_simulationParameters;
}

template <class V, class M>
const V &
GPMSAFactory<V, M>::simulationOutput(
    unsigned int simulationId) const
{
  queso_require_less_msg(simulationId, m_simulationOutputs.size(), "simulationId is too large");

  queso_require_msg(m_simulationOutputs[simulationId], "vector is NULL");

  return *(this->m_simulationOutputs[simulationId]);
}

template <class V, class M>
const std::vector<typename SharedPtr<V>::Type> &
GPMSAFactory<V, M>::simulationOutputs() const
{
  return this->m_simulationOutputs;
}

template <class V, class M>
const V &
GPMSAFactory<V, M>::experimentScenario(
    unsigned int experimentId) const
{
  queso_require_less_msg(experimentId, (this->m_experimentScenarios).size(), "experimentId is too large");

  queso_require_msg(this->m_experimentScenarios[experimentId], "vector is NULL");

  return *(this->m_experimentScenarios[experimentId]);
}

template <class V, class M>
const std::vector<typename SharedPtr<V>::Type> &
GPMSAFactory<V, M>::experimentScenarios() const
{
  return this->m_experimentScenarios;
}

template <class V, class M>
const V &
GPMSAFactory<V, M>::experimentOutput(
    unsigned int experimentId) const
{
  queso_require_less_msg(experimentId, (this->m_experimentOutputs).size(), "experimentId is too large");

  queso_require_msg(this->m_experimentOutputs[experimentId], "vector is NULL");

  return *(this->m_experimentOutputs[experimentId]);
}

template <class V, class M>
const std::vector<typename SharedPtr<V>::Type> &
GPMSAFactory<V, M>::experimentOutputs() const
{
  return this->m_experimentOutputs;
}

template <class V, class M>
const BaseEnvironment &
GPMSAFactory<V, M>::env() const
{
  return this->m_env;
}

template <class V, class M>
const GPMSAEmulator<V, M> &
GPMSAFactory<V, M>::getGPMSAEmulator() const
{
  return *(this->gpmsaEmulator);
}

template <class V, class M>
void
GPMSAFactory<V, M>::addSimulationMesh
  (typename SharedPtr<SimulationOutputMesh<V> >::Type simulationMesh)
{
  // Make sure the new mesh solution coefficients don't overlap with
  // coefficients from any previously-added meshes
  if (!m_simulationMeshes.empty())
    {
      const SimulationOutputMesh<V> & mesh = *m_simulationMeshes.back();
      queso_require_equal_to
        (mesh.first_solution_index() + mesh.n_outputs(),
         simulationMesh->first_solution_index());
      queso_require_greater(mesh.n_outputs(), 0);
    }

  m_simulationMeshes.push_back(simulationMesh);
}

template <class V, class M>
void
GPMSAFactory<V, M>::addSimulation(typename SharedPtr<V>::Type simulationScenario,
                                  typename SharedPtr<V>::Type simulationParameter,
                                  typename SharedPtr<V>::Type simulationOutput)
{
  queso_require_less_msg(this->m_numSimulationAdds, this->m_numSimulations, "too many simulation adds...");

  this->m_simulationScenarios[this->m_numSimulationAdds] = simulationScenario;
  this->m_simulationParameters[this->m_numSimulationAdds] = simulationParameter;
  this->m_simulationOutputs[this->m_numSimulationAdds] = simulationOutput;
  this->m_numSimulationAdds++;

  if ((this->m_numSimulationAdds == this->m_numSimulations) &&
      (this->m_numExperimentAdds == this->m_numExperiments) &&
      (this->m_constructedGP == false)) {
    this->setUpEmulator();
  }
}

template <class V, class M>
void
GPMSAFactory<V, M>::addSimulations(
    const std::vector<typename SharedPtr<V>::Type> & simulationScenarios,
    const std::vector<typename SharedPtr<V>::Type> & simulationParameters,
    const std::vector<typename SharedPtr<V>::Type> & simulationOutputs)
{
  for (unsigned int i = 0; i < this->m_numSimulations; i++) {
    this->addSimulation(simulationScenarios[i],
                        simulationParameters[i],
                        simulationOutputs[i]);
  }
}



template <class V, class M>
void
GPMSAFactory<V, M>::setUpDiscrepancyBases()
{
  // If our "placeholder" basis is still there, then the user hasn't
  // requested anything different, and we need to autogenerate the
  // default basis.
  if ((m_discrepancyBases.size() == 1) &&
      !m_discrepancyBases[0].get())
    {
      const Map & output_map = m_simulationOutputs[0]->map();
      const BaseEnvironment &env = m_simulationOutputs[0]->env();

      // Replace our placeholder with the default discrepancy basis:
      m_discrepancyBases.clear();

      // Start by generating functional data gaussian discrepancy
      // bases, then generate multivariate delta-function discrepancy
      // bases.
      unsigned int first_multivariate_index = 0;
      for (unsigned int m=0; m != m_simulationMeshes.size(); ++m)
        {
          const SimulationOutputMesh<V> & mesh = *m_simulationMeshes[m];
          const unsigned int mesh_n_outputs = mesh.n_outputs();
          queso_assert_greater(mesh_n_outputs, 0);
          queso_assert_equal_to(mesh.first_solution_index(),
                                first_multivariate_index);
          first_multivariate_index += mesh_n_outputs;

          std::vector<typename SharedPtr<V>::Type> mesh_discrepancy_bases;
          mesh.generateDiscrepancyBases
            (this->options(), m, mesh_discrepancy_bases);
          m_discrepancyBases.insert (m_discrepancyBases.end(),
                                     mesh_discrepancy_bases.begin(),
                                     mesh_discrepancy_bases.end());
        }

      const unsigned int numSimulationOutputs =
        this->m_simulationOutputSpace.dimLocal();
      const unsigned int n_multivariate_indices =
        numSimulationOutputs - first_multivariate_index;
      const unsigned int n_variables = m_simulationMeshes.size() + n_multivariate_indices;
      const unsigned int variable_index_to_output_index = first_multivariate_index - m_simulationMeshes.size();
      for (unsigned int i=0; i != n_variables; ++i)
        {
          typename SharedPtr<V>::Type standard_basis(new V(env, output_map));
          // De-normalize the basis so it will be re-normalized later.
          (*standard_basis)[i+variable_index_to_output_index] =
            this->m_opts->output_scale(i);
          m_discrepancyBases.push_back(standard_basis);
        }
    }
}



template <class V, class M>
void
GPMSAFactory<V, M>::setUpEmulator()
{
  this->m_opts->template set_final_scaling<V>
    (m_simulationScenarios,
     m_simulationParameters,
     m_simulationOutputs,
     m_experimentScenarios,
     m_experimentOutputs,
     m_simulationMeshes);

  this->setUpDiscrepancyBases();

  const unsigned int numSimulationOutputs =
    this->m_simulationOutputSpace.dimLocal();

  const Map & output_map = m_simulationOutputs[0]->map();

  const MpiComm & comm = output_map.Comm();

  Map serial_map(m_numSimulations, 0, comm);

  const BaseEnvironment &env = m_simulationOutputs[0]->env();

  simulationOutputMeans.reset
    (new V (env, output_map));

  for (unsigned int i=0; i != m_numSimulations; ++i)
    for (unsigned int j=0; j != numSimulationOutputs; ++j)
      (*simulationOutputMeans)[j] +=
        this->m_opts->normalized_output(j, (*m_simulationOutputs[i])[j]);

  for (unsigned int j=0; j != numSimulationOutputs; ++j)
    (*simulationOutputMeans)[j] /= m_numSimulations;

  // Each group of functional data will use a *single* mean
  for (unsigned int m = 0; m != this->m_simulationMeshes.size(); ++m)
    {
      const SimulationOutputMesh<V> & mesh = *(this->m_simulationMeshes[m]);
      const unsigned int begin = mesh.first_solution_index();
      const unsigned int end = begin + mesh.n_outputs();
      double mean = 0;
      for (unsigned int i = begin; i != end; ++i)
        mean += (*simulationOutputMeans)[i];
      mean /= (end - begin);

      for (unsigned int i = begin; i != end; ++i)
        (*simulationOutputMeans)[i] = mean;
    }

  M simulation_matrix(env, serial_map, numSimulationOutputs);

  for (unsigned int i=0; i != m_numSimulations; ++i)
    for (unsigned int j=0; j != numSimulationOutputs; ++j)
      simulation_matrix(i,j) =
        this->m_opts->normalized_output(j, (*m_simulationOutputs[i])[j]) -
        (*simulationOutputMeans)[j];

  // GSL only finds left singular vectors if n_rows>=n_columns, so we need to
  // calculate them indirectly from the eigenvalues of M^T*M

  M S_trans(simulation_matrix.transpose());

  M SM_squared(S_trans*simulation_matrix);

  M SM_singularVectors(env, SM_squared.map(), numSimulationOutputs);
  V SM_singularValues(env, SM_squared.map());

  SM_squared.eigen(SM_singularValues, &SM_singularVectors);

  // Check the eigenvalues are in ascending order
  for (unsigned int i = 0; i < numSimulationOutputs-1; i++) {
    queso_assert_less_equal(SM_singularValues[i], SM_singularValues[i+1]);
  }

  // Count the number of "nonzero" eigenvalues (eigenvalues over a tolerance)
  double eigenvalue_tolerance = 1e-10;
  for (unsigned int i = 0; i < SM_singularValues.sizeLocal(); i++) {
    // All eigenvalues should be positive, so we don't take fabs
    if (SM_singularValues[i] > eigenvalue_tolerance) {
      num_nonzero_eigenvalues++;
    }
  }

  // Is this the right way to enforce the scalar case of num_svd_terms
  if (numSimulationOutputs > 1) {
    this->num_svd_terms = this->m_opts->m_maxEmulatorBasisVectors ?
      std::min((unsigned int)(this->m_opts->m_maxEmulatorBasisVectors),
          num_nonzero_eigenvalues) : num_nonzero_eigenvalues;
  }
  else {
    // We do this so logic like 'num_svd_terms < num_nonzero_eigenvalues' is
    // false in the scalar case.
    this->num_svd_terms = 1;
    this->num_nonzero_eigenvalues = 1;
  }

  queso_require_greater_equal(num_svd_terms, 0);

  // Copy only those vectors we want into K_eta
  m_TruncatedSVD_simulationOutputs.resize(num_svd_terms, V(env, output_map));

  // The singular values are in ascending order (with associated vectors ordered
  // accordingly).  Therefore we want to pull the singular vectors from the
  // back, not the front.
  for (unsigned int k = 0; k != num_svd_terms; ++k)
    m_TruncatedSVD_simulationOutputs[k] =
      SM_singularVectors.getColumn(numSimulationOutputs-1-k);

  Map copied_map(numSimulationOutputs * m_numSimulations, 0, comm);

  K.reset
    (new M(env, copied_map, m_numSimulations * num_svd_terms));
  for (unsigned int k=0; k != num_svd_terms; ++k)
    for (unsigned int i1=0; i1 != m_numSimulations; ++i1)
      for (unsigned int i2=0; i2 != numSimulationOutputs; ++i2)
        {
          const unsigned int i = i1 * numSimulationOutputs + i2;
          const unsigned int j = k * m_numSimulations + i1;
          (*K)(i,j) = m_TruncatedSVD_simulationOutputs[k][i2];
        }

  KT_K_inv.reset
    (new M((K->transpose() * *K).inverse()));

  // Create the giant matrices!

  // We build B from two parts, the diag(D_i)*P_D^T block and
  // the diag(K_i)*P_K^T block, but we build simultaneously so we only
  // need one loop over experiments.
  // Since P_D and P_K are permutation matrices we'll just apply the
  // permutations as we insert.

  // W_y is simple block diagonal.

  // Neither of these ought to be dense matrices but we're sacrificing
  // efficiency for clarity for now.

  const unsigned int num_discrepancy_bases = m_discrepancyBases.size();
  unsigned int Brows = 0;
  for (unsigned int i=0; i != m_numExperiments; ++i)
    Brows += m_experimentOutputs[i]->sizeLocal();
  const unsigned int Bcols =
    m_numExperiments * (num_discrepancy_bases + num_svd_terms);

  const Map B_row_map(Brows, 0, comm);

  m_BMatrix.reset
    (new M(env, B_row_map, Bcols));

  const Map Wy_row_map(Brows, 0, comm);

  M& B = *m_BMatrix;

  // Observation precision matrix
  M Wy(env, Wy_row_map, Brows);

  m_observationErrorMatrix.reset(new M(env, Wy_row_map, Brows));

  const unsigned int first_multivariate_index =
    m_simulationMeshes.empty() ?  0 :
    (m_simulationMeshes.back()->first_solution_index() +
     m_simulationMeshes.back()->n_outputs());
  const unsigned int variable_index_to_output_index =
    first_multivariate_index - m_simulationMeshes.size();

  // i = the current experimental output, used as a row index
  // This is easier iterated than calculated when each experiment
  // might have a different number of outputs.
  for (unsigned int ex = 0, i = 0; ex != m_numExperiments; ++ex)
    {
      const unsigned int numExperimentOutputs =
        m_experimentOutputs[ex]->sizeGlobal();
      Map serial_output_map(numExperimentOutputs, 0, comm);

      M D_i(env, serial_output_map,
            (unsigned int)(m_discrepancyBases.size()));

      for (unsigned int j=0; j != numExperimentOutputs; ++j)
        for (unsigned int k=0; k != m_discrepancyBases.size(); ++k)
          {
            const V & discrepancy_basis = *m_discrepancyBases[k];
            // If there's just multivariate data here, then no
            // interpolation is necessary.
            double discrepancy_basis_value = discrepancy_basis[j];

            if (this->m_experimentVariables.size())
              {
                const unsigned int var = this->m_experimentVariables[ex][j];
                const SimulationOutputPoint & pt = this->m_experimentPoints[ex][j];

                if (var < this->m_simulationMeshes.size())
                  {
                    SimulationOutputMesh<V> & mesh =
                      *this->m_simulationMeshes[var];
                    discrepancy_basis_value =
                      mesh.interpolateOutput(discrepancy_basis, pt);
                  }
                else
                  discrepancy_basis_value =
                    discrepancy_basis[var + variable_index_to_output_index];
              }

            D_i(j,k) = discrepancy_basis_value /
                       this->m_opts->output_scale(j);
          }

      const M & Sigma_i = *m_observationErrorMatrices[ex];

#ifndef NDEBUG
      // Each error covariance matrix had better be SPD
      M Sigma_i_copy(Sigma_i);
      int rv = Sigma_i_copy.chol();
      queso_assert_msg(!rv, "Observation error matrix Sigma_" << ex <<
                       " was not SPD!");
#endif

      // I think this computes the inverse from the LU factorisation.  Perhaps
      // we should have it done through the cholesky factorisation.
      //
      // inverse() returns a new matrix, so a reference is not enough here.
      M W_i(Sigma_i.inverse());

      // For the multivariate case, the bases K_eta computed from
      // simulator outputs are the same as the bases K_i which apply
      // to quantities of interest, because simulator outputs are QoIs
      // alone.

      // Wj = the current experimental output, used as a col index for
      // cross-correlations
      for (unsigned int outi = 0, Wj = i;
           outi != m_experimentOutputs[ex]->sizeLocal(); ++outi, ++i)
        {
          for (unsigned int outj = 0; outj != num_discrepancy_bases; ++outj)
            {
              unsigned int j = ex + m_numExperiments * outj;

              B(i,j) = D_i(outi,outj);
            }

          for (unsigned int outj = 0; outj != num_svd_terms; ++outj)
            {
              unsigned int j = ex +
                m_numExperiments * (num_discrepancy_bases + outj);

              V & singular_vector = m_TruncatedSVD_simulationOutputs[outj];
              // If there's just multivariate data here, then no
              // interpolation is necessary.
              double Kvalue = singular_vector[outi];

              if (this->m_experimentVariables.size())
                {
                  const unsigned int var = this->m_experimentVariables[ex][outi];
                  const SimulationOutputPoint & pt = this->m_experimentPoints[ex][outi];

                  if (var < this->m_simulationMeshes.size())
                    {
                      SimulationOutputMesh<V> & mesh =
                        *this->m_simulationMeshes[var];
                      Kvalue =

                        mesh.interpolateOutput(singular_vector, pt);
                    }
                  else
                    Kvalue =
                      singular_vector[var + variable_index_to_output_index];
                }

              B(i,j) = Kvalue;
            }

          // No fancy perturbation for Wj, just make sure it gets
          // reset properly
          for (unsigned int outj = 0;
               outj != m_experimentOutputs[ex]->sizeLocal();
               ++outj, ++Wj)
            {
              Wy(i,Wj) = W_i(outi,outj) *
                (this->m_opts->output_scale(outi) *
                 this->m_opts->output_scale(outj));

              (*m_observationErrorMatrix)(i,Wj) = Sigma_i(outi,outj) /
                (this->m_opts->output_scale(outi) *
                 this->m_opts->output_scale(outj));

            }
          Wj -= m_experimentOutputs[ex]->sizeLocal();
        }
    }

  M BT_Wy_B (B.transpose() * Wy * B);

  // Adding a "small ridge" to make sure this is invertible, as on
  // p.577 - defaulted to 1e-4 from discussion notes, but may be
  // overridden by user options.
  //
  // This is "DKridge" in the MATLAB prototype
  //
  // The number of rows BT_Wy_B has is Bcols
  for (unsigned int i=0; i != Bcols; ++i)
    BT_Wy_B(i,i) +=
      this->m_opts->m_observationalPrecisionRidge;

  BT_Wy_B_inv.reset(new M(BT_Wy_B.inverse()));

  // Add a ridge to the inverse even, if the user requested one.
  //
  // The number of rows BT_Wy_B_inv has is Bcols
  for (unsigned int i=0; i != Bcols; ++i)
    (*BT_Wy_B_inv)(i,i) +=
      this->m_opts->m_observationalCovarianceRidge;

  queso_assert_equal_to(BT_Wy_B_inv->numCols(), BT_Wy_B_inv->numRowsGlobal());
  queso_assert_equal_to(BT_Wy_B_inv->numCols(), Bcols);

  this->setUpHyperpriors(Wy);

  this->m_constructedGP = true;
  this->gpmsaEmulator.reset
    (new GPMSAEmulator<V, M>(
      this->prior().imageSet(),
      this->m_scenarioSpace,
      this->m_parameterSpace,
      this->m_simulationOutputSpace,
      this->m_numSimulations,
      this->m_numExperiments,
      this->m_simulationScenarios,
      this->m_simulationParameters,
      this->m_simulationOutputs,
      this->m_experimentScenarios,
      this->m_experimentOutputs,
      this->m_discrepancyBases,
      this->m_observationErrorMatrices,
      this->m_observationErrorMatrix,
      *(this->m_totalPrior),
      *this->residual,
      *this->BT_Wy_B_inv,
      *this->KT_K_inv,
      *this->m_opts,
      this->m_simulationMeshes,
      this->m_experimentPoints,
      this->m_experimentVariables));

  // FIXME: epic hack.  pls fix.
  this->gpmsaEmulator->num_svd_terms = num_svd_terms;
  this->gpmsaEmulator->num_nonzero_eigenvalues = num_nonzero_eigenvalues;
}



template <class V, class M>
void
GPMSAFactory<V, M>::addExperiments(
    const std::vector<typename SharedPtr<V>::Type> & experimentScenarios,
    const std::vector<typename SharedPtr<V>::Type> & experimentOutputs,
    const typename SharedPtr<M>::Type experimentErrors)
{
  queso_deprecated();

  queso_require_less_equal_msg(experimentScenarios.size(), this->m_numExperiments, "too many experiments...");

  unsigned int offset = 0;
  for (unsigned int i = 0; i < this->m_experimentScenarios.size(); i++) {
    this->m_experimentScenarios[i] = experimentScenarios[i];
    this->m_experimentOutputs[i] = experimentOutputs[i];

    const unsigned int outsize =
      this->m_experimentOutputs[i]->sizeGlobal();

    const BaseEnvironment & output_env = this->m_experimentOutputs[i]->env();
    const Map & output_map = this->m_experimentOutputs[i]->map();
    typename SharedPtr<M>::Type new_matrix(new M(output_env, output_map, 0.0));
    m_observationErrorMatrices.push_back(new_matrix);

    for (unsigned int outi = 0; outi != outsize; ++outi)
      for (unsigned int outj = 0; outj != outsize; ++outj)
        (*this->m_observationErrorMatrices[i])(outi,outj) =
          (*experimentErrors)(offset+outi, offset+outj);

    offset += outsize;
  }

  this->m_numExperimentAdds += experimentScenarios.size();

  if ((this->m_numSimulationAdds == this->m_numSimulations) &&
      (this->m_numExperimentAdds == this->m_numExperiments) &&
      (this->m_constructedGP == false)) {
    this->setUpEmulator();
  }
}


template <class V, class M>
void
GPMSAFactory<V, M>::addExperiments(
    const std::vector<typename SharedPtr<V>::Type> & experimentScenarios,
    const std::vector<typename SharedPtr<V>::Type> & experimentOutputs,
    const std::vector<typename SharedPtr<M>::Type> & experimentErrors,
    const std::vector<std::vector<SimulationOutputPoint> > * experimentPoints,
    const std::vector<std::vector<unsigned int> > * experimentVariables)
{
  queso_require_less_equal_msg(experimentScenarios.size(), this->m_numExperiments, "too many experiments...");
  queso_require_equal_to(experimentScenarios.size(),
                         experimentOutputs.size());
  queso_require_equal_to(experimentScenarios.size(),
                         experimentErrors.size());

  if (experimentPoints)
    {
      queso_require_equal_to(experimentScenarios.size(),
                             experimentPoints->size());
      queso_require(experimentVariables);
      queso_require_equal_to(experimentScenarios.size(),
                             experimentVariables->size());

      this->m_experimentPoints = *experimentPoints;
      this->m_experimentVariables = *experimentVariables;
    }

  for (unsigned int i = 0; i != this->m_experimentScenarios.size(); ++i)
    {
      queso_require_equal_to(experimentOutputs[i]->sizeGlobal(),
                             experimentErrors[i]->numCols());
      queso_require_equal_to(experimentOutputs[i]->sizeGlobal(),
                             experimentErrors[i]->numRowsGlobal());
      if (experimentPoints)
        {
          queso_require_equal_to(experimentOutputs[i]->sizeGlobal(),
                                 (*experimentPoints)[i].size());
          queso_require_equal_to(experimentOutputs[i]->sizeGlobal(),
                                 (*experimentVariables)[i].size());
        }
    }

  this->m_observationErrorMatrices.resize
    (this->m_experimentScenarios.size());

  for (unsigned int i = 0; i < this->m_experimentScenarios.size(); i++) {
    this->m_experimentScenarios[i] = experimentScenarios[i];
    this->m_experimentOutputs[i] = experimentOutputs[i];
    this->m_observationErrorMatrices[i] = experimentErrors[i];
  }
  this->m_numExperimentAdds += experimentScenarios.size();

  if ((this->m_numSimulationAdds == this->m_numSimulations) &&
      (this->m_numExperimentAdds == this->m_numExperiments) &&
      (this->m_constructedGP == false)) {
    this->setUpEmulator();
  }
}


template <class V, class M>
void
GPMSAFactory<V, M>::setDiscrepancyBases(
    const std::vector<typename SharedPtr<V>::Type> & discrepancyBases)
{
  m_discrepancyBases = discrepancyBases;

  // We should not yet have constructed the underlying GP model
  queso_assert_equal_to(this->m_constructedGP, false);
}



template <class V, class M>
M &
GPMSAFactory<V, M>::getObservationErrorCovariance
  (unsigned int simulationNumber)
{
  queso_assert_less(simulationNumber, m_numSimulations);
  queso_assert_equal_to(m_observationErrorMatrices.size(), m_numSimulations);

  return *m_observationErrorMatrices[simulationNumber];
}


template <class V, class M>
const M &
GPMSAFactory<V, M>::getObservationErrorCovariance
  (unsigned int simulationNumber) const
{
  queso_assert_less(simulationNumber, m_numSimulations);
  queso_assert_equal_to(m_observationErrorMatrices.size(), m_numSimulations);

  return *m_observationErrorMatrices[simulationNumber];
}



template <class V, class M>
const ConcatenatedVectorRV<V, M> &
GPMSAFactory<V, M>::prior() const
{
  return *(this->m_totalPrior);
}

template <class V, class M>
void
GPMSAFactory<V, M>::print(std::ostream& /* os */) const
{
  // Do nothing
}

// Private methods follow
template <class V, class M>
void
GPMSAFactory<V, M>::setUpHyperpriors(const M & Wy)
{
  const unsigned int numSimulationOutputs =
    this->m_simulationOutputSpace.dimLocal();

  const unsigned int num_discrepancy_bases = m_discrepancyBases.size();

  const unsigned int first_multivariate_index = m_simulationMeshes.empty() ?
    0 : (m_simulationMeshes.back()->first_solution_index() +
         m_simulationMeshes.back()->n_outputs());
  const unsigned int n_multivariate_indices =
    numSimulationOutputs - first_multivariate_index;
  const unsigned int n_variables = m_simulationMeshes.size() + n_multivariate_indices;
  const unsigned int num_discrepancy_groups = n_variables;

  const MpiComm & comm = m_simulationOutputs[0]->map().Comm();

  unsigned int rank_B;
  if (m_BMatrix->numRowsGlobal() > m_BMatrix->numCols())
    rank_B = m_BMatrix->rank(0, 1.e-4);
  else
    rank_B = m_BMatrix->transpose().rank(0, 1.e-4);

  double truncationErrorPrecisionShape = this->m_opts->m_truncationErrorPrecisionShape;
  double truncationErrorPrecisionScale = this->m_opts->m_truncationErrorPrecisionScale;

  double observationalPrecisionShape = this->m_opts->m_observationalPrecisionShape;
  double observationalPrecisionScale = this->m_opts->m_observationalPrecisionScale;

  if (numSimulationOutputs > 1)
    {
      unsigned int totalExperimentOutputs = 0;
      for (unsigned int i=0; i != m_numExperiments; ++i)
        totalExperimentOutputs += m_experimentOutputs[i]->sizeLocal();

      Map y_map(totalExperimentOutputs, 0, comm);
      Map eta_map(m_numSimulations * numSimulationOutputs, 0, comm);

      const unsigned int yhat_size =
        m_numExperiments * (num_discrepancy_bases + num_svd_terms);

      const unsigned int etahat_size =
        m_numSimulations * num_svd_terms;

      Map zhat_map(yhat_size + etahat_size, 0, comm);

      V y(this->m_env, y_map);
      V eta(this->m_env, eta_map);

      for (unsigned int i = 0, yindex = 0; i < this->m_numExperiments; i++) {
        for (unsigned int k = 0; k != this->m_experimentOutputs[i]->sizeLocal(); ++k)
          // FIXME - this won't normalize properly in the functional
          // case
          y[yindex++] =
            this->m_opts->normalized_output(k, (*((this->m_experimentOutputs)[i]))[k]) -
            (*simulationOutputMeans)[k];
      }

      for (unsigned int i = 0; i < this->m_numSimulations; i++) {
        for (unsigned int k = 0; k != numSimulationOutputs; ++k)
          eta[i*numSimulationOutputs+k] =
            this->m_opts->normalized_output(k, (*((this->m_simulationOutputs)[i]))[k]) -
            (*simulationOutputMeans)[k];
      }

      M& B = *m_BMatrix;

      V yhat(*BT_Wy_B_inv * (B.transpose() * (Wy * y)));

      queso_assert_equal_to(yhat.sizeGlobal(), yhat_size);

      V etahat(*KT_K_inv * (K->transpose() * eta));

      residual.reset(new V(this->m_env, zhat_map));
      for (unsigned int i = 0; i < yhat_size; ++i)
        (*residual)[i] = yhat[i];

      for (unsigned int i = 0; i < etahat_size; ++i)
        (*residual)[yhat_size+i] = etahat[i];

      truncationErrorPrecisionShape +=
        (this->m_numSimulations * (numSimulationOutputs - num_svd_terms)) / 2.0;

      V eta_temp(eta);
      eta_temp -= *K * etahat;

      double truncationErrorPrecisionRate = 1.0 / truncationErrorPrecisionScale;
      truncationErrorPrecisionRate +=
        scalarProduct(eta, eta_temp) / 2.0;
      truncationErrorPrecisionScale = 1.0 / truncationErrorPrecisionRate;

      observationalPrecisionShape +=
        (totalExperimentOutputs - rank_B) / 2.0;

      V y_temp(Wy * y);
      y_temp -= Wy * B * yhat;

      // At this point everything has already been normalized, so no
      // more normalization to do?
      double observationalPrecisionRate = 1.0 / observationalPrecisionScale;
      observationalPrecisionRate +=
        scalarProduct(y, y_temp) / 2.0;
      observationalPrecisionScale = 1.0 / observationalPrecisionRate;
    }
  else
    {
      const unsigned int totalRuns = this->m_numExperiments + this->m_numSimulations;
      Map z_map(totalRuns, 0, comm);
      residual.reset(new V (this->m_env, z_map));

      // Form residual = D - mean // = D - mu*1 in (3)
      // We currently use the mean of the simulation data, not a free
      // hyperparameter mean
      for (unsigned int i = 0; i < this->m_numExperiments; i++) {
        (*residual)[i] =
          this->m_opts->normalized_output(0, (*((this->m_experimentOutputs)[i]))[0]) -
          (*simulationOutputMeans)[0];
      }
      for (unsigned int i = 0; i < this->m_numSimulations; i++) {
        (*residual)[i+this->m_numExperiments] =
          this->m_opts->normalized_output(0, (*((this->m_simulationOutputs)[i]))[0]) -
          (*simulationOutputMeans)[0];
      }
    }


  double emulatorCorrelationStrengthAlpha = this->m_opts->m_emulatorCorrelationStrengthAlpha;
  double emulatorCorrelationStrengthBeta = this->m_opts->m_emulatorCorrelationStrengthBeta;
  double discrepancyPrecisionShape = this->m_opts->m_discrepancyPrecisionShape;
  double discrepancyPrecisionScale = this->m_opts->m_discrepancyPrecisionScale;
  double discrepancyCorrelationStrengthAlpha = this->m_opts->m_discrepancyCorrelationStrengthAlpha;
  double discrepancyCorrelationStrengthBeta = this->m_opts->m_discrepancyCorrelationStrengthBeta;
  double emulatorDataPrecisionShape = this->m_opts->m_emulatorDataPrecisionShape;
  double emulatorDataPrecisionScale = this->m_opts->m_emulatorDataPrecisionScale;

  this->oneDSpace.reset
    (new VectorSpace<V, M>(this->m_env, "", 1, NULL));

  this->numDiscrepancyGroupsDSpace.reset
    (new VectorSpace<V, M>(this->m_env, "", num_discrepancy_groups, NULL));

  // Truncation error precision
  if (this->num_svd_terms < num_nonzero_eigenvalues)
    {
      this->truncationErrorPrecisionSpace.reset
        (new VectorSpace<V, M>
         (this->m_env,
          "",
          1,
          NULL));

      this->truncationErrorPrecisionMin.reset(new V(this->oneDSpace->zeroVector()));
      this->truncationErrorPrecisionMax.reset(new V(this->oneDSpace->zeroVector()));
      this->m_truncationErrorPrecisionShapeVec.reset
        (new V(this->truncationErrorPrecisionSpace->zeroVector()));
      this->m_truncationErrorPrecisionScaleVec.reset
        (new V(this->truncationErrorPrecisionSpace->zeroVector()));
      this->truncationErrorPrecisionMin->cwSet(0);
      this->truncationErrorPrecisionMax->cwSet(INFINITY);

      this->m_truncationErrorPrecisionShapeVec->cwSet(truncationErrorPrecisionShape);
      this->m_truncationErrorPrecisionScaleVec->cwSet(truncationErrorPrecisionScale);

      this->truncationErrorPrecisionDomain.reset
        (new BoxSubset<V, M>
          ("",
           *(this->oneDSpace),
           *(this->truncationErrorPrecisionMin),
           *(this->truncationErrorPrecisionMax)));

      this->m_truncationErrorPrecision.reset
        (new GammaVectorRV<V, M>
         ("",
          *(this->truncationErrorPrecisionDomain),
          *(this->m_truncationErrorPrecisionShapeVec),
          *(this->m_truncationErrorPrecisionScaleVec)));
    }

  // Emulator precision
  this->emulatorPrecisionSpace.reset
    (new VectorSpace<V, M>
     (this->m_env,
      "",
      num_svd_terms,
      NULL));

  this->emulatorPrecisionMin.reset
    (new V(this->emulatorPrecisionSpace->zeroVector()));
  this->emulatorPrecisionMax.reset
    (new V(this->emulatorPrecisionSpace->zeroVector()));
  this->m_emulatorPrecisionShapeVec.reset
    (new V(this->emulatorPrecisionSpace->zeroVector()));
  this->m_emulatorPrecisionScaleVec.reset
    (new V(this->emulatorPrecisionSpace->zeroVector()));
  this->emulatorPrecisionMin->cwSet(0.3);
  this->emulatorPrecisionMax->cwSet(INFINITY);
  this->m_emulatorPrecisionShapeVec->cwSet(m_opts->m_emulatorPrecisionShape);
  this->m_emulatorPrecisionScaleVec->cwSet(m_opts->m_emulatorPrecisionScale);

  this->emulatorPrecisionDomain.reset
    (new BoxSubset<V, M>
     ("",
      *(this->emulatorPrecisionSpace),
      *(this->emulatorPrecisionMin),
      *(this->emulatorPrecisionMax)));

  this->m_emulatorPrecision.reset
    (new GammaVectorRV<V, M>
     ("",
      *(this->emulatorPrecisionDomain),
      *(this->m_emulatorPrecisionShapeVec),
      *(this->m_emulatorPrecisionScaleVec)));

  // Emulator correlation strength
  unsigned int dimScenario = (this->scenarioSpace()).dimLocal();
  unsigned int dimParameter = (this->parameterSpace()).dimLocal();
  this->emulatorCorrelationSpace.reset
    (new VectorSpace<V, M>
     (this->m_env,
      "",
      dimScenario + dimParameter,
      NULL));

  this->emulatorCorrelationMin.reset
    (new V(this->emulatorCorrelationSpace->zeroVector()));
  this->emulatorCorrelationMax.reset
    (new V(this->emulatorCorrelationSpace->zeroVector()));
  this->m_emulatorCorrelationStrengthAlphaVec.reset
    (new V(this->emulatorCorrelationSpace->zeroVector()));
  this->m_emulatorCorrelationStrengthBetaVec.reset
    (new V(this->emulatorCorrelationSpace->zeroVector()));
  this->emulatorCorrelationMin->cwSet(0);
  this->emulatorCorrelationMax->cwSet(1);
  this->m_emulatorCorrelationStrengthAlphaVec->cwSet(emulatorCorrelationStrengthAlpha);
  this->m_emulatorCorrelationStrengthBetaVec->cwSet(emulatorCorrelationStrengthBeta);

  this->emulatorCorrelationDomain.reset
    (new BoxSubset<V, M>
     ("",
      *(this->emulatorCorrelationSpace),
      *(this->emulatorCorrelationMin),
      *(this->emulatorCorrelationMax)));

  this->m_emulatorCorrelationStrength.reset
    (new BetaVectorRV<V, M>
     ("",
      *(this->emulatorCorrelationDomain),
      *(this->m_emulatorCorrelationStrengthAlphaVec),
      *(this->m_emulatorCorrelationStrengthBetaVec)));

  // Observation precision
  this->observationalPrecisionSpace.reset
    (new VectorSpace<V, M>
     (this->m_env,
      "",
      1,
      NULL));

  this->observationalPrecisionMin.reset
    (new V(this->observationalPrecisionSpace->zeroVector()));
  this->observationalPrecisionMax.reset
    (new V(this->observationalPrecisionSpace->zeroVector()));
  this->m_observationalPrecisionShapeVec.reset
    (new V(this->observationalPrecisionSpace->zeroVector()));
  this->m_observationalPrecisionScaleVec.reset
    (new V(this->observationalPrecisionSpace->zeroVector()));
  this->observationalPrecisionMin->cwSet(0.3);
  this->observationalPrecisionMax->cwSet(INFINITY);
  this->m_observationalPrecisionShapeVec->cwSet
    (this->m_opts->m_observationalPrecisionShape);
  this->m_observationalPrecisionScaleVec->cwSet
    (this->m_opts->m_observationalPrecisionScale);

  this->observationalPrecisionDomain.reset
    (new BoxSubset<V, M>
     ("",
      *(this->observationalPrecisionSpace),
      *(this->observationalPrecisionMin),
      *(this->observationalPrecisionMax)));

  this->m_observationalPrecision.reset
    (new GammaVectorRV<V, M>
     ("",
      *(this->observationalPrecisionDomain),
      *(this->m_observationalPrecisionShapeVec),
      *(this->m_observationalPrecisionScaleVec)));

  // Discrepancy precision
  this->discrepancyPrecisionMin.reset
    (new V(this->numDiscrepancyGroupsDSpace->zeroVector()));
  this->discrepancyPrecisionMax.reset
    (new V(this->numDiscrepancyGroupsDSpace->zeroVector()));
  this->m_discrepancyPrecisionShapeVec.reset
    (new V(this->numDiscrepancyGroupsDSpace->zeroVector()));
  this->m_discrepancyPrecisionScaleVec.reset
    (new V(this->numDiscrepancyGroupsDSpace->zeroVector()));
  this->discrepancyPrecisionMin->cwSet(0);
  this->discrepancyPrecisionMax->cwSet(INFINITY);
  this->m_discrepancyPrecisionShapeVec->cwSet(discrepancyPrecisionShape);
  this->m_discrepancyPrecisionScaleVec->cwSet(discrepancyPrecisionScale);

  this->discrepancyPrecisionDomain.reset
    (new BoxSubset<V, M>
     ("",
      *(this->numDiscrepancyGroupsDSpace),
      *(this->discrepancyPrecisionMin),
      *(this->discrepancyPrecisionMax)));

  this->m_discrepancyPrecision.reset
    (new GammaVectorRV<V, M>
     ("",
      *(this->discrepancyPrecisionDomain),
      *(this->m_discrepancyPrecisionShapeVec),
      *(this->m_discrepancyPrecisionScaleVec)));

  // Discrepancy correlation strength
  this->discrepancyCorrelationSpace.reset
    (new VectorSpace<V, M>
     (this->m_env,
      "",
      num_discrepancy_groups*dimScenario,
      NULL));

  this->discrepancyCorrelationMin.reset
    (new V(this->discrepancyCorrelationSpace->zeroVector()));
  this->discrepancyCorrelationMax.reset
    (new V(this->discrepancyCorrelationSpace->zeroVector()));
  this->m_discrepancyCorrelationStrengthAlphaVec.reset
    (new V(this->discrepancyCorrelationSpace->zeroVector()));
  this->m_discrepancyCorrelationStrengthBetaVec.reset
    (new V(this->discrepancyCorrelationSpace->zeroVector()));
  this->discrepancyCorrelationMin->cwSet(0);
  this->discrepancyCorrelationMax->cwSet(1);
  this->m_discrepancyCorrelationStrengthAlphaVec->cwSet(discrepancyCorrelationStrengthAlpha);
  this->m_discrepancyCorrelationStrengthBetaVec->cwSet(discrepancyCorrelationStrengthBeta);

  this->discrepancyCorrelationDomain.reset
    (new BoxSubset<V, M>
     ("",
      *(this->discrepancyCorrelationSpace),
      *(this->discrepancyCorrelationMin),
      *(this->discrepancyCorrelationMax)));

  this->m_discrepancyCorrelationStrength.reset
    (new BetaVectorRV<V, M>
     ("",
      *(this->discrepancyCorrelationDomain),
      *(this->m_discrepancyCorrelationStrengthAlphaVec),
      *(this->m_discrepancyCorrelationStrengthBetaVec)));

  // Emulator data precision
  this->emulatorDataPrecisionMin.reset
    (new V(this->oneDSpace->zeroVector()));
  this->emulatorDataPrecisionMax.reset
    (new V(this->oneDSpace->zeroVector()));
  this->m_emulatorDataPrecisionShapeVec.reset
    (new V(this->oneDSpace->zeroVector()));
  this->m_emulatorDataPrecisionScaleVec.reset
    (new V(this->oneDSpace->zeroVector()));
  this->emulatorDataPrecisionMin->cwSet(60.0);
  this->emulatorDataPrecisionMax->cwSet(1e5);
  this->m_emulatorDataPrecisionShapeVec->cwSet(emulatorDataPrecisionShape);
  this->m_emulatorDataPrecisionScaleVec->cwSet(emulatorDataPrecisionScale);

  this->emulatorDataPrecisionDomain.reset
    (new BoxSubset<V, M>
     ("",
      *(this->oneDSpace),
      *(this->emulatorDataPrecisionMin),
      *(this->emulatorDataPrecisionMax)));

  this->m_emulatorDataPrecision.reset
    (new GammaVectorRV<V, M>
     ("",
      *(this->emulatorDataPrecisionDomain),
      *(this->m_emulatorDataPrecisionShapeVec),
      *(this->m_emulatorDataPrecisionScaleVec)));

  // Now form full prior
  const unsigned int dimHyper =
    1 +
    (this->num_svd_terms < num_nonzero_eigenvalues) +
    this->m_opts->m_calibrateObservationalPrecision +
    num_svd_terms +
    dimParameter +
    dimScenario +
    num_discrepancy_groups +
    (num_discrepancy_groups*dimScenario);  // yum

  const unsigned int dimSum = dimParameter + dimHyper;

  this->hyperparamSpace.reset
    (new VectorSpace<V, M>
     (this->m_env,
      "",
      dimHyper,
      NULL));

  this->totalSpace.reset
    (new VectorSpace<V, M>
     (this->m_env,
      "",
      dimSum,
      NULL));
  this->hyperparamMins.reset(new V(this->hyperparamSpace->zeroVector()));
  this->hyperparamMaxs.reset(new V(this->hyperparamSpace->zeroVector()));

  this->hyperparamMins->cwSet(0.0);
  this->hyperparamMaxs->cwSet(1.0);

  for (unsigned int basis = 0; basis != num_svd_terms; ++basis) {
    // Min emulator precision
    (*(this->hyperparamMins))[(num_svd_terms<num_nonzero_eigenvalues)+basis] = 0.3;
    // Max emulator precision
    (*(this->hyperparamMaxs))[(num_svd_terms<num_nonzero_eigenvalues)+basis] = INFINITY;
  }

  // Starting index of the discrepancy precisions in the hyperparameter vector
  unsigned int discrepancyPrecisionIdx = (num_svd_terms<num_nonzero_eigenvalues) +
                                         num_svd_terms +
                                         dimScenario +
                                         dimParameter;

  for (unsigned int disc_grp = 0; disc_grp < num_discrepancy_groups; disc_grp++) {
    // Min discrepancy precision
    (*(this->hyperparamMins))[discrepancyPrecisionIdx+disc_grp] = 0;
    // Max discrepancy precision
    (*(this->hyperparamMaxs))[discrepancyPrecisionIdx+disc_grp] = INFINITY;
  }

  const int emulator_data_precision_index =
    dimHyper - 1 - this->m_opts->m_calibrateObservationalPrecision;
  (*(this->hyperparamMins))[emulator_data_precision_index] = 60.0;  // Min emulator data precision
  (*(this->hyperparamMaxs))[emulator_data_precision_index] = 1e5;   // Max emulator data precision

  if (this->m_opts->m_calibrateObservationalPrecision) {
    (*(this->hyperparamMins))[dimHyper-1] = 0.3;      // Min observation error precision
    (*(this->hyperparamMaxs))[dimHyper-1] = INFINITY; // Max observation error precision
  }

  this->hyperparamDomain.reset
    (new BoxSubset<V, M>
     ("",
      *(this->hyperparamSpace),
      *(this->hyperparamMins),
      *(this->hyperparamMaxs)));

  this->totalDomain.reset
    (new ConcatenationSubset<V, M>
     ("",
      *(this->totalSpace),
      this->m_parameterPrior.imageSet(),
      *(this->hyperparamDomain)));

  this->priors.push_back(&(this->m_parameterPrior));

  if (this->num_svd_terms < num_nonzero_eigenvalues)
    this->priors.push_back(this->m_truncationErrorPrecision.get());

  this->priors.push_back(this->m_emulatorPrecision.get());
  this->priors.push_back(this->m_emulatorCorrelationStrength.get());
  this->priors.push_back(this->m_discrepancyPrecision.get());
  this->priors.push_back(this->m_discrepancyCorrelationStrength.get());
  this->priors.push_back(this->m_emulatorDataPrecision.get());
  if (this->m_opts->m_calibrateObservationalPrecision)
    this->priors.push_back(this->m_observationalPrecision.get());

  // Finally
  this->m_totalPrior.reset
    (new ConcatenatedVectorRV<V, M>
     ("",
      this->priors,
      *(this->totalDomain)));
}

}  // End namespace QUESO

template class QUESO::GPMSAFactory<QUESO::GslVector, QUESO::GslMatrix>;
