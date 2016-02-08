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

#include <queso/GPMSA.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

template <class V, class M>
GPMSAEmulator<V, M>::GPMSAEmulator(
    const VectorSet<V, M> & domain,
    const VectorSpace<V, M> & m_scenarioSpace,
    const VectorSpace<V, M> & m_parameterSpace,
    const VectorSpace<V, M> & m_simulationOutputSpace,
    const VectorSpace<V, M> & m_experimentOutputSpace,
    const unsigned int m_numSimulations,
    const unsigned int m_numExperiments,
    const std::vector<V *> & m_simulationScenarios,
    const std::vector<V *> & m_simulationParameters,
    const std::vector<V *> & m_simulationOutputs,
    const std::vector<V *> & m_experimentScenarios,
    const std::vector<V *> & m_experimentOutputs,
    const std::vector<V>   & m_discrepancyBases,
    const std::vector<M>   & m_observationErrorMatrices,
    const M & m_experimentErrors,
    const ConcatenatedVectorRV<V, M> & m_totalPrior)
  :
  BaseScalarFunction<V, M>("", m_totalPrior.imageSet()),
  m_scenarioSpace(m_scenarioSpace),
  m_parameterSpace(m_parameterSpace),
  m_simulationOutputSpace(m_simulationOutputSpace),
  m_experimentOutputSpace(m_experimentOutputSpace),
  m_numSimulations(m_numSimulations),
  m_numExperiments(m_numExperiments),
  m_simulationScenarios(m_simulationScenarios),
  m_simulationParameters(m_simulationParameters),
  m_simulationOutputs(m_simulationOutputs),
  m_experimentScenarios(m_experimentScenarios),
  m_experimentOutputs(m_experimentOutputs),
  m_discrepancyBases(m_discrepancyBases),
  m_observationErrorMatrices(m_observationErrorMatrices),
  m_experimentErrors(m_experimentErrors),
  m_totalPrior(m_totalPrior)
{
  queso_assert_greater(m_numSimulations, 0);

  const MpiComm & comm = m_simulationOutputs[0]->map().Comm();

  queso_assert_equal_to(comm.NumProc(), 1);

  Map serial_map(m_numSimulations, 0, comm);

  unsigned int numOutputs = this->m_experimentOutputSpace.dimLocal();

  const BaseEnvironment &env = m_simulationOutputs[0]->env();

  M simulation_matrix(env, serial_map, numOutputs);

  for (unsigned int i=0; i != m_numSimulations; ++i)
    for (unsigned int j=0; j != numOutputs; ++j)
      simulation_matrix(i,j) = 
        (*m_simulationOutputs[i])[j];

  // GSL only finds left singular vectors if n_rows>=n_columns, so we need to
  // calculate them indirectly from the eigenvalues of M^T*M

  M S_trans(simulation_matrix.transpose());

  M SM_squared(S_trans*simulation_matrix);

  M SM_singularVectors(env, SM_squared.map(), numOutputs);
  V SM_singularValues(env, SM_squared.map());

  SM_squared.eigen(SM_singularValues, &SM_singularVectors);

  const unsigned int MAX_SVD_TERMS =
    std::min(m_numSimulations,(unsigned int)(5));
  num_svd_terms = std::min(MAX_SVD_TERMS, numOutputs);

  // Copy only those vectors we want into K_eta
  m_TruncatedSVD_simulationOutputs.reset
    (new M(m_simulationOutputs[0]->env(),
           m_simulationOutputs[0]->map(),
           num_svd_terms));

  for (unsigned int i=0; i != numOutputs; ++i)
    for (unsigned int k = 0; k != num_svd_terms; ++k)
      (*m_TruncatedSVD_simulationOutputs)(i,k) = SM_singularVectors(i,k);

  Map copied_map(numOutputs * m_numSimulations, 0,
                 m_simulationOutputs[0]->map().Comm());

  K.reset
    (new M(m_simulationOutputs[0]->env(), copied_map,
           m_numSimulations));
  for (unsigned int k=0; k != num_svd_terms; ++k)
    for (unsigned int i1=0; i1 != m_numSimulations; ++i1)
      for (unsigned int i2=0; i2 != numOutputs; ++i2)
        {
          const unsigned int i = i1 * numOutputs + i2;
          const unsigned int j = k * m_numSimulations + i1;
          (*K)(i,j) = SM_singularVectors(i2,k);
        }

  KKT.reset
    (new M(m_simulationOutputs[0]->env(), copied_map,
           numOutputs * m_numSimulations));

  K->multiply(K->transpose(), *KKT);

  Map outputs_map(numOutputs, 0, comm);

  for (unsigned int i = 0; i != m_numExperiments; ++i)
    {
      M D_i(m_simulationOutputs[0]->env(), outputs_map,
            (unsigned int)(m_discrepancyBases.size()));

      for (unsigned int j=0; j != numOutputs; ++j)
        for (unsigned int k=0; k != m_discrepancyBases.size(); ++k)
          D_i(j,k) = m_discrepancyBases[k][j];

      m_discrepancyMatrices.push_back(D_i);
    }

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
  const unsigned int Brows = m_numExperiments * numOutputs;
  const unsigned int Bcols =
    m_numExperiments * (num_discrepancy_bases + num_svd_terms);

  const Map B_row_map(Brows, 0, comm);

  m_BMatrix.reset
    (new M(m_simulationOutputs[0]->env(), B_row_map, Bcols));

  const unsigned int Wyrows = m_numExperiments * numOutputs;

  const Map Wy_row_map(Wyrows, 0, comm);

  m_observationErrorMatrix.reset
    (new M(m_simulationOutputs[0]->env(), Wy_row_map, Wyrows));

  M& B = *m_BMatrix;

  M& Wy = *m_observationErrorMatrix;

  for (unsigned int ex = 0; ex != m_numExperiments; ++ex)
    {
      const M & D_i = m_discrepancyMatrices[ex];

      for (unsigned int outi = 0; outi != numOutputs; ++outi)
        {
          unsigned int i = ex*numOutputs+outi;
          for (unsigned int outj = 0; outj != num_discrepancy_bases; ++outj)
            {
              unsigned int j = ex + m_numExperiments * outj;

              B(i,j) = D_i(outi,outj);
            }

          for (unsigned int outj = 0; outj != num_svd_terms; ++outj)
            {
              unsigned int j = ex +
                m_numExperiments * (num_discrepancy_bases + outj);

              B(i,j) = D_i(outi,outj);
            }

          // For the multivariate case, the bases K_eta computed from
          // simulator outputs are the same as the bases K_i which
          // apply to quantities of interest, because simulator
          // outputs are QoIs alone.
          //
          // FIXME - we need to interpolate K_i in the functional case.

          for (unsigned int outj = 0; outj != numOutputs; ++outj)
            {
              // No fancy perturbation here
              unsigned int j = ex*numOutputs+outj;

              Wy(i,j) = m_observationErrorMatrices[ex](outi,outj);
            }
        }
    }

  M BT_Wy_B = B.transpose() * Wy * B;
}

template <class V, class M>
GPMSAEmulator<V, M>::~GPMSAEmulator()
{
  // heap allocations get deleted via ScopedPtr
}

template <class V, class M>
double
GPMSAEmulator<V, M>::lnValue(const V & domainVector,
                                       const V * domainDirection,
                                       V * gradVector,
                                       M * hessianMatrix,
                                       V * hessianEffect) const
{
  // Components of domainVector:
  // theta(1)                     // = "theta", "t" in Higdon et. al. 2008
  // theta(2)
  // ...
  // theta(dimParameter)
  // emulator_mean                // = "mu"
  // emulator_precision(1)        // = "lambda_eta" in scalar case,
  // ...                          //   "lambda_{wi}" in vector
  // emulator_precision(num_svd_terms)
  // emulator_corr_strength(1)    // = "rho_{eta k}"
  // ...                          // dimScenario = "p_x", dimParameter = "p_t"
  // emulator_corr_strength(dimScenario + dimParameter)
  // discrepancy_precision(1)     // = "lambda_delta" in scalar case,
  //                              //   "lambda_{vi}" in vector
  // ...
  // discrepancy_precision(F)     // FIXME: F = 1, |G_1| = p_delta, for now.
  // discrepancy_corr_strength(1) // = "rho_delta"
  // ...
  // discrepancy_corr_strength(dimScenario)
  // emulator_data_precision(1)   // = "small white noise", "small ridge"

  // Other variables:
  // m_numSimulations             // = "m"
  // m_numExperiments             // = "n"
  // m_simulationScenarios        // = "eta"
  // m_experimentScenarios        // = "y"
  // m_experimentErrors           // = "Sigma_y"
  // dimScenario                  // = "p_x", "n_eta"
  //                              // "n_y" := sum(n_y_i)
  //                              //      (== n*n_eta for us for now)
  // num_svd_terms                // = "p_eta"
  // num_discrepancy_bases        // = "p_delta"
  // m_TruncatedSVD_simulationOutputs  // = "K_eta"
  // covMatrix                    // = "Sigma_D" in scalar case,
  //                                   "Sigma_zhat" in vector
  // m_discrepancyMatrices        // = "D_i"
  // m_observationErrorMatrices   // = "W_i"
  // m_observationErrorMatrix     // = "W_y"
  // m_BMatrix                    // = "B"
  // num_discrepancy_groups       // = "F"

  // Construct covariance matrix
  const unsigned int totalRuns = this->m_numExperiments + this->m_numSimulations;
  const unsigned int numOutputs = this->m_experimentOutputSpace.dimLocal();
  const unsigned int totalOutputs = totalRuns * numOutputs;
  const unsigned int num_discrepancy_bases = m_discrepancyBases.size();
  const unsigned int residualSize = (numOutputs == 1) ?  totalOutputs :
    totalOutputs * num_svd_terms + m_numExperiments * num_discrepancy_bases;

  double prodScenario = 1.0;
  double prodParameter = 1.0;
  double prodDiscrepancy = 1.0;
  unsigned int dimScenario = (this->m_scenarioSpace).dimLocal();
  unsigned int dimParameter = (this->m_parameterSpace).dimLocal();

  // This is cumbersome.  All I want is a matrix.
  VectorSpace<V, M> gpSpace(this->m_scenarioSpace.env(), "",
                            residualSize, NULL);
  V residual(gpSpace.zeroVector());
  M covMatrix(residual);

  V domainVectorParameter(*(this->m_simulationParameters[0]));
  for (unsigned int k = 0; k < dimParameter; k++) {
    queso_assert (!isnan(domainVector[k]));
    domainVectorParameter[k] = domainVector[k];
  }

  // This for loop is a disaster and could do with a *lot* of optimisation
  for (unsigned int i = 0; i < totalRuns; i++) {

     const V* scenario1;
     const V* parameter1;

    // Decide whether to do experiment part of the covariance matrix
    // Get i-th simulation-parameter pair
    if (i < this->m_numExperiments) {
      // Experiment scenario (known)
      scenario1 = (this->m_experimentScenarios)[i];

      // Experiment parameter (unknown)
      parameter1 = &domainVectorParameter;
    }
    else {
      scenario1 =
        (this->m_simulationScenarios)[i-this->m_numExperiments];
      parameter1 =
        (this->m_simulationParameters)[i-this->m_numExperiments];
    }

    for (unsigned int j = 0; j < totalRuns; j++) {

       const V* scenario2;
       const V* parameter2;

      if (j < this->m_numExperiments) {
        scenario2 = (this->m_experimentScenarios)[j];
        parameter2 = &domainVectorParameter;
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
        dimParameter + 1 + num_svd_terms;
      for (unsigned int k = 0; k < dimScenario; k++) {
        const double & emulator_corr_strength =
          domainVector[emulatorCorrStrStart+k];
        prodScenario *= std::pow(emulator_corr_strength,
                                 4.0 * ((*scenario1)[k] - (*scenario2)[k]) *
                                       ((*scenario1)[k] - (*scenario2)[k]));
      }

      queso_assert (!isnan(prodScenario));

      // = second term in (1)
      prodParameter = 1.0;
      for (unsigned int k = 0; k < dimParameter; k++) {
        queso_assert (!isnan(domainVector[emulatorCorrStrStart+dimScenario+k]));
        queso_assert (!isnan((*parameter1)[k]));
        queso_assert (!isnan((*parameter2)[k]));
        const double & emulator_corr_strength =
          domainVector[emulatorCorrStrStart+dimScenario+k];
        prodParameter *= std::pow(
            emulator_corr_strength,
            4.0 * ((*parameter1)[k] - (*parameter2)[k]) *
                  ((*parameter1)[k] - (*parameter2)[k]));
      }

      queso_assert (!isnan(prodParameter));

      // Sigma_eta in scalar case,
      // Sigma_u in vector case
      for (unsigned int basis = 0; basis != num_svd_terms; ++basis)
        {
          // Offset in vector case
          const unsigned int offset = (numOutputs == 1) ?
            0 : m_numExperiments * num_discrepancy_bases;

          // coefficient in (1)
          const double emulator_precision =
            domainVector[dimParameter+basis+1];
          queso_assert_greater(emulator_precision, 0);

          covMatrix(offset+basis*m_numSimulations+i,
                    offset+basis*m_numSimulations+j) =
            prodScenario * prodParameter / emulator_precision;
        }

      // If we're in the experiment cross correlation part, need extra
      // foo: Sigma_delta/Sigma_v and Sigma_y
      if (i < this->m_numExperiments && j < this->m_numExperiments) {
        V* cross_scenario1 = (this->m_simulationScenarios)[i];
        V* cross_scenario2 = (this->m_simulationScenarios)[j];
        prodDiscrepancy = 1.0;
        unsigned int discrepancyCorrStrStart = dimParameter +
                                               num_svd_terms +
                                               dimParameter +
                                               dimScenario + 2;
        for (unsigned int k = 0; k < dimScenario; k++) {
          const double & discrepancy_corr_strength = 
            domainVector[discrepancyCorrStrStart+k];
          prodDiscrepancy *=
            std::pow(discrepancy_corr_strength, 4.0 *
                     ((*cross_scenario1)[k] - (*cross_scenario2)[k]) *
                     ((*cross_scenario1)[k] - (*cross_scenario2)[k]));
        }

        const double discrepancy_precision =
          domainVector[discrepancyCorrStrStart-1];
        queso_assert_greater(discrepancy_precision, 0);
        queso_assert (!isnan(prodDiscrepancy));

            // Sigma_delta term from below (3) in univariate case
            // Sigma_v term from p. 576 in multivariate case
        const double R_v = prodDiscrepancy / discrepancy_precision;
        for (unsigned int disc = 0; disc != num_discrepancy_bases;
             ++disc)
          covMatrix(disc*m_numExperiments+i,
                    disc*m_numExperiments+j) += R_v;

        // Experimental error comes in via K in the multivariate
        // case
        if (numOutputs == 1)
          {
            // Sigma_y term from below (3)
            const double experimentalError =
              (this->m_experimentErrors)(i,j);

            queso_assert_greater_equal (experimentalError, 0);

            covMatrix(i,j) += experimentalError;
          }
      }
    }

    // Add small white noise component to diagonal to make stuff +ve def
    // = "small ridge"
    unsigned int dimSum = 3 +
                          num_svd_terms +
                          dimParameter +
                          dimParameter +
                          dimScenario +
                          dimScenario;  // yum

    queso_assert_greater(domainVector[dimSum-1], 0);
    double nugget = 1.0 / domainVector[dimSum-1];

    for (unsigned int disc = 0; disc != num_discrepancy_bases;
         ++disc)
      covMatrix(disc*m_numExperiments+i,
                disc*m_numExperiments+j) += nugget;
  }

  // Form residual = D - mean // = D - mu*1 in (3)
  // We don't subtract off mean here because we expect normalized data
  for (unsigned int i = 0; i < this->m_numExperiments; i++) {
    for (unsigned int k = 0; k != numOutputs; ++k)
      residual[i*numOutputs+k] =
        (*((this->m_experimentOutputs)[i]))[k];
  }
  for (unsigned int i = 0; i < this->m_numSimulations; i++) {
    for (unsigned int k = 0; k != numOutputs; ++k)
      residual[(i+this->m_numExperiments)*numOutputs+k] =
        (*((this->m_simulationOutputs)[i]))[k];
  }

  // Solve covMatrix * sol = residual
  // = Sigma_D^-1 * (D - mu 1) from (3)
  V sol(covMatrix.invertMultiply(residual));

  // Premultiply by residual^T as in (3)
  double minus_2_log_lhd = 0.0;
  // There's no dot product function in GslVector.
  for (unsigned int i = 0; i < totalOutputs; i++) {
    minus_2_log_lhd += sol[i] * residual[i];
  }

if (isnan(minus_2_log_lhd))
  for (unsigned int i = 0; i < totalOutputs; i++) {
    if (isnan(sol[i]))
      std::cout << "NaN sol[" << i << ']' << std::endl;
    if (isnan(residual[i]))
      std::cout << "NaN residual[" << i << ']' << std::endl;

    std::cout << "Covariance Matrix:" << std::endl;
    covMatrix.print(std::cout);
  }

// std::cout << "minus_2_log_lhd = " << minus_2_log_lhd << std::endl;

  queso_assert_greater(minus_2_log_lhd, 0);

  // FIXME: missing ln(det(Sigma_D)) term!
  minus_2_log_lhd += std::log(covMatrix.determinant());

  // Multiply by -1/2 coefficient from (3)
  return -0.5 * minus_2_log_lhd;
}

template <class V, class M>
double
GPMSAEmulator<V, M>::actualValue(const V & domainVector,
                                           const V * domainDirection,
                                           V * gradVector,
                                           M * hessianMatrix,
                                           V * hessianEffect) const
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
    const VectorSpace<V, M> & experimentOutputSpace,
    unsigned int numSimulations,
    unsigned int numExperiments)
  :
    m_env(env),
    m_parameterPrior(parameterPrior),
    m_scenarioSpace(scenarioSpace),
    m_parameterSpace(parameterSpace),
    m_simulationOutputSpace(simulationOutputSpace),
    m_experimentOutputSpace(experimentOutputSpace),
    m_numSimulations(numSimulations),
    m_numExperiments(numExperiments),
    m_simulationScenarios(numSimulations, (V *)NULL),
    m_simulationParameters(numSimulations, (V *)NULL),
    m_simulationOutputs(numSimulations, (V *)NULL),
    m_experimentScenarios(numExperiments, (V *)NULL),
    m_experimentOutputs(numExperiments, (V *)NULL),
    m_numSimulationAdds(0),
    m_numExperimentAdds(0),
    priors(7, (const BaseVectorRV<V, M> *)NULL),  // Needed for gcc 4.3.2
    m_constructedGP(false)
{
  // We should have the same number of outputs from both simulations
  // and experiments
  queso_assert_equal_to(simulationOutputSpace.dimGlobal(),
                        experimentOutputSpace.dimGlobal());

  {
    const unsigned int numOutputs =
      this->m_experimentOutputSpace.dimLocal();
    const Map & output_map = experimentOutputSpace.map();

    // Set up the default discrepancy basis:
    V all_ones_basis(env, output_map);
    for (unsigned int i=0; i != numOutputs; ++i)
      all_ones_basis[i] = 1;
    m_discrepancyBases.push_back(all_ones_basis);

    // Set up the default observation error covariance matrix:

    M identity_matrix(env, output_map, 1.0);

    for (unsigned int i = 0; i != numExperiments; ++i)
      m_observationErrorMatrices.push_back(identity_matrix);
  }

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

  this->setUpHyperpriors();

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
const VectorSpace<V, M> &
GPMSAFactory<V, M>::experimentOutputSpace() const
{
  return this->m_experimentOutputSpace;
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
const std::vector<V *> &
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
const std::vector<V *> &
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
const std::vector<V *> &
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
const std::vector<V *> &
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
const std::vector<V *> &
GPMSAFactory<V, M>::experimentOutputs() const
{
  return this->m_experimentOutputs;
}

template <class V, class M>
const M &
GPMSAFactory<V, M>::experimentErrors() const
{
  return *(this->m_experimentErrors);
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
GPMSAFactory<V, M>::addSimulation(V & simulationScenario,
                                            V & simulationParameter,
                                            V & simulationOutput)
{
  queso_require_less_msg(this->m_numSimulationAdds, this->m_numSimulations, "too many simulation adds...");

  this->m_simulationScenarios[this->m_numSimulationAdds] = &simulationScenario;
  this->m_simulationParameters[this->m_numSimulationAdds] = &simulationParameter;
  this->m_simulationOutputs[this->m_numSimulationAdds] = &simulationOutput;
  this->m_numSimulationAdds++;

  if ((this->m_numSimulationAdds == this->m_numSimulations) &&
      (this->m_numExperimentAdds == this->m_numExperiments) &&
      (this->m_constructedGP == false)) {
    this->m_constructedGP = true;
    this->gpmsaEmulator.reset
      (new GPMSAEmulator<V, M>(
        this->prior().imageSet(),
        this->m_scenarioSpace,
        this->m_parameterSpace,
        this->m_simulationOutputSpace,
        this->m_experimentOutputSpace,
        this->m_numSimulations,
        this->m_numExperiments,
        this->m_simulationScenarios,
        this->m_simulationParameters,
        this->m_simulationOutputs,
        this->m_experimentScenarios,
        this->m_experimentOutputs,
        this->m_discrepancyBases,
        this->m_observationErrorMatrices,
        *(this->m_experimentErrors),
        *(this->m_totalPrior)));
  }
}

template <class V, class M>
void
GPMSAFactory<V, M>::addSimulations(
    const std::vector<V *> & simulationScenarios,
    const std::vector<V *> & simulationParameters,
    const std::vector<V *> & simulationOutputs)
{
  for (unsigned int i = 0; i < this->m_numSimulations; i++) {
    this->addSimulation(*(simulationScenarios[i]), *(simulationParameters[i]),
        *(simulationOutputs[i]));
  }
}

template <class V, class M>
void
GPMSAFactory<V, M>::addExperiments(
    const std::vector<V *> & experimentScenarios,
    const std::vector<V *> & experimentOutputs,
    const M * experimentErrors)
{
  queso_require_less_equal_msg(experimentScenarios.size(), this->m_numExperiments, "too many experiments...");

  for (unsigned int i = 0; i < this->m_experimentScenarios.size(); i++) {
    this->m_experimentScenarios[i] = experimentScenarios[i];
    this->m_experimentOutputs[i] = experimentOutputs[i];
  }
  this->m_experimentErrors = experimentErrors;
  this->m_numExperimentAdds += experimentScenarios.size();

  if ((this->m_numSimulationAdds == this->m_numSimulations) &&
      (this->m_numExperimentAdds == this->m_numExperiments) &&
      (this->m_constructedGP == false)) {
    this->m_constructedGP = true;
    this->gpmsaEmulator.reset
      (new GPMSAEmulator<V, M>(
        this->prior().imageSet(),
        this->m_scenarioSpace,
        this->m_parameterSpace,
        this->m_simulationOutputSpace,
        this->m_experimentOutputSpace,
        this->m_numSimulations,
        this->m_numExperiments,
        this->m_simulationScenarios,
        this->m_simulationParameters,
        this->m_simulationOutputs,
        this->m_experimentScenarios,
        this->m_experimentOutputs,
        this->m_discrepancyBases,
        this->m_observationErrorMatrices,
        *(this->m_experimentErrors),
        *(this->m_totalPrior)));
  }
}


template <class V, class M>
void
GPMSAFactory<V, M>::setDiscrepancyBases(
    const std::vector<V *> & discrepancyBases)
{
  m_discrepancyBases.clear();

  for (unsigned int i = 0; i < discrepancyBases.size(); i++) {
    m_discrepancyBases.push_back(*(discrepancyBases[i]));
  }

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

  return m_observationErrorMatrices[simulationNumber];
}


template <class V, class M>
const M &
GPMSAFactory<V, M>::getObservationErrorCovariance
  (unsigned int simulationNumber) const
{
  queso_assert_less(simulationNumber, m_numSimulations);
  queso_assert_equal_to(m_observationErrorMatrices.size(), m_numSimulations);

  return m_observationErrorMatrices[simulationNumber];
}



template <class V, class M>
const ConcatenatedVectorRV<V, M> &
GPMSAFactory<V, M>::prior() const
{
  return *(this->m_totalPrior);
}

template <class V, class M>
void
GPMSAFactory<V, M>::print(std::ostream& os) const
{
  // Do nothing
}

// Private methods follow
template <class V, class M>
void
GPMSAFactory<V, M>::setUpHyperpriors()
{
  double emulatorPrecisionShape = this->m_opts->m_emulatorPrecisionShape;
  double emulatorPrecisionScale = this->m_opts->m_emulatorPrecisionScale;
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

  // Emulator mean
  this->emulatorMeanMin.reset(new V(this->oneDSpace->zeroVector()));
  this->emulatorMeanMax.reset(new V(this->oneDSpace->zeroVector()));
  this->emulatorMeanMin->cwSet(-INFINITY);
  this->emulatorMeanMax->cwSet(INFINITY);

  this->emulatorMeanDomain.reset
    (new BoxSubset<V, M>
      ("",
       *(this->oneDSpace),
       *(this->emulatorMeanMin),
       *(this->emulatorMeanMax)));

  this->m_emulatorMean.reset
    (new UniformVectorRV<V, M>
     ("",
      *(this->emulatorMeanDomain)));

  // Emulator precision
  this->emulatorPrecisionMin.reset(new V(this->oneDSpace->zeroVector()));
  this->emulatorPrecisionMax.reset(new V(this->oneDSpace->zeroVector()));
  this->m_emulatorPrecisionShapeVec.reset(new V(this->oneDSpace->zeroVector()));
  this->m_emulatorPrecisionScaleVec.reset(new V(this->oneDSpace->zeroVector()));
  this->emulatorPrecisionMin->cwSet(0.3);
  this->emulatorPrecisionMax->cwSet(INFINITY);
  this->m_emulatorPrecisionShapeVec->cwSet(emulatorPrecisionShape);
  this->m_emulatorPrecisionScaleVec->cwSet(emulatorPrecisionScale);

  this->emulatorPrecisionDomain.reset
    (new BoxSubset<V, M>
     ("",
      *(this->oneDSpace),
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

  // Discrepancy precision
  this->discrepancyPrecisionMin.reset
    (new V(this->oneDSpace->zeroVector()));
  this->discrepancyPrecisionMax.reset 
    (new V(this->oneDSpace->zeroVector()));
  this->m_discrepancyPrecisionShapeVec.reset
    (new V(this->oneDSpace->zeroVector()));
  this->m_discrepancyPrecisionScaleVec.reset
    (new V(this->oneDSpace->zeroVector()));
  this->discrepancyPrecisionMin->cwSet(0);
  this->discrepancyPrecisionMax->cwSet(INFINITY);
  this->m_discrepancyPrecisionShapeVec->cwSet(discrepancyPrecisionShape);
  this->m_discrepancyPrecisionScaleVec->cwSet(discrepancyPrecisionScale);

  this->discrepancyPrecisionDomain.reset
    (new BoxSubset<V, M>
     ("",
      *(this->oneDSpace),
      *(this->discrepancyPrecisionMin),
      *(this->emulatorPrecisionMax)));

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
      dimScenario,
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
  unsigned int dimSum = 4 +
                        dimParameter +
                        dimParameter +
                        dimScenario +
                        dimScenario;  // yum

  this->totalSpace.reset
    (new VectorSpace<V, M>
     (this->m_env,
      "",
      dimSum,
      NULL));
  this->totalMins.reset(new V(this->totalSpace->zeroVector()));
  this->totalMaxs.reset(new V(this->totalSpace->zeroVector()));

  // Hackety hack McHackington.  There's no better way to do this unfortunately
  this->totalMins->cwSet(0);
  this->totalMaxs->cwSet(1);

  (*(this->totalMins))[dimParameter] = -INFINITY;  // Min mean
  (*(this->totalMaxs))[dimParameter] = INFINITY;  // Max mean
  (*(this->totalMins))[dimParameter+1] = 0.3;  // Min emulator precision
  (*(this->totalMaxs))[dimParameter+1] = INFINITY;  // Max emulator precision
  (*(this->totalMins))[dimSum-1] = 60.0;  // Min emulator data precision
  (*(this->totalMaxs))[dimSum-1] = 1e5;  // Max emulator data precision

  // Min discrepancy precision
  (*(this->totalMins))[dimScenario+dimParameter+dimParameter+2] = 0;
  // Max discrepancy precision
  (*(this->totalMaxs))[dimScenario+dimParameter+dimParameter+2] = INFINITY;

  this->totalDomain.reset
    (new BoxSubset<V, M>
     ("",
      *(this->totalSpace),
      *(this->totalMins),
      *(this->totalMaxs)));

  this->priors[0] = &(this->m_parameterPrior);
  this->priors[1] = this->m_emulatorMean.get();
  this->priors[2] = this->m_emulatorPrecision.get();
  this->priors[3] = this->m_emulatorCorrelationStrength.get();
  this->priors[4] = this->m_discrepancyPrecision.get();
  this->priors[5] = this->m_discrepancyCorrelationStrength.get();
  this->priors[6] = this->m_emulatorDataPrecision.get();

  // Finally
  this->m_totalPrior.reset
    (new ConcatenatedVectorRV<V, M>
     ("",
      this->priors,
      *(this->totalDomain)));
}

}  // End namespace QUESO

template class QUESO::GPMSAFactory<QUESO::GslVector, QUESO::GslMatrix>;
