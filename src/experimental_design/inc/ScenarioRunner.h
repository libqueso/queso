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

#ifndef QUESO_SCENARIO_RUNNER
#define QUESO_SCENARIO_RUNNER

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorRV.h>
#include <queso/VectorSet.h>
#include <queso/ExperimentalLikelihoodWrapper.h>
#include <queso/ExperimentMetricBase.h>

namespace QUESO
{
  //! Helper class that will run a statistical inverse problem for a 
  //! given set of experimental scenario parameters
  /*!
   *  Intended for use as part of an experimental design loop (e.g. GridSearchExperimentalDesign),
      this class automates the running of a statistical inverse problem and providing
      a metric by which the given scenario can be judged.
   */
  template<class V = GslVector, class M = GslMatrix>
  class ScenarioRunner
  {
  public:

    ScenarioRunner( std::shared_ptr<BaseVectorRV<V,M>> & prior,
                    std::shared_ptr<ExperimentalLikelihoodWrapper<V,M>> & wrapper,
                    std::shared_ptr<ExperimentMetricBase<V,M>> & metric)
    : m_prior(prior),
      m_wrapper(wrapper),
      m_experiment_metric(metric) {}

    //! Passes the scenario parameters to the likelihood function,
    //! creates and runs a SIP object, and sets the m_exp_metric_value
    //! based on the provided experiment metric
    void runExperiment(std::vector<double> & scenario_params);

    double getExperimentMetricValue()
    {
      return m_exp_metric_value;
    }

  protected:
    std::shared_ptr<BaseVectorRV<V,M>> m_prior;
    std::shared_ptr<ExperimentalLikelihoodWrapper<V,M>> m_wrapper;
    std::shared_ptr<ExperimentMetricBase<V,M>> m_experiment_metric;

    // The "value" of the given scenario based on the provided experiment metric
    double m_exp_metric_value;
  };

}

#endif //QUESO_SCENARIO_RUNNER

