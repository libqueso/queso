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

#ifndef QUESO_GRID_SEARCH_EXPERIMENTAL_DESIGN
#define QUESO_GRID_SEARCH_EXPERIMENTAL_DESIGN

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/BoxSubset.h>
#include <queso/ScenarioRunner.h>

namespace QUESO
{
  //! Experimental design class that will iterate over a discrete domain of
  //! experimental scenario parameters (e.g. temperature, pressure, mesh length)
  //! and chose the "best" experiment based on the provided metric
  template<class V = GslVector, class M = GslMatrix>
  class GridSearchExperimentalDesign
  {
  public:

    GridSearchExperimentalDesign( const BoxSubset<V,M> & scenario_domain,
                                  std::vector<unsigned int> & n_points,
                                  std::shared_ptr<ScenarioRunner<V,M>> & runner);

    //! Iterate over all possible experimental scenarios and calculate a metric for each.
    //!
    //! @return The scenario parameter values of the experiment with the highest metric value
    void run(V & experimental_params, std::string & filename_prefix);

    GridSearchExperimentalDesign() = delete;

  private:
    //! Experimental scenario parameter domain
    const BoxSubset<V,M> & m_scenario_domain;

    //! The number of points by which each scenario parameter should be discretized
    std::vector<unsigned int> & m_n_points;

    //! A ScenarioRunner to facilitate the solution and measurement of each scenario
    std::shared_ptr<ScenarioRunner<V,M>> m_scenario_runner;

    //! The total number of scenarios to evaluate
    unsigned int m_total_scenarios;

    //! Helper function to get scenarion parameter values for the global coordinate n
    void get_params_from_global_coord(unsigned int n, std::vector<double> & param_values);

    double get_param_value(unsigned int index, unsigned int scenario_param)
    {
      double min  = this->get_param_min(scenario_param);
      double step = this->get_param_step_size(scenario_param);
      
      return min + index*step;
    }

    double get_param_min(unsigned int scenario_param)
    {
      return this->m_scenario_domain.minValues()[scenario_param];
    }

    double get_param_max(unsigned int scenario_param)
    {
      return this->m_scenario_domain.maxValues()[scenario_param];
    }

    unsigned int get_n_params()
    {
      return this->m_n_points.size();
    }

    double get_param_step_size(unsigned int scenario_param)
    {
      double min = this->get_param_min(scenario_param);
      double max = this->get_param_max(scenario_param);
      unsigned int n_steps = m_n_points[scenario_param] - 1;
      
      return (max-min)/n_steps;
    }

  };
}

#endif //QUESO_GRID_SEARCH_EXPERIMENTAL_DESIGN

