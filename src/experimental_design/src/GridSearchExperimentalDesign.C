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

// This class
#include <queso/GridSearchExperimentalDesign.h>

// QUESO
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/BoxSubset.h>
#include <queso/ScenarioRunner.h>
#include <queso/MultiDimensionalIndexing.h>

namespace QUESO
{
  template<class V, class M>
  GridSearchExperimentalDesign<V,M>::GridSearchExperimentalDesign(const BoxSubset<V,M> & scenario_domain,
                                                              std::vector<unsigned int> & n_points,
                                                              std::shared_ptr<ScenarioRunner<V,M>> & runner)
  : m_scenario_domain(scenario_domain),
    m_n_points(n_points),
    m_scenario_runner(runner)
  {
    unsigned int n_param = m_scenario_domain.vectorSpace().dimGlobal();

    queso_assert_equal_to(n_param,n_points.size());

    unsigned int total_points = 1;
    for (unsigned int d = 0; d < n_param; ++d)
      total_points *= m_n_points[d];

    m_total_scenarios = total_points;
  }


  template<class V, class M>
  void GridSearchExperimentalDesign<V,M>::run(V & experimental_params, std::string & filename_prefix)
  {

    // don't need to store all the metric values,
    // just the largest one
    double max_value = -1.0; // dummy initial value
    unsigned int max_value_index = 0;

    for (unsigned int n = 0; n < m_total_scenarios; ++n)
      {
        std::vector<double> scenario(this->get_n_params());

        this->get_params_from_global_coord(n,scenario);

        if (this->m_scenario_domain.env().fullRank() == 0)
          {
            std::stringstream ss;
            ss <<"\n-------------------------------\n"
                <<"Running scenario " <<n <<std::endl
                <<"Parameters: ";
            for (unsigned int s = 0; s < (this->get_n_params()-1); ++s)
              ss <<scenario[s] <<",";

            ss  <<scenario[this->get_n_params()-1]
                <<"\n-------------------------------\n";

            std::cout <<ss.str();
          }

        std::string filename = "";
        if (filename_prefix != "")
          filename = "posterior_"+filename_prefix+"_"+std::to_string(n);

        m_scenario_runner->runExperiment(scenario,filename);

        double value = m_scenario_runner->getExperimentMetricValue();

        if (this->m_scenario_domain.env().fullRank() == 0)
          {
            std::stringstream ss1;
            ss1 <<"-------------------------------\n"
                <<"Scenario " <<n <<" metric: " <<value
                <<"\n-------------------------------\n";

            std::cout <<ss1.str();
          }

        // check if n==0 because we always want to set max_value on the first iteration
        if ( (value > max_value) || (n == 0) )
          {
            max_value = value;
            max_value_index = n;
          }

      }

    std::vector<double> params(this->get_n_params());

    this->get_params_from_global_coord(max_value_index,params);
    
    for (unsigned int p = 0; p < this->get_n_params(); ++p)
      experimental_params[p] = params[p];

  }


  template<class V, class M>
  void GridSearchExperimentalDesign<V,M>::get_params_from_global_coord(unsigned int n, std::vector<double> & param_values)
  {
    // copied fromm InterpolationSurrogateBuilder::set_domain_vector()

    std::vector<unsigned int> indices(m_total_scenarios);
    MultiDimensionalIndexing::globalToCoord(n,m_n_points,indices);

    for(unsigned int d = 0; d < this->get_n_params(); ++d)
      param_values[d] = this->get_param_value(indices[d],d);

  }

}

template class QUESO::GridSearchExperimentalDesign<QUESO::GslVector,QUESO::GslMatrix>;

