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
#include <queso/ScenarioRunner.h>

// QUESO
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorRV.h>
#include <queso/VectorSet.h>
#include <queso/ExperimentalLikelihoodWrapper.h>
#include <queso/ExperimentMetricBase.h>
#include <queso/StatisticalInverseProblem.h>


namespace QUESO
{
  template<class V, class M>
  void ScenarioRunner<V,M>::runExperiment(std::vector<double> & scenario_params, std::string & filename)
  {
    m_wrapper->reinit(scenario_params);

    QUESO::GenericVectorRV<V,M> post("post_", m_prior->imageSet());

    QUESO::StatisticalInverseProblem<V,M> sip("", NULL,
                                              *(m_prior.get()),
                                              m_wrapper->get_likelihood(),
                                              post);
    
    m_experiment_metric->run(sip);

    if (m_output_rawchain)
      {
        queso_assert_not_equal_to_msg(filename,"","ERROR, you must provide a prefix to output the rawchain file");
        sip.chain().unifiedWriteContents(filename,m_output_type);
      }

    m_exp_metric_value = m_experiment_metric->evaluate(sip);
  }

}


template class QUESO::ScenarioRunner<QUESO::GslVector,QUESO::GslMatrix>;

