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

#ifndef QUESO_EXPERIMENTAL_LIKELIHOOD_WRAPPER
#define QUESO_EXPERIMENTAL_LIKELIHOOD_WRAPPER

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/ScalarFunction.h>
#include <queso/ExperimentalLikelihoodInterface.h>

namespace QUESO
{

  //! Base class for likelihood functions that are used in experimental design
  /*!
   *  This class allows for a modifiable evaluateModel() function based
   *  for a particular set of experimental scenario parameters as part
   *  of an experimental design loop.
   *
   *  The intended use is for a user create their own likelihood class
   *  that is derived from both this class and an appropriate
   *  Likelihood class. By doing so, a user need only implement
   *  the reinit() and evaluateModel() functions themselves.
   */
  template<class V = GslVector, class M = GslMatrix>
  class ExperimentalLikelihoodWrapper
  {
  public:

    ExperimentalLikelihoodWrapper(std::shared_ptr<BaseScalarFunction<V,M>> & likelihood, std::shared_ptr<ExperimentalLikelihoodInterface<V,M>> & interface)
    : m_likelihood(likelihood),
      m_interface(interface)
    {}

    void reinit(std::vector<double> & scenario_params)
    {
      m_interface->reinit(scenario_params,this->get_likelihood());
    }
    
    BaseScalarFunction<V,M> & get_likelihood()
    {
      return *(m_likelihood.get());
    }

  private:
    
    std::shared_ptr<BaseScalarFunction<V,M>> m_likelihood;
    std::shared_ptr<ExperimentalLikelihoodInterface<V,M>> m_interface;

  };

}

#endif //QUESO_EXPERIMENTAL_LIKELIHOOD_WRAPPER

