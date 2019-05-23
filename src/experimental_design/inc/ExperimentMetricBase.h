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

#ifndef QUESO_EXPERIMENT_METRIC_BASE
#define QUESO_EXPERIMENT_METRIC_BASE

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/StatisticalInverseProblem.h>

namespace QUESO
{
  //! Base class for providing a numerical metric to an experimental scenario
  /*!
   *  This class provides a way to measure the "usefulness" of a given
   *  experiment compared to other proposed experiments (see GridSearchExperimentalDesign).
   *  Subclasses will specify the solution type for the statistical inverse problem
   *  as well as the metric by which the given experiment can be measured and compared.
   */
  template<class V = GslVector, class M = GslMatrix>
  class ExperimentMetricBase
  {
  public:

    ExperimentMetricBase() {}

    virtual ~ExperimentMetricBase(){}

    //! Runs the statistical inverse problem. Pure virtual
    //! so the user can select what type of solver they want to use
    virtual void run(StatisticalInverseProblem<V,M> & sip) =0;

    //! Provides a numerical metric for evaluating the experiment
    virtual double evaluate(StatisticalInverseProblem<V,M> & sip) =0;

  };

}

#endif //QUESO_EXPERIMENT_METRIC_BASE

