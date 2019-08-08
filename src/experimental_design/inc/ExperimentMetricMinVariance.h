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

#ifndef QUESO_EXPERIMENT_METRIC_MIN_VARIANCE
#define QUESO_EXPERIMENT_METRIC_MIN_VARIANCE

// Base class
#include <queso/ExperimentMetricBase.h>

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/StatisticalInverseProblem.h>

namespace QUESO
{
  //! Measures the Expected Information Gain (EIG) of the given experiment
  template<class V = GslVector, class M = GslMatrix>
  class ExperimentMetricMinVariance : public ExperimentMetricBase<V,M>
  {
  public:

    //! Constructor
    /*! */
    ExperimentMetricMinVariance(GslVector & initial, GslMatrix & cov)
    : m_paramInitials(initial),
      m_propCovMatrix(cov)
    {}

    //! In order to get the EIG, we need to solve the SIP with MLSampling
    virtual void run(StatisticalInverseProblem<V,M> & sip);

    //! Return the EIG calculated by the SIP object
    virtual double evaluate(StatisticalInverseProblem<V,M> & sip);

  private:
    GslVector & m_paramInitials;
    GslMatrix & m_propCovMatrix; 

  };

}

#endif //QUESO_EXPERIMENT_METRIC_MIN_VARIANCE

