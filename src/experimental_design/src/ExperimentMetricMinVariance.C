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
#include <queso/ExperimentMetricMinVariance.h>

// QUESO
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/StatisticalInverseProblem.h>

namespace QUESO
{
  template<class V, class M>
  void ExperimentMetricMinVariance<V,M>::run(StatisticalInverseProblem<V,M> & sip)
  {
    sip.solveWithBayesMetropolisHastings(NULL,m_paramInitials,&m_propCovMatrix);
  }

  template<class V, class M>
  double ExperimentMetricMinVariance<V,M>::evaluate(StatisticalInverseProblem<V,M> & sip)
  {
    GslVector var = sip.chain().unifiedSampleVariancePlain();

    double minvar = 0.0;
    for (unsigned int i=0; i<var.sizeGlobal(); ++i)
      if ( (i == 0) || (var[i] < minvar) )
        minvar = var[i];

    return -1.0*minvar; // by returning the negative variance, the lowest variance becomes the max metric value
  }

}

template class QUESO::ExperimentMetricMinVariance<QUESO::GslVector,QUESO::GslMatrix>;
