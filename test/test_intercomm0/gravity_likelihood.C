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

#include <cmath>

#include <queso/GenericScalarFunction.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/UniformVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/ScalarFunction.h>
#include <queso/VectorSet.h>

#include <gravity_likelihood.h>

template<class V, class M>
Likelihood<V, M>::Likelihood(const char * prefix,
    const QUESO::VectorSet<V, M> & domain)
  : QUESO::BaseScalarFunction<V, M>(prefix, domain),
    m_heights(0),
    m_times(0),
    m_stdDevs(0)
{
  double const heights[] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110,
                            120, 130, 140};

  double const times  [] = {1.41, 2.14, 2.49, 2.87, 3.22, 3.49, 3.81, 4.07,
                            4.32, 4.47, 4.75, 4.99, 5.16, 5.26};

  double const stdDevs[] = {0.020, 0.120, 0.020, 0.010, 0.030, 0.010, 0.030,
                            0.030, 0.030, 0.050, 0.010, 0.040, 0.010, 0.09};

  std::size_t const n = sizeof(heights) / sizeof(*heights);
  m_heights.assign(heights, heights + n);
  m_times.assign(times, times + n);
  m_stdDevs.assign(stdDevs, stdDevs + n);
}

template<class V, class M>
Likelihood<V, M>::~Likelihood()
{
  // Deconstruct here
}

template<class V, class M>
double
Likelihood<V, M>::lnValue(const V & domainVector, const V * domainDirection,
    V * gradVector, M * hessianMatrix, V * hessianEffect) const
{
  double g = domainVector[0];

  double misfitValue = 0.0;
  for (unsigned int i = 0; i < m_heights.size(); ++i) {
    double modelTime = std::sqrt(2.0 * m_heights[i] / g);
    double ratio = (modelTime - m_times[i]) / m_stdDevs[i];
    misfitValue += ratio * ratio;
  }

  return -0.5 * misfitValue;
}

template<class V, class M>
double
Likelihood<V, M>::actualValue(const V & domainVector,
    const V * domainDirection, V * gradVector, M * hessianMatrix,
    V * hessianEffect) const
{
  return std::exp(this->lnValue(domainVector, domainDirection, gradVector,
        hessianMatrix, hessianEffect));
}

template class Likelihood<QUESO::GslVector, QUESO::GslMatrix>;
