/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * $Id$
 *
 * Brief description of this file: 
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __UQ_1D_1D_QUADRATURE_H__
#define __UQ_1D_1D_QUADRATURE_H__

#include <uqEnvironment.h>
#include <uqDefines.h>
#include <vector>
#include <math.h>
#include <fstream>

//*****************************************************
// Base 1D quadrature class
//*****************************************************
class
uqBase1DQuadratureClass {
public:
           uqBase1DQuadratureClass(double minDomainValue,
                                   double maxDomainValue,
                                   unsigned int order);
  virtual ~uqBase1DQuadratureClass();

           double minDomainValue() const;
           double maxDomainValue() const;
           const std::vector<double>& positions() const;
           const std::vector<double>& weights  () const;

  virtual  void dumbRoutine() const = 0;
protected:
  double              m_minDomainValue;
  double              m_maxDomainValue;
  unsigned int        m_order;
  std::vector<double> m_positions;
  std::vector<double> m_weights;
};

//*****************************************************
// Uniform/Legendre 1D quadrature class
//*****************************************************
class uqUniformLegendre1DQuadratureClass : public uqBase1DQuadratureClass {
public:
  uqUniformLegendre1DQuadratureClass(double minDomainValue,
                                     double maxDomainValue,
                                     unsigned int order);
 ~uqUniformLegendre1DQuadratureClass();

  void dumbRoutine() const;

protected:
  using uqBase1DQuadratureClass::m_minDomainValue;
  using uqBase1DQuadratureClass::m_maxDomainValue;
  using uqBase1DQuadratureClass::m_order;
  using uqBase1DQuadratureClass::m_positions;
  using uqBase1DQuadratureClass::m_weights;
};

//*****************************************************
// Gaussian/Hermite 1D quadrature class
//*****************************************************
class uqGaussianHermite1DQuadratureClass : public uqBase1DQuadratureClass {
public:
  uqGaussianHermite1DQuadratureClass(double mean,
                                     double stddev,
                                     unsigned int order);
 ~uqGaussianHermite1DQuadratureClass();

  void dumbRoutine() const;

protected:
  using uqBase1DQuadratureClass::m_minDomainValue;
  using uqBase1DQuadratureClass::m_maxDomainValue;
  using uqBase1DQuadratureClass::m_order;
  using uqBase1DQuadratureClass::m_positions;
  using uqBase1DQuadratureClass::m_weights;

  double m_mean;
  double m_stddev;
};
#endif // __UQ_1D_1D_QUADRATURE_H__

