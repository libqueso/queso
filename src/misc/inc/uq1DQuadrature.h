//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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
// 
// $Id:$
//
//--------------------------------------------------------------------------

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

           double                     minDomainValue() const;
           double                     maxDomainValue() const;
           unsigned int               order         () const;
           const std::vector<double>& positions     () const;
           const std::vector<double>& weights       () const;
  virtual  void                       dumbRoutine   () const = 0;

protected:
  double              m_minDomainValue;
  double              m_maxDomainValue;
  unsigned int        m_order;
  std::vector<double> m_positions;
  std::vector<double> m_weights;
};

//*****************************************************
// Generic 1D quadrature class
//*****************************************************
class uqGeneric1DQuadratureClass : public uqBase1DQuadratureClass {
public:
  uqGeneric1DQuadratureClass(double minDomainValue,
                             double maxDomainValue,
                             const std::vector<double>& positions,
                             const std::vector<double>& weights);
 ~uqGeneric1DQuadratureClass();

  void dumbRoutine() const;

protected:
  using uqBase1DQuadratureClass::m_minDomainValue;
  using uqBase1DQuadratureClass::m_maxDomainValue;
  using uqBase1DQuadratureClass::m_order;
  using uqBase1DQuadratureClass::m_positions;
  using uqBase1DQuadratureClass::m_weights;
};

//*****************************************************
// Uniform/Legendre 1D quadrature class
//*****************************************************
class uqUniformLegendre1DQuadratureClass : public uqBase1DQuadratureClass {
public:
  uqUniformLegendre1DQuadratureClass(double       minDomainValue,
                                     double       maxDomainValue,
                                     unsigned int order,
                                     bool         densityIsNormalized);
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
  uqGaussianHermite1DQuadratureClass(double       mean,
                                     double       stddev,
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

//*****************************************************
// Wigner/Chebyshev1st 1D quadrature class
//*****************************************************
class uqWignerInverseChebyshev1st1DQuadratureClass : public uqBase1DQuadratureClass {
public:
  uqWignerInverseChebyshev1st1DQuadratureClass(double       minDomainValue,
                                               double       maxDomainValue,
                                               unsigned int order);
 ~uqWignerInverseChebyshev1st1DQuadratureClass();

  void dumbRoutine() const;

protected:
  using uqBase1DQuadratureClass::m_minDomainValue;
  using uqBase1DQuadratureClass::m_maxDomainValue;
  using uqBase1DQuadratureClass::m_order;
  using uqBase1DQuadratureClass::m_positions;
  using uqBase1DQuadratureClass::m_weights;
};

//*****************************************************
// Wigner/Chebyshev2nd 1D quadrature class
//*****************************************************
class uqWignerChebyshev2nd1DQuadratureClass : public uqBase1DQuadratureClass {
public:
  uqWignerChebyshev2nd1DQuadratureClass(double       minDomainValue,
                                        double       maxDomainValue,
                                        unsigned int order);
 ~uqWignerChebyshev2nd1DQuadratureClass();

  void dumbRoutine() const;

protected:
  using uqBase1DQuadratureClass::m_minDomainValue;
  using uqBase1DQuadratureClass::m_maxDomainValue;
  using uqBase1DQuadratureClass::m_order;
  using uqBase1DQuadratureClass::m_positions;
  using uqBase1DQuadratureClass::m_weights;
};
#endif // __UQ_1D_1D_QUADRATURE_H__

