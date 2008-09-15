/* uq/libs/mcmc/inc/uqMonteCarloDC.h
 *
 * Copyright (C) 2008 The PECOS Team, http://queso.ices.utexas.edu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_MC_DC_H__
#define __UQ_MC_DC_H__

/*! A templated class that implements a Monte Carlo Distribution Calculator
 */
template <class P_V,class P_M,class Q_V,class Q_M>
class uqMonteCarloDCClass
{
public:
  uqMonteCarloDCClass(const uqEnvironmentClass&                           env,                     /*! The QUESO toolkit environment.   */
                      const char*                                         prefix,                  /*! Prefix for the validation phase. */
                      const uqParamSpaceClass          <P_V,P_M>&         paramSpace,              /*! The parameter space.             */
                      const uqQoISpaceClass            <Q_V,Q_M>&         qoiSpace,                /*! The QoI space.                   */
                      const uqProbDensity_BaseClass    <P_V,P_M>*         propagParamDensityObj,   /*! */
                      const uqSampleGenerator_BaseClass<P_V,P_M>*         propagParamGeneratorObj, /*! */
                      const uqQoIFunction_BaseClass    <P_V,P_M,Q_V,Q_M>& qoiFunctionObj);         /*! */
 ~uqMonteCarloDCClass();

  void propagateUncertainty();

  void print               (std::ostream& os) const;
};

template<class P_V,class P_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const uqMonteCarloDCClass<P_V,P_M,Q_V,Q_M>& obj);

template <class P_V,class P_M,class Q_V,class Q_M>
uqMonteCarloDCClass<P_V,P_M,Q_V,Q_M>::uqMonteCarloDCClass(
  const uqEnvironmentClass&                           env,
  const char*                                         prefix,
  const uqParamSpaceClass          <P_V,P_M>&         paramSpace,
  const uqQoISpaceClass            <Q_V,Q_M>&         qoiSpace,
  const uqProbDensity_BaseClass    <P_V,P_M>*         propagParamDensityObj,
  const uqSampleGenerator_BaseClass<P_V,P_M>*         propagParamGeneratorObj,
  const uqQoIFunction_BaseClass    <P_V,P_M,Q_V,Q_M>& qoiFunctionObj)
{
}

template <class P_V,class P_M,class Q_V,class Q_M>
uqMonteCarloDCClass<P_V,P_M,Q_V,Q_M>::~uqMonteCarloDCClass()
{
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqMonteCarloDCClass<P_V,P_M,Q_V,Q_M>::propagateUncertainty()
{
  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqMonteCarloDCClass<P_V,P_M,Q_V,Q_M>::print(std::ostream& os) const
{
  return;
}

template<class P_V,class P_M,class Q_V,class Q_M> 
std::ostream& operator<<(std::ostream& os, const uqMonteCarloDCClass<P_V,P_M,Q_V,Q_M>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_MC_DC_H__
