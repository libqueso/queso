/* uq/libs/mcmc/inc/uqStochasticProblem.h
 *
 * Copyright (C) 2008 The PECOS Team, http://www.ices.utexas.edu/centers/pecos
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

#ifndef __UQ_STOCHASTIC_PROBLEM_H__
#define __UQ_STOCHASTIC_PROBLEM_H__

#include <uqRandomVariable.h>

//*****************************************************
// Base class
//*****************************************************
template<class V, class M>
class uqStochasticProblem_BaseClass {
public:
  uqStochasticProblem_BaseClass(const uqEnvironmentCLass& env,
                                  char*                     prefix);
  virtual ~uqStochasticProblem_BaseClass();

  virtual       void                   prepareSolver() = 0;
  virtual       void                   solve        () = 0;
  virtual const uqRandomVariableClass& solution     () = 0;

protected:
  const uqEnvironmentClass& m_env;
        std::string         m_prefix;

        uqRandomVariableClass* m_solution;
};

template<class V, class M>
uqStochasticProblem_BaseClass<V,M>::uqStochasticProblem_BaseClass(
  double (*routinePtr)(const V& paramValues, const void* routineDataPtr),
  const void* routineDataPtr)
  :
  m_routinePtr        (routinePtr),
  m_routineDataPtr    (routineDataPtr),
  m_priorDensity      (NULL),
  m_likelihoodFunction(NULL)
{
}

template<class V, class M>
uqStochasticProblem_BaseClass<V,M>::~uqStochasticProblem_BaseClass()
{
}

#endif // __UQ_STOCHASTIC_PROBLEM_H__
