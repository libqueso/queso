/* uq/libs/mcmc/inc/uqProblem.h
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

#ifndef __UQ_PROBLEM_H__
#define __UQ_PROBLEM_H__

#include <uqProblemSlice.h>

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
class uqProblemClass
{
public:
  uqProblemClass(const uqEnvironmentClass& env);
 ~uqProblemClass();

 void setSlice(unsigned int sliceId);

 void quantify();

 void print   (std::ostream& os) const;

private:
  const uqEnvironmentClass&                                  m_env;

  po::options_description*                                   m_optionsDesc;
  std::vector<uqProblemSliceClass<P_V,P_M,L_V,L_M,Q_V,Q_M>*> m_slices;
};

template<class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const uqProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>& obj);

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
uqProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::uqProblemClass(const uqEnvironmentClass& env)
  :
  m_env        (env),
  m_optionsDesc(new po::options_description("Uncertainty Quantification Problem")),
  m_slices     (0)
{
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
uqProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::~uqProblemClass()
{
  for (unsigned int i = 0; i < m_slices.size(); ++i) {
    if (m_slices[i] != NULL) delete m_slices[i];
  }
  if (m_optionsDesc) delete m_optionsDesc;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
void
uqProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::setSlice(unsigned int sliceId)
{
  return;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
void
uqProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::quantify()
{
  return;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
void
uqProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::print(std::ostream& os) const
{
}

template<class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const uqProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_PROBLEM_H__
