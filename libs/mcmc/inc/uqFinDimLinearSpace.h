/* uq/libs/mcmc/inc/uqFinDimLinearSpace.h
 *
 * Copyright (C) 2008 The PECOS Team, http://www.ices.utexas.edu/centers/pecos

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

#ifndef __UQ_FIN_DIM_LINEAR_SPACE_H__
#define __UQ_FIN_DIM_LINEAR_SPACE_H__

#include <uqEnvironment.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <uqDefines.h>

template <class V, class M>
class uqFinDimLinearSpaceClass
{
public:
           uqFinDimLinearSpaceClass();
           uqFinDimLinearSpaceClass(const uqEnvironmentClass& env,
                                    const char*               prefix);
  virtual ~uqFinDimLinearSpaceClass();

  virtual unsigned int      dim                       ()                 const = 0;
          V*                newVector                 ()                 const; // See template specialization
          V*                newVector                 (const V& v)       const;
          M*                newMatrix                 ()                 const; // See template specialization
          M*                newDiagMatrix             (const V& v)       const;
          M*                newDiagMatrix             (double diagValue) const; // See template specialization

  virtual void              print                     (std::ostream& os) const;

protected:
  const uqEnvironmentClass& m_env;
  std::string               m_prefix;
  unsigned int              m_dim;
};

template <class V, class M>
uqFinDimLinearSpaceClass<V,M>::uqFinDimLinearSpaceClass()
  :
  m_env(*(new uqEnvironmentClass()))
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqFinDimLinearSpaceClass<V,M>::constructor(), default",
                      "should not be used by user");
}

template <class V, class M>
uqFinDimLinearSpaceClass<V,M>::uqFinDimLinearSpaceClass(
  const uqEnvironmentClass& env,
  const char*               prefix)
  :
  m_env   (env),
  m_prefix(""),
  m_dim   (0)
{
  //std::cout << "Entering uqFinDimLinearSpaceClass<V,M>::constructor()"
  //          << std::endl;

  if ((prefix         != NULL) && 
      (strlen(prefix) != 0   )) {
    std::string tmpString(prefix);
    m_prefix = tmpString + "_";
  }

  //std::cout << "Leaving uqFinDimLinearSpaceClass<V,M>::constructor()"
  //          << std::endl;
}

template <class V, class M>
uqFinDimLinearSpaceClass<V,M>::~uqFinDimLinearSpaceClass()
{
  //std::cout << "Entering uqFinDimLinearSpaceClass<V,M>::destructor()"
  //          << std::endl;

  //std::cout << "Leaving uqFinDimLinearSpaceClass<V,M>::destructor()"
  //          << std::endl;
}

template <class V, class M>
unsigned int
uqFinDimLinearSpaceClass<V,M>::dim() const
{
  return m_dim;
}

template <class V, class M>
V*
uqFinDimLinearSpaceClass<V,M>::newVector(const V& v) const
{
  if (v.size() != m_dim) return NULL;

  return new V(v);
}

template <class V, class M>
M*
  uqFinDimLinearSpaceClass<V,M>::newDiagMatrix(const V& v) const
{
  if (v.size() != m_dim) return NULL;

  return new M(v);
}

template <class V, class M>
void
uqFinDimLinearSpaceClass<V,M>::print(std::ostream& os) const
{
  return;
}
#endif // __UQ_FIN_DIM_LINEAR_SPACE_H__

