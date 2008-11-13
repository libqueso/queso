/* uq/libs/queso/inc/uqMarkovChainPosition.h
 *
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_CHAIN_POSITION_H__
#define __UQ_CHAIN_POSITION_H__

#include <uqEnvironment.h>

template <class V>
class uqMarkovChainPositionClass
{
public:
  uqMarkovChainPositionClass(const uqBaseEnvironmentClass& env);
  uqMarkovChainPositionClass(const uqBaseEnvironmentClass& env,
                             const V& vecValues,
                             bool     outOfBounds,
                             double   logTarget);
  uqMarkovChainPositionClass(const uqMarkovChainPositionClass<V>& rhs);
 ~uqMarkovChainPositionClass();

  uqMarkovChainPositionClass<V>& operator= (const uqMarkovChainPositionClass<V>& rhs);

  const V& vecValues  () const;
  bool     outOfBounds() const;
  double   logTarget  () const;

  void     set        (const V& vecValues,
                       bool     outOfBounds,
                       double   logTarget);

  void     print      (std::ostream& os) const;

private:
  const uqBaseEnvironmentClass& m_env;
  V*     m_vecValues;
  bool   m_outOfBounds;
  double m_logTarget;
};

template <class V>
uqMarkovChainPositionClass<V>::uqMarkovChainPositionClass(const uqBaseEnvironmentClass& env)
  :
  m_env        (env),
  m_vecValues(NULL),
  m_outOfBounds(false),
  m_logTarget  (0.)
{
}

template <class V>
uqMarkovChainPositionClass<V>::uqMarkovChainPositionClass(
  const uqBaseEnvironmentClass& env,
  const V& vecValues,
  bool     outOfBounds,
  double   logTarget)
  :
  m_env        (env),
  m_vecValues(new V(vecValues)),
  m_outOfBounds(outOfBounds),
  m_logTarget  (logTarget)
{
}

template <class V>
uqMarkovChainPositionClass<V>::uqMarkovChainPositionClass(const uqMarkovChainPositionClass<V>& rhs)
  :
  m_env        (rhs.m_env                ),
  m_vecValues(new V(*rhs.m_vecValues)),
  m_outOfBounds(rhs.m_outOfBounds        ),
  m_logTarget  (rhs.m_logTarget          )
{
}

template <class V>
uqMarkovChainPositionClass<V>::~uqMarkovChainPositionClass()
{
  if (m_vecValues) delete m_vecValues;
}

template <class V>
uqMarkovChainPositionClass<V>&
uqMarkovChainPositionClass<V>::operator=(const uqMarkovChainPositionClass<V>& rhs)
{
  if (m_vecValues == NULL) m_vecValues = new V(*rhs.m_vecValues);
  else                    *m_vecValues = *rhs.m_vecValues;
  m_outOfBounds = rhs.m_outOfBounds;
  m_logTarget   = rhs.m_logTarget;

  return *this;
}

template <class V>
const V&
uqMarkovChainPositionClass<V>::vecValues() const
{
  UQ_FATAL_TEST_MACRO((m_vecValues == NULL),
                      m_env.rank(),
                      "uqMarkovChainPositionClass<V>::vecValues()",
                      "m_vecValues is NULL");
  return *m_vecValues;
}

template <class V>
bool
uqMarkovChainPositionClass<V>::outOfBounds() const
{
  return m_outOfBounds;
}

template <class V>
double
uqMarkovChainPositionClass<V>::logTarget() const
{
  return m_logTarget;
}

template <class V>
void
uqMarkovChainPositionClass<V>::set(
  const V& vecValues,
  bool     outOfBounds,
  double   logTarget)
{
  if (m_vecValues == NULL) m_vecValues = new V(vecValues);
  else                    *m_vecValues = vecValues;
  m_outOfBounds = outOfBounds;
  m_logTarget   = logTarget;

  return;
}

template <class V>
std::ostream& operator<<(std::ostream& os, const uqMarkovChainPositionClass<V>& obj)
{
  obj.print(os);

  return os;
}

#endif // __UQ_CHAIN_POSITION_H__
