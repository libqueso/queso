/* uq/libs/mcmc/inc/uqChainPosition.h
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

#ifndef __UQ_CHAIN_POSITION_H__
#define __UQ_CHAIN_POSITION_H__

#include <uqEnvironment.h>

template <class V>
class uqChainPositionClass
{
public:
  uqChainPositionClass(const uqEnvironmentClass& env);
  uqChainPositionClass(const uqEnvironmentClass& env,
                       const V& paramValues,
                       bool     outOfBounds,
                       double   logPosterior);
  uqChainPositionClass(const uqChainPositionClass<V>& rhs);
 ~uqChainPositionClass();

  uqChainPositionClass<V>& operator= (const uqChainPositionClass<V>& rhs);

  const V& paramValues () const;
  bool     outOfBounds () const;
  double   logPosterior() const;

  void     set         (const V& paramValues,
                        bool     outOfBounds,
                        double   logPosterior);

  void     print       (std::ostream& os) const;

private:
  const uqEnvironmentClass& m_env;
  V*     m_paramValues;
  bool   m_outOfBounds;
  double m_logPosterior;
};

template <class V>
uqChainPositionClass<V>::uqChainPositionClass(const uqEnvironmentClass& env)
  :
  m_env         (env),
  m_paramValues (NULL),
  m_outOfBounds (false),
  m_logPosterior(0.)
{
}

template <class V>
uqChainPositionClass<V>::uqChainPositionClass(
  const uqEnvironmentClass& env,
  const V& paramValues,
  bool     outOfBounds,
  double   logPosterior)
  :
  m_env         (env),
  m_paramValues (new V(paramValues)),
  m_outOfBounds (outOfBounds),
  m_logPosterior(logPosterior)
{
}

template <class V>
uqChainPositionClass<V>::uqChainPositionClass(const uqChainPositionClass<V>& rhs)
  :
  m_env         (rhs.m_env                ),
  m_paramValues (new V(*rhs.m_paramValues)),
  m_outOfBounds (rhs.m_outOfBounds        ),
  m_logPosterior(rhs.m_logPosterior       )
{
}

template <class V>
uqChainPositionClass<V>::~uqChainPositionClass()
{
  if (m_paramValues) delete m_paramValues;
}

template <class V>
uqChainPositionClass<V>&
uqChainPositionClass<V>::operator=(const uqChainPositionClass<V>& rhs)
{
  if (m_paramValues == NULL) m_paramValues = new V(*rhs.m_paramValues);
  else                      *m_paramValues = *rhs.m_paramValues;
  m_outOfBounds  = rhs.m_outOfBounds;
  m_logPosterior = rhs.m_logPosterior;

  return *this;
}

template <class V>
const V&
uqChainPositionClass<V>::paramValues() const
{
  UQ_FATAL_TEST_MACRO((m_paramValues == NULL),
                      m_env.rank(),
                      "uqChainPositionClass<V>::paramValues()",
                      "m_paramValues is NULL");
  return *m_paramValues;
}

template <class V>
bool
uqChainPositionClass<V>::outOfBounds() const
{
  return m_outOfBounds;
}

template <class V>
double
uqChainPositionClass<V>::logPosterior() const
{
  return m_logPosterior;
}

template <class V>
void
uqChainPositionClass<V>::set(
  const V& paramValues,
  bool     outOfBounds,
  double   logPosterior)
{
  if (m_paramValues == NULL) m_paramValues = new V(paramValues);
  else                      *m_paramValues = paramValues;
  m_outOfBounds  = outOfBounds;
  m_logPosterior = logPosterior;

  return;
}

template <class V>
std::ostream& operator<<(std::ostream& os, const uqChainPositionClass<V>& obj)
{
  obj.print(os);

  return os;
}

#endif // __UQ_CHAIN_POSITION_H__
