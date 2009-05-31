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

#ifndef __UQ_CHAIN_POSITION_DATA_H__
#define __UQ_CHAIN_POSITION_DATA_H__

#include <uqEnvironment.h>

template <class V>
class uqMarkovChainPositionDataClass
{
public:
  uqMarkovChainPositionDataClass(const uqBaseEnvironmentClass& env);
  uqMarkovChainPositionDataClass(const uqBaseEnvironmentClass& env,
                                 const V& vecValues,
                                 bool     outOfTargetSupport,
                                 double   logTarget);
  uqMarkovChainPositionDataClass(const uqMarkovChainPositionDataClass<V>& rhs);
 ~uqMarkovChainPositionDataClass();

  uqMarkovChainPositionDataClass<V>& operator= (const uqMarkovChainPositionDataClass<V>& rhs);

  const V& vecValues         () const;
  bool     outOfTargetSupport() const;
  double   logTarget         () const;

  void     set               (const V& vecValues,
                              bool     outOfTargetSupport,
                              double   logTarget);

  void     print             (std::ostream& os) const;

private:
  const uqBaseEnvironmentClass& m_env;
  V*     m_vecValues;
  bool   m_outOfTargetSupport;
  double m_logTarget;
};

template <class V>
uqMarkovChainPositionDataClass<V>::uqMarkovChainPositionDataClass(const uqBaseEnvironmentClass& env)
  :
  m_env               (env),
  m_vecValues         (NULL),
  m_outOfTargetSupport(false),
  m_logTarget         (0.)
{
}

template <class V>
uqMarkovChainPositionDataClass<V>::uqMarkovChainPositionDataClass(
  const uqBaseEnvironmentClass& env,
  const V& vecValues,
  bool     outOfTargetSupport,
  double   logTarget)
  :
  m_env               (env),
  m_vecValues         (new V(vecValues)),
  m_outOfTargetSupport(outOfTargetSupport),
  m_logTarget         (logTarget)
{
}

template <class V>
uqMarkovChainPositionDataClass<V>::uqMarkovChainPositionDataClass(const uqMarkovChainPositionDataClass<V>& rhs)
  :
  m_env               (rhs.m_env               ),
  m_vecValues         (new V(*rhs.m_vecValues )),
  m_outOfTargetSupport(rhs.m_outOfTargetSupport),
  m_logTarget         (rhs.m_logTarget         )
{
}

template <class V>
uqMarkovChainPositionDataClass<V>::~uqMarkovChainPositionDataClass()
{
  if (m_vecValues) delete m_vecValues;
}

template <class V>
uqMarkovChainPositionDataClass<V>&
uqMarkovChainPositionDataClass<V>::operator=(const uqMarkovChainPositionDataClass<V>& rhs)
{
  if (m_vecValues == NULL) m_vecValues = new V(*rhs.m_vecValues);
  else                    *m_vecValues = *rhs.m_vecValues;
  m_outOfTargetSupport = rhs.m_outOfTargetSupport;
  m_logTarget          = rhs.m_logTarget;

  return *this;
}

template <class V>
const V&
uqMarkovChainPositionDataClass<V>::vecValues() const
{
  UQ_FATAL_TEST_MACRO((m_vecValues == NULL),
                      m_env.fullRank(),
                      "uqMarkovChainPositionDataClass<V>::vecValues()",
                      "m_vecValues is NULL");
  return *m_vecValues;
}

template <class V>
bool
uqMarkovChainPositionDataClass<V>::outOfTargetSupport() const
{
  return m_outOfTargetSupport;
}

template <class V>
double
uqMarkovChainPositionDataClass<V>::logTarget() const
{
  return m_logTarget;
}

template <class V>
void
uqMarkovChainPositionDataClass<V>::set(
  const V& vecValues,
  bool     outOfTargetSupport,
  double   logTarget)
{
  if (m_vecValues == NULL) m_vecValues = new V(vecValues);
  else                    *m_vecValues = vecValues;
  m_outOfTargetSupport = outOfTargetSupport;
  m_logTarget          = logTarget;

  return;
}

template <class V>
std::ostream& operator<<(std::ostream& os, const uqMarkovChainPositionDataClass<V>& obj)
{
  obj.print(os);

  return os;
}

#endif // __UQ_CHAIN_POSITION_DATA_H__
