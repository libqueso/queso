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

#include <uqParameter.h>
#include <uqEnvironment.h>

template <class V>
class uqChainPositionClass
{
public:
  uqChainPositionClass(const uqEnvironmentClass& env);
  uqChainPositionClass(const uqEnvironmentClass& env,
                       const V& paramValues,
                       bool     outOfBounds,
                       double   m2lPrior,
                       const V& misfitVector,
                       const V& lrVarianceVector,
                       const V& m2lLikelihoodVector,
                       double   logPosterior);
  uqChainPositionClass(const uqChainPositionClass<V>& rhs);
 ~uqChainPositionClass();

  uqChainPositionClass<V>& operator= (const uqChainPositionClass<V>& rhs);

  const V& paramValues        () const;
  bool     outOfBounds        () const;
  double   m2lPrior           () const;
  const V& misfitVector       () const;
  const V& lrVarianceVector   () const;
  const V& m2lLikelihoodVector() const;
  double   logPosterior       () const;

  void     set                (const V& paramValues,
                               bool     outOfBounds,
                               double   m2lPrior,
                               const V& misfitVector,
                               const V& lrVarianceVector,
                               const V& m2lLikelihoodVector,
                               double   logPosterior);

  void     print              (std::ostream& os) const;

private:
  const uqEnvironmentClass& m_env;
  V*     m_paramValues;
  bool   m_outOfBounds;
  double m_m2lPrior;
  V*     m_misfitVector;
  V*     m_lrVarianceVector;
  V*     m_m2lLikelihoodVector;
  double m_logPosterior;
};

template <class V>
uqChainPositionClass<V>::uqChainPositionClass(const uqEnvironmentClass& env)
  :
  m_env                (env),
  m_paramValues        (NULL),
  m_outOfBounds        (false),
  m_m2lPrior           (0.),
  m_misfitVector       (NULL),
  m_lrVarianceVector   (NULL),
  m_m2lLikelihoodVector(NULL),
  m_logPosterior       (0.)
{
}

template <class V>
uqChainPositionClass<V>::uqChainPositionClass(
  const uqEnvironmentClass& env,
  const V& paramValues,
  bool     outOfBounds,
  double   m2lPrior,
  const V& misfitVector,
  const V& lrVarianceVector,
  const V& m2lLikelihoodVector,
  double   logPosterior)
  :
  m_env                (env),
  m_paramValues        (new V(paramValues)),
  m_outOfBounds        (outOfBounds),
  m_m2lPrior           (m2lPrior),
  m_misfitVector       (new V(misfitVector)),
  m_lrVarianceVector   (new V(lrVarianceVector)),
  m_m2lLikelihoodVector(new V(m2lLikelihoodVector)),
  m_logPosterior       (logPosterior)
{
}

template <class V>
uqChainPositionClass<V>::uqChainPositionClass(const uqChainPositionClass<V>& rhs)
  :
  m_env                (rhs.m_env                        ),
  m_paramValues        (new V(*rhs.m_paramValues)        ),
  m_outOfBounds        (rhs.m_outOfBounds                ),
  m_m2lPrior           (rhs.m_m2lPrior                   ),
  m_misfitVector       (new V(*rhs.m_misfitVector)       ),
  m_lrVarianceVector   (new V(*rhs.m_lrVarianceVector)   ),
  m_m2lLikelihoodVector(new V(*rhs.m_m2lLikelihoodVector)),
  m_logPosterior       (rhs.m_logPosterior               )
{
}

template <class V>
uqChainPositionClass<V>::~uqChainPositionClass()
{
  if (m_m2lLikelihoodVector) delete m_m2lLikelihoodVector;
  if (m_lrVarianceVector)    delete m_lrVarianceVector;
  if (m_misfitVector)        delete m_misfitVector;
  if (m_paramValues)         delete m_paramValues;
}

template <class V>
uqChainPositionClass<V>&
uqChainPositionClass<V>::operator=(const uqChainPositionClass<V>& rhs)
{
  if (m_paramValues == NULL) m_paramValues = new V(*rhs.m_paramValues);
  else                      *m_paramValues = *rhs.m_paramValues;
  m_outOfBounds   = rhs.m_outOfBounds;
  m_m2lPrior      = rhs.m_m2lPrior;
  if (m_misfitVector          == NULL) m_misfitVector        = new V(*rhs.m_misfitVector);
  else                                *m_misfitVector        = *rhs.m_misfitVector;
  if (m_lrVarianceVector      == NULL) m_lrVarianceVector    = new V(*rhs.m_lrVarianceVector);
  else                                *m_lrVarianceVector    = *rhs.m_lrVarianceVector;
  if (m_m2lLikelihoodVector   == NULL) m_m2lLikelihoodVector = new V(*rhs.m_m2lLikelihoodVector);
  else                                *m_m2lLikelihoodVector = *rhs.m_m2lLikelihoodVector;
  m_logPosterior  = rhs.m_logPosterior;

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
uqChainPositionClass<V>::m2lPrior() const
{
  return m_m2lPrior;
}

template <class V>
const V&    
uqChainPositionClass<V>::misfitVector() const
{
  return *m_misfitVector;
}

template <class V>
const V&
uqChainPositionClass<V>::lrVarianceVector() const
{
  return *m_lrVarianceVector;
}

template <class V>
const V&    
uqChainPositionClass<V>::m2lLikelihoodVector() const
{
  return *m_m2lLikelihoodVector;
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
  double   m2lPrior,
  const V& misfitVector,
  const V& lrVarianceVector,
  const V& m2lLikelihoodVector,
  double   logPosterior)
{
  if (m_paramValues == NULL) m_paramValues = new V(paramValues);
  else                      *m_paramValues = paramValues;
  m_outOfBounds   = outOfBounds;
  m_m2lPrior      = m2lPrior;
  if (m_m2lLikelihoodVector == NULL) m_m2lLikelihoodVector  = new V(m2lLikelihoodVector);
  else                               *m_m2lLikelihoodVector = m2lLikelihoodVector;
  if (m_lrVarianceVector    == NULL) m_lrVarianceVector     = new V(lrVarianceVector);
  else                               *m_lrVarianceVector    = lrVarianceVector;
  if (m_misfitVector        == NULL) m_misfitVector         = new V(misfitVector);
  else                               *m_misfitVector        = misfitVector;
  m_logPosterior  = logPosterior;

  return;
}

template <class V>
std::ostream& operator<<(std::ostream& os, const uqChainPositionClass<V>& obj)
{
  obj.print(os);

  return os;
}

#endif // __UQ_CHAIN_POSITION_H__
