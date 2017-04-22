//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2017 The PECOS Development Team
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

#include <queso/MarkovChainPositionData.h>
#include <queso/GslVector.h>

namespace QUESO {

// Constructor 1 -----------------------------------
template <class V>
MarkovChainPositionData<V>::MarkovChainPositionData(const BaseEnvironment& env)
  :
  m_env               (env),
  m_vecValues         (NULL),
  m_outOfTargetSupport(false),
  m_logLikelihood     (0.),
  m_logTarget         (0.)
{
}
// Constructor 2 -----------------------------------
template <class V>
MarkovChainPositionData<V>::MarkovChainPositionData(
  const BaseEnvironment& env,
  const V& vecValues,
  bool     outOfTargetSupport,
  double   logLikelihood,
  double   logTarget)
  :
  m_env               (env),
  m_vecValues         (new V(vecValues)),
  m_outOfTargetSupport(outOfTargetSupport),
  m_logLikelihood     (logLikelihood),
  m_logTarget         (logTarget)
{
}
// Copy constructor---------------------------------
template <class V>
MarkovChainPositionData<V>::MarkovChainPositionData(const MarkovChainPositionData<V>& rhs)
  :
  m_env               (rhs.m_env               ),
  m_vecValues         (new V(*rhs.m_vecValues )),
  m_outOfTargetSupport(rhs.m_outOfTargetSupport),
  m_logLikelihood     (rhs.m_logLikelihood     ),
  m_logTarget         (rhs.m_logTarget         )
{
}
// Destructor---------------------------------------
template <class V>
MarkovChainPositionData<V>::~MarkovChainPositionData()
{
  if (m_vecValues) delete m_vecValues;
}
// Set methods--------------------------------------
template <class V>
MarkovChainPositionData<V>&
MarkovChainPositionData<V>::operator=(const MarkovChainPositionData<V>& rhs)
{
  if (m_vecValues == NULL) m_vecValues = new V(*rhs.m_vecValues);
  else                    *m_vecValues = *rhs.m_vecValues;
  m_outOfTargetSupport = rhs.m_outOfTargetSupport;
  m_logLikelihood      = rhs.m_logLikelihood;
  m_logTarget          = rhs.m_logTarget;

  return *this;
}

// Statistical methods-------------------------------
template <class V>
const V&
MarkovChainPositionData<V>::vecValues() const
{
  queso_require_msg(m_vecValues, "m_vecValues is NULL");
  return *m_vecValues;
}

//--------------------------------------------------
template <class V>
bool
MarkovChainPositionData<V>::outOfTargetSupport() const
{
  return m_outOfTargetSupport;
}

//--------------------------------------------------
template <class V>
double
MarkovChainPositionData<V>::logLikelihood() const
{
  return m_logLikelihood;
}

//--------------------------------------------------
template <class V>
double
MarkovChainPositionData<V>::logTarget() const
{
  return m_logTarget;
}
//--------------------------------------------------
template <class V>
void
MarkovChainPositionData<V>::set(
  const V& vecValues,
  bool     outOfTargetSupport,
  double   logLikelihood,
  double   logTarget)
{
  if (m_vecValues == NULL) m_vecValues = new V(vecValues);
  else                    *m_vecValues = vecValues;
  m_outOfTargetSupport = outOfTargetSupport;
  m_logLikelihood      = logLikelihood;
  m_logTarget          = logTarget;

  return;
}

}  // End namespace QUESO

template class QUESO::MarkovChainPositionData<QUESO::GslVector>;
