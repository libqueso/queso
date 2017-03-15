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

#include <queso/TKGroup.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Default constructor ------------------------------
template<class V, class M>
BaseTKGroup<V,M>::BaseTKGroup()
  :
  m_emptyEnv             (new EmptyEnvironment()),
  m_env                  (*m_emptyEnv),
  m_prefix               (""),
  m_vectorSpace          (NULL),
  m_scales               (),
  m_preComputingPositions(),
  m_rvs                  (),
  m_stageId              (0)
{
}
// Constructor with values---------------------------
template<class V, class M>
BaseTKGroup<V,M>::BaseTKGroup(
  const char*                    prefix,
  const VectorSpace<V,M>& vectorSpace,
  const std::vector<double>&     scales)
  :
  m_emptyEnv             (NULL),
  m_env                  (vectorSpace.env()),
  m_prefix               (prefix),
  m_vectorSpace          (&vectorSpace),
  m_scales               (scales.size(),1.),
  m_preComputingPositions(scales.size()+1,NULL), // Yes, +1
  m_rvs                  (scales.size(),NULL), // IMPORTANT: it stays like this for scaledTK, but it will be overwritten to '+1' by hessianTK constructor
  m_stageId              (0)
{
  for (unsigned int i = 0; i < m_scales.size(); ++i) {
    m_scales[i] = scales[i];
  }
}
// Destructor ---------------------------------------
template<class V, class M>
BaseTKGroup<V,M>::~BaseTKGroup()
{
  for (unsigned int i = 0; i < m_rvs.size(); ++i) {
    if (m_rvs[i]) delete m_rvs[i];
  }
  for (unsigned int i = 0; i < m_preComputingPositions.size(); ++i) {
    if (m_preComputingPositions[i]) delete m_preComputingPositions[i];
  }
  if (m_emptyEnv) delete m_emptyEnv;
}
// Misc methods--------------------------------------
template<class V, class M>
const BaseEnvironment&
BaseTKGroup<V,M>::env() const
{
  return m_env;
}
//---------------------------------------------------
template<class V, class M>
const V&
BaseTKGroup<V,M>::preComputingPosition(unsigned int stageId) const
{
  queso_require_greater_msg(m_preComputingPositions.size(), stageId, "m_preComputingPositions.size() <= stageId");

  queso_require_msg(m_preComputingPositions[stageId], "m_preComputingPositions[stageId] == NULL");

  return *m_preComputingPositions[stageId];
}
//---------------------------------------------------
template<class V, class M>
bool
BaseTKGroup<V,M>::setPreComputingPosition(const V& position, unsigned int stageId)
{
  queso_require_greater_msg(m_preComputingPositions.size(), stageId, "m_preComputingPositions.size() <= stageId");

  queso_require_msg(!(m_preComputingPositions[stageId]), "m_preComputingPositions[stageId] != NULL");

  m_preComputingPositions[stageId] = new V(position);

  return true;
}
//---------------------------------------------------
template<class V, class M>
void
BaseTKGroup<V,M>::clearPreComputingPositions()
{
  for (unsigned int i = 0; i < m_preComputingPositions.size(); ++i) {
    if (m_preComputingPositions[i]) {
      delete m_preComputingPositions[i];
      m_preComputingPositions[i] = NULL;
    }
  }

  return;
}

template <class V, class M>
unsigned int
BaseTKGroup<V, M>::set_dr_stage(unsigned int stageId)
{
  return this->m_stageId;
}

// I/O methods---------------------------------------
template<class V, class M>
void
BaseTKGroup<V,M>::print(std::ostream& os) const
{
  os << "In BaseTKGroup<V,M>::print()"
     << ": nothing to be printed" << std::endl;
  return;
}

}  // End namespace QUESO

template class QUESO::BaseTKGroup<QUESO::GslVector, QUESO::GslMatrix>;
