//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
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

#include <queso/SimulationStorage.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::SimulationStorage(
  const VectorSpace<S_V,S_M>& scenarioSpace,
  const VectorSpace<P_V,P_M>& parameterSpace,
  const VectorSpace<Q_V,Q_M>& outputSpace,
  unsigned int                       numSimulations)
  :
  m_env                   (scenarioSpace.env()),
  m_scenarioSpace         (scenarioSpace),
  m_parameterSpace        (parameterSpace),
  m_outputSpace           (outputSpace),
  m_paper_m               (numSimulations),
  m_paper_n_eta           (outputSpace.dimLocal()),
  m_addId                 (0),
  m_scenarioVecs_original (m_paper_m, (S_V*) NULL),
  m_parameterVecs_original(m_paper_m, (P_V*) NULL),
  m_outputVecs_original   (m_paper_m, (Q_V*) NULL),
  m_eta_space             (NULL),
  m_etaVec_original       (NULL)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << "\n  m_paper_m = "     << m_paper_m
                            << "\n  m_paper_n_eta = " << m_paper_n_eta
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::constructor()"
                            << std::endl;
  }
}

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::~SimulationStorage()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::destructor()..."
                            << std::endl;
  }

  delete m_etaVec_original;
  delete m_eta_space;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::destructor()"
                            << std::endl;
  }
}

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
void
SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::addSimulation(const S_V& scenarioVec, const P_V& parameterVec, const Q_V& outputVec)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::addSimulation()"
                            << ": m_addId = " << m_addId
                            << std::endl;
  }

  queso_require_less_msg(m_addId, m_paper_m, "too many adds...");

  m_scenarioVecs_original [m_addId] = &scenarioVec;
  m_parameterVecs_original[m_addId] = &parameterVec;
  m_outputVecs_original   [m_addId] = &outputVec;
  m_addId++;

  if (m_addId == m_paper_m) {
    //***********************************************************************
    // Form 'etaVec_original'
    //***********************************************************************
    m_eta_space = new VectorSpace<Q_V,Q_M>(m_env, "m_eta_simul_storage", m_paper_m * m_paper_n_eta, NULL),
    m_etaVec_original = new Q_V(m_eta_space->zeroVector());
    m_etaVec_original->cwSetConcatenated(m_outputVecs_original);

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
      *m_env.subDisplayFile() << "KEY In SimulationStorage<S_V,S_M,D_V,D_M>::addSimulation()"
                              << ": m_addId = " << m_addId
                              << ", populated etaVec_original of size " << m_etaVec_original->sizeLocal()
  //<< ", m_etaVec_original = "  << m_etaVec_original
  //<< ", *m_etaVec_original = " << *m_etaVec_original
                              << std::endl;
    }
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::addSimulation()"
                            << ": m_addId = " << m_addId
                            << std::endl;
  }

  return;
}

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
unsigned int
SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::numSimulations() const
{
  return m_paper_m;
}

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
const std::vector<const S_V* >&
SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::xs_asterisks_original() const
{
  return m_scenarioVecs_original;
}

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
const std::vector<const P_V* >&
SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::ts_asterisks_original() const
{
  return m_parameterVecs_original;
}

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
const VectorSpace<S_V,S_M>&
SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::scenarioSpace() const
{
  return m_scenarioSpace;
}

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
const VectorSpace<P_V,P_M>&
SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::parameterSpace() const
{
  return m_parameterSpace;
}

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
const VectorSpace<Q_V,Q_M>&
SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::outputSpace() const
{
  return m_outputSpace;
}

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
const S_V&
SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::scenarioVec_original(unsigned int simulationId) const
{
  queso_require_less_msg(simulationId, m_scenarioVecs_original.size(), "simulationId is too large");

  queso_require_msg(m_scenarioVecs_original[simulationId], "vector is NULL");

  return *(m_scenarioVecs_original[simulationId]);
}

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
const P_V&
SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::parameterVec_original(unsigned int simulationId) const
{
  queso_require_less_msg(simulationId, m_parameterVecs_original.size(), "simulationId is too large");

  queso_require_msg(m_parameterVecs_original[simulationId], "vector is NULL");

  return *(m_parameterVecs_original[simulationId]);
}

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
const Q_V&
SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::outputVec_original(unsigned int simulationId) const
{
  queso_require_less_msg(simulationId, m_outputVecs_original.size(), "simulationId is too large");

  queso_require_msg(m_outputVecs_original[simulationId], "vector is NULL");

  return *(m_outputVecs_original[simulationId]);
}

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
const Q_V&
SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::etaVec_original() const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
    *m_env.subDisplayFile() << "Entering SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::etaVec_original()"
                            << ": m_etaVec_original = " << m_etaVec_original
                            << std::endl;
  }

  queso_require_msg(m_etaVec_original, "'m_etaVec_original' is NULL");

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
    *m_env.subDisplayFile() << "In SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::etaVec_original()"
                            << ": *m_etaVec_original = " << *m_etaVec_original
                            << std::endl;
  }

  return *m_etaVec_original;
}

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
const BaseEnvironment&
  SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::env() const
{
  return m_env;
}

template<class S_V,class S_M,class P_V,class P_M,class Q_V,class Q_M>
void
SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::print(std::ostream& os) const
{
  return;
}

}  // End namespace QUESO

template class QUESO::SimulationStorage<QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix>;
