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

#ifndef UQ_SIMULATION_STORAGE_H
#define UQ_SIMULATION_STORAGE_H

#include <queso/Environment.h>
#include <queso/VectorSpace.h>

namespace QUESO {

class GslVector;
class GslMatrix;

template <class S_V = GslVector, class S_M = GslMatrix, class P_V = GslVector, class P_M = GslMatrix, class Q_V = GslVector, class Q_M = GslMatrix>
class SimulationStorage
{
public:
  SimulationStorage(const VectorSpace<S_V,S_M>& scenarioSpace,
                           const VectorSpace<P_V,P_M>& parameterSpace,
                           const VectorSpace<Q_V,Q_M>& outputSpace,
                                 unsigned int                 numSimulations);
 ~SimulationStorage();

        void                         addSimulation        (const S_V& scenarioVec, const P_V& parameterVec, const Q_V& outputVec);
        unsigned int                 numSimulations       () const;
  const std::vector<const S_V* >&    xs_asterisks_original() const;
  const std::vector<const P_V* >&    ts_asterisks_original() const;
  const VectorSpace<S_V,S_M>& scenarioSpace        () const;
  const VectorSpace<P_V,P_M>& parameterSpace       () const;
  const VectorSpace<Q_V,Q_M>& outputSpace          () const;
  const S_V&                         scenarioVec_original (unsigned int simulationId) const;
  const P_V&                         parameterVec_original(unsigned int simulationId) const;
  const Q_V&                         outputVec_original   (unsigned int simulationId) const;
  const Q_V&                         etaVec_original      () const;

  const BaseEnvironment&      env                  () const;
        void                         print                (std::ostream& os) const;
  friend std::ostream& operator<<(std::ostream& os, const SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>& obj)
  {
    obj.print(os);

    return os;
  }


private:
  // Private variables
  const BaseEnvironment&      m_env;
  const VectorSpace<S_V,S_M>& m_scenarioSpace;
  const VectorSpace<P_V,P_M>& m_parameterSpace;
  const VectorSpace<Q_V,Q_M>& m_outputSpace;

        unsigned int                 m_paper_m;
        unsigned int                 m_paper_n_eta;

        unsigned int                 m_addId;
        std::vector<const S_V* >     m_scenarioVecs_original;
        std::vector<const P_V* >     m_parameterVecs_original;
        std::vector<const Q_V* >     m_outputVecs_original;

        VectorSpace<Q_V,Q_M>* m_eta_space;
        Q_V*                         m_etaVec_original;
};

}  // End namespace QUESO

#endif // UQ_SIMULATION_STORAGE_H
