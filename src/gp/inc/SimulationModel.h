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

#ifndef UQ_SIMULATION_MODEL_H
#define UQ_SIMULATION_MODEL_H

#include <queso/SimulationModelOptions.h>
#include <queso/SimulationStorage.h>
#include <queso/SequenceOfVectors.h>
#include <queso/Environment.h>

namespace QUESO {

class GslVector;
class GslMatrix;

template <class S_V = GslVector, class S_M = GslMatrix, class P_V = GslVector, class P_M = GslMatrix, class Q_V = GslVector, class Q_M = GslMatrix>
class SimulationModel
{
public:
  SimulationModel(const char*                                              prefix,
                         const SmOptionsValues*                            alternativeOptionsValues, // dakota
                         const SimulationStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>& simulStorage);
 ~SimulationModel();

        unsigned int                   numBasis             () const;
  const std::vector<const S_V* >&      xs_asterisks_standard() const;
  const S_V&                           xSeq_original_mins   () const;
  const S_V&                           xSeq_original_ranges () const;
  const std::vector<const P_V* >&      ts_asterisks_standard() const;
  const Q_V&                           etaSeq_original_mean () const;
#ifdef _GPMSA_CODE_TREATS_SIMULATION_VECTORS_IN_CHUNKS
        double                         etaSeq_chunkStd      (unsigned int chunkId) const;
#else
        double                         etaSeq_allStd        () const;
#endif
  const Q_V&                           etaVec_transformed   (const std::string& debugString) const;
  const Q_V&                           basisVec             (unsigned int basisId) const;
  const Q_M&                           Kmat_eta             () const;
  const Q_M&                           Kmat                 () const;

  const SimulationModelOptions& optionsObj           () const;
        void                           print                (std::ostream& os) const;
  friend std::ostream& operator<<(std::ostream& os,
      const SimulationModel<S_V,S_M,P_V,P_M,Q_V,Q_M>& obj) {
    obj.print(os);
    return os;
  }

private:
        unsigned int                   computePEta          (const Q_V& svdS_vec);
  // Private variables
  const BaseEnvironment&           m_env;
  const SmOptionsValues *          m_optionsObj;
        SimulationModelOptions *   m_simulationModelOptions;

        unsigned int                      m_paper_p_x;
        unsigned int                      m_paper_p_t;
        unsigned int                      m_paper_m;
        unsigned int                      m_paper_n_eta;

        VectorSpace<S_V,S_M>       m_p_x_space;
        SequenceOfVectors<S_V,S_M> m_xSeq_original;
        S_V                               m_xSeq_original_mins;
        S_V                               m_xSeq_original_maxs;
        S_V                               m_xSeq_original_ranges;
        SequenceOfVectors<S_V,S_M> m_xSeq_standard;
        S_V                               m_xSeq_standard_mins;
        S_V                               m_xSeq_standard_maxs;
        S_V                               m_xSeq_standard_ranges;
        std::vector<const S_V* >          m_xs_asterisks_standard; // to be deleted on destructor

        VectorSpace<P_V,P_M>       m_p_t_space;
        SequenceOfVectors<P_V,P_M> m_tSeq_original;
        P_V                               m_tSeq_mins;
        P_V                               m_tSeq_maxs;
        P_V                               m_tSeq_ranges;
        SequenceOfVectors<P_V,P_M> m_tSeq_standard;
        std::vector<const P_V* >          m_ts_asterisks_standard; // to be deleted on destructor

        VectorSpace<Q_V,Q_M>       m_n_eta_space;
        SequenceOfVectors<Q_V,Q_M> m_etaSeq_original;
        Q_V                               m_etaSeq_original_mean;
        Q_V                               m_etaSeq_original_std;
#ifdef _GPMSA_CODE_TREATS_SIMULATION_VECTORS_IN_CHUNKS
        std::vector<double>               m_etaSeq_chunkMeans;
        std::vector<double>               m_etaSeq_chunkStds;
#else
        double                            m_etaSeq_allMean;
        double                            m_etaSeq_allStd;
#endif
        SequenceOfVectors<Q_V,Q_M> m_etaSeq_transformed;
        Q_V                               m_etaSeq_transformed_mean;
        Q_V                               m_etaSeq_transformed_std;
        VectorSpace<Q_V,Q_M>       m_eta_space;
        Q_V                               m_etaVec_transformed;
        Q_M                               m_etaMat_transformed;

        VectorSpace<Q_V,Q_M>       m_m_space;
        Q_V                               m_m_unitVec;
        Q_M                               m_m_Imat;

        unsigned int                      m_paper_p_eta;
        VectorSpace<Q_V,Q_M>*      m_p_eta_space; // to be deleted on destructor
        Q_M*                              m_Kmat_eta;    // to be deleted on destructor
        std::vector<Q_V* >                m_kvec_is;     // to be deleted on destructor
        std::vector<Q_M* >                m_Kmat_is;     // to be deleted on destructor
        Q_M*                              m_Kmat;        // to be deleted on destructor
};

}  // End namespace QUESO

#endif  // UQ_SIMULATION_MODEL_H
