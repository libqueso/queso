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

#ifndef UQ_EXPERIMENT_STORAGE_H
#define UQ_EXPERIMENT_STORAGE_H

#include <queso/VectorSpace.h>

namespace QUESO {

class GslVector;
class GslMatrix;

template <class S_V = GslVector, class S_M = GslMatrix, class D_V = GslVector, class D_M = GslMatrix>
class ExperimentStorage
{
public:
  ExperimentStorage(const VectorSpace<S_V,S_M>& scenarioSpace, unsigned int numExperiments);
 ~ExperimentStorage();

        void                         addExperiment       (const S_V& scenarioVec_standard, const D_V& dataVec_transformed, const D_M& covMat_transformed_inv);
        unsigned int                 numExperiments      () const;
  const VectorSpace<S_V,S_M>& scenarioSpace              () const;
  const std::vector<const S_V* >&    xs_standard         () const;
  const std::vector<unsigned int>&   n_ys_transformed    () const;
        unsigned int                 n_y                 () const;
  const S_V&                         scenarioVec_standard(unsigned int experimentId) const;
  const D_V&                         dataVec_transformed (unsigned int experimentId) const;
  const D_V&                         yVec_transformed    () const;
  const D_M&                         Wy                  () const;

  const BaseEnvironment&             env                 () const;
        void                         print               (std::ostream& os) const;
  friend std::ostream& operator<<(std::ostream& os,
      const ExperimentStorage<S_V,S_M,D_V,D_M>& obj) {
    obj.print(os);
    return os;
  }


private:
  // Private variables
  const BaseEnvironment&      m_env;
  const VectorSpace<S_V,S_M>& m_scenarioSpace;
        unsigned int                 m_paper_n;
        std::vector<unsigned int>    m_paper_n_ys_transformed;
        unsigned int                 m_paper_n_y;

        unsigned int                 m_addId;
        std::vector<const S_V* >     m_scenarioVecs_standard;
        std::vector<const D_V* >     m_dataVecs_transformed;
        std::vector<const D_M* >     m_covMats_transformed_inv;  // = W_i's
        VectorSpace<D_V,D_M>* m_y_space;
        D_V*                         m_yVec_transformed;
        D_M*                         m_Wy;
};

}  // End namespace QUESO

#endif // UQ_EXPERIMENT_STORAGE_H
