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

#ifndef UQ_EXPERIMENT_MODEL_H
#define UQ_EXPERIMENT_MODEL_H

#include <queso/ExperimentModelOptions.h>
#include <queso/ExperimentStorage.h>
#include <queso/SequenceOfVectors.h>
#include <queso/Environment.h>

namespace QUESO {

class GslVector;
class GslMatrix;

template <class S_V = GslVector, class S_M = GslMatrix, class D_V = GslVector, class D_M = GslMatrix>
class ExperimentModel
{
public:
  ExperimentModel(const char*                                      prefix,
                         const EmOptionsValues*                    alternativeOptionsValues, // dakota
                         const ExperimentStorage<S_V,S_M,D_V,D_M>& experimentStorage,
                         const std::vector<D_M* >&                        Dmats,
                         const std::vector<D_M* >&                        Kmats_interp);
 ~ExperimentModel();

        unsigned int                   numBasis      () const;
        unsigned int                   numBasisGroups() const;
  const std::vector<unsigned int>&     Gs            () const;
  const D_M&                           Dmat          (unsigned int basisId) const;
  const D_M&                           Dmat_BlockDiag() const;
  const std::vector<D_M* >&            Kmats_interp  () const;

  const ExperimentModelOptions& optionsObj    () const;
        void                           print         (std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os,
      const ExperimentModel<S_V,S_M,D_V,D_M>& obj) {
    obj.print(os);
    return os;
  }


private:
  // Private variables
  const BaseEnvironment & m_env;
  const EmOptionsValues * m_optionsObj;
  ExperimentModelOptions * m_experimentModelOptions;

  unsigned int m_paper_p_x;
  unsigned int m_paper_n;
  unsigned int m_paper_p_delta;
  unsigned int m_paper_n_y;
  std::vector<D_M* > m_Dmats;          // NOT to be deleted on destructor
  std::vector<D_M* > m_Kmats_interp;   // NOT to be deleted on destructor
  VectorSpace<D_V,D_M>* m_n_y_space;      // to be deleted on destructor
  D_M* m_Dmat_BlockDiag; // to be deleted on destructor
};

}  // End namespace QUESO

#endif // UQ_EXPERIMENT_MODEL_H
