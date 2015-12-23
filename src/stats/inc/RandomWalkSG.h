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

#ifndef UQ_RANDOM_WALK_SG_H
#define UQ_RANDOM_WALK_SG_H

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

#include <queso/SequenceGenerator.h>

namespace QUESO {

template <class V = GslVector, class M = GslMatrix>
class RandomWalkSG : public SequenceGenerator<V, M>
{
public:
  RandomWalkSG(const char * prefix,
      const MhOptionsValues * options,
      const BaseVectorRV<V, M> & sourceRv,
      const V & initialPosition,
      const M * inputProposalCovMatrix);

  virtual ~RandomWalkSG();

  virtual void propose(unsigned int positionId,
      const BaseVectorSequence<V, M> & workingChain,
      V & proposedState);

  virtual void updateProposalCovariance(unsigned int positionId,
      const BaseVectorSequence<V, M> & workingChain,
      M & covarianceMatrix) = 0;
};

}  // End namespace QUESO

#endif  // UQ_RANDOM_WALK_SG_H
