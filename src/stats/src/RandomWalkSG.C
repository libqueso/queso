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

#include <queso/RandomWalkSG.h>
#include <queso/ScaledCovMatrixTKGroup.h>

namespace QUESO {

template <class V, class M>
RandomWalkSG<V, M>::RandomWalkSG(const char * prefix,
    const MhOptionsValues * options,
    const BaseVectorRV<V, M> & sourceRv,
    const V & initialPosition,
    const M * inputProposalCovMatrix)
: SequenceGenerator<V, M>(prefix,
                          options,
                          sourceRv,
                          initialPosition,
                          inputProposalCovMatrix)
{
  std::vector<double> scales(1, 1.0);  // One element equal to 1.
  this->m_tk = new ScaledCovMatrixTKGroup<V, M>(
      "", sourceRv.imageSet().vectorSpace(), scales, *inputProposalCovMatrix);
}

template <class V, class M>
RandomWalkSG<V, M>::~RandomWalkSG()
{
}

template <class V, class M>
void
RandomWalkSG<V, M>::propose(unsigned int positionId,
    const BaseVectorSequence<V, M> & workingChain,
    V & proposedState)
{
  // Get current state
  workingChain.getPositionValues(positionId, proposedState);

  // Get any user-provided covariance matrix
  M covMatrix(this->m_vectorSpace.zeroVector());
  this->updateProposalCovariance(positionId, workingChain, covMatrix);

  (dynamic_cast<ScaledCovMatrixTKGroup<V, M>* >(this->m_tk))->updateLawCovMatrix(covMatrix);

  // This updates the mean of the transition kernel
  this->m_tk->setPreComputingPosition(proposedState, 0);

  // Get proposed state
  this->m_tk->rv(0).realizer().realization(proposedState);
}

}  // End namespace QUESO

template class QUESO::RandomWalkSG<QUESO::GslVector, QUESO::GslMatrix>;
