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

#ifndef UQ_HESSIANCOV_TRANSITION_KERNEL_GROUP_H
#define UQ_HESSIANCOV_TRANSITION_KERNEL_GROUP_H

#include <queso/TKGroup.h>
#include <queso/VectorRV.h>
#include <queso/ScalarFunctionSynchronizer.h>

namespace QUESO {

class GslVector;
class GslMatrix;

//*****************************************************
// TK with Hessians
//*****************************************************
/*! \class HessianCovMatricesTKGroup
 *  \brief This class allows the representation of a transition kernel with Hessians. */

template <class V = GslVector, class M = GslMatrix>
class HessianCovMatricesTKGroup : public BaseTKGroup<V,M> {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  HessianCovMatricesTKGroup(const char*                                   prefix,
                                   const VectorSpace<V,M>&                vectorSpace,
                                   const std::vector<double>&                    scales,
                                   const ScalarFunctionSynchronizer<V,M>& targetPdfSynchronizer);

   //! Destructor.
  ~HessianCovMatricesTKGroup();
  //@}

  //! @name Statistical/Mathematical methods
  //@{
  //! Whether or not the matrix is symmetric. Always 'false'.
  /*! \todo: It only returns 'false', thus a test for its symmetricity must be included.*/
  bool                          symmetric                 () const;

  //! Gaussian increment property to construct a transition kernel.
  const GaussianVectorRV<V,M>& rv                        (unsigned int                     stageId ) const;

  //! Gaussian increment property to construct a transition kernel.
  const GaussianVectorRV<V,M>& rv                        (const std::vector<unsigned int>& stageIds);

  virtual const GaussianVectorRV<V, M> & rv(const V & position) const;
   //@}

  //! @name Misc methods
  //@{
  //! Sets the pre-computing positions \c m_preComputingPositions[stageId] with a new vector of size \c position.
  bool                          setPreComputingPosition   (const V& position, unsigned int  stageId );

  //! Clears the pre-computing positions \c m_preComputingPositions[stageId]
  void                          clearPreComputingPositions();

  virtual unsigned int set_dr_stage(unsigned int stageId);

  virtual bool covMatrixIsDirty() { return false; }
  virtual void cleanCovMatrix() { }
  //@}

  //! @name I/O methods
  //@{
  //! TODO: Prints the transition kernel.
  /*! \todo: implement me!*/
  void                          print                     (std::ostream& os) const;
  //@}
private:
  using BaseTKGroup<V,M>::m_env;
  using BaseTKGroup<V,M>::m_prefix;
  using BaseTKGroup<V,M>::m_vectorSpace;
  using BaseTKGroup<V,M>::m_scales;
  using BaseTKGroup<V,M>::m_preComputingPositions;
  using BaseTKGroup<V,M>::m_rvs;

  const ScalarFunctionSynchronizer<V,M>& m_targetPdfSynchronizer;
  std::vector<V*>                               m_originalNewtonSteps;
  std::vector<M*>                               m_originalCovMatrices;
};

}  // End namespace QUESO

#endif // UQ_HESSIANCOV_TRANSITION_KERNEL_GROUP_H
