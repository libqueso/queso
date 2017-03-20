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

#ifndef UQ_TRANSITION_KERNEL_GROUP_H
#define UQ_TRANSITION_KERNEL_GROUP_H

#include <queso/VectorRV.h>
#include <queso/GaussianVectorRV.h>
#include <queso/ScalarFunctionSynchronizer.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*! \file TKGroup.h
 *  \brief Class to set up a Transition Kernel.
 *
 *  \class BaseTKGroup
 *  \brief This base class allows the representation of a transition kernel.
 *
 * The transition kernel is a conditional distribution function that represents the probability of
 * moving from a point \f$ x \in R^n \f$ to a point in the set \f$ A \in B \f$, where \f$ B \f$ is
 * a Borel set over \f$ R^n \f$. Since it is a distribution function, \f$ P(x,R^n)=1 \f$, where it
 * is permitted that the chain can make a transition from the point \f$ x \f$ to \f$ x \f$, that is
 * \f$ P(x, {x}) \f$ is not necessarily zero. */

template <class V = GslVector, class M = GslMatrix>
class BaseTKGroup {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  BaseTKGroup();

  //! Constructor.
  BaseTKGroup(const char*                    prefix,
                     const VectorSpace<V,M>& vectorSpace,
                     const std::vector<double>&     scales);
  //! Destructor.
  virtual ~BaseTKGroup();
  //@}

  //! @name Statistical/Mathematical methods
  //@{
  //! Whether or not the matrix is symmetric. See template specialization.
  virtual       bool                          symmetric                 () const = 0;

  //! Gaussian increment property to construct a transition kernel. See template specialization.
  virtual const BaseVectorRV<V,M>& rv                        (unsigned int                     stageId ) const = 0;

  //! Gaussian increment property to construct a transition kernel. See template specialization.
  virtual const BaseVectorRV<V,M>& rv                        (const std::vector<unsigned int>& stageIds) = 0;

  //! Constructs transition kernel pdf based on internal \c m_stageId variable
  virtual const BaseVectorRV<V, M> & rv(const V & position) const = 0;
  //@}

  //! @name Misc methods
  //@{
  //! QUESO's environment.
  const BaseEnvironment& env() const;

  //! Pre-computing position; access to protected attribute *m_preComputingPositions[stageId].
  const V&                                    preComputingPosition      (unsigned int stageId) const;

  //! Sets the pre-computing positions \c m_preComputingPositions[stageId] with a new vector of size \c position.
  virtual       bool                          setPreComputingPosition   (const V& position, unsigned int stageId);

  //! Clears the pre-computing positions \c m_preComputingPositions[stageId]
  virtual       void                          clearPreComputingPositions();

  //! Does nothing.  Subclasses may re-implement.  Returns the current stage id.
  virtual unsigned int set_dr_stage(unsigned int stageId);

  //! The user can override this method to implement state-dependent
  //! transition kernel behaviour.
  /*!
   * Default behaviour is a no-op.
   */
  virtual void updateTK() { };

  //! This method determines whether or not the user has 'dirtied' the
  //! covariance matrix, thereby necessitating an AM reset.
  /*!
   * Should return true if the cov matrix was modified in such a way to affect
   * the Adaptive Metropolis sample covariance calculation.
   */
  virtual bool covMatrixIsDirty() = 0;

  //! Performs whatever cleanup is needed (usually just resetting an internal
  //! bool flag) after the covariance matrix was modified.
  /*!
   * After this method is called, covMatrixIsDirty should return \c true.
   */
  virtual void cleanCovMatrix() = 0;
  //@}

  //! @name I/O methods
  //@{
  //! TODO: Prints the transition kernel.
  /*! \todo: implement me!*/
  virtual       void                          print                     (std::ostream& os) const;
  //@}
protected:
  const EmptyEnvironment*                m_emptyEnv;
  const BaseEnvironment&                 m_env;
        std::string                      m_prefix;
  const VectorSpace<V,M>*                m_vectorSpace;
        std::vector<double>              m_scales;
        std::vector<const V*>            m_preComputingPositions;
        std::vector<BaseVectorRV<V,M>* > m_rvs; // Gaussian, not Base... And nothing const...
  unsigned int m_stageId;
};

}  // End namespace QUESO

#endif // UQ_TRANSITION_KERNEL_GROUP_H
