//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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
// 
// $Id$
//
//--------------------------------------------------------------------------

#ifndef __UQ_EXPERIMENT_STORAGE_H__
#define __UQ_EXPERIMENT_STORAGE_H__

#include <queso/VectorSpace.h>

namespace QUESO {

template <class S_V,class S_M,class D_V,class D_M>
class ExperimentStorage
{
public:
  ExperimentStorage(const VectorSpace<S_V,S_M>& scenarioSpace, unsigned int numExperiments);
 ~ExperimentStorage();

        void                         addExperiment         (const S_V& scenarioVec_standard, const D_V& dataVec_transformed, const D_M& Wmat_transformed);
        unsigned int                 numExperiments        () const;
  const VectorSpace<S_V,S_M>& scenarioSpace         () const;
  const std::vector<const S_V* >&    xs_standard           () const;
  const std::vector<unsigned int>&   n_ys_transformed      () const;
        unsigned int                 n_y                   () const;
  const S_V&                         scenarioVec_standard  (unsigned int experimentId) const;
  const D_V&                         dataVec_transformed   (unsigned int experimentId) const;
  const D_V&                         yVec_transformed      () const;
  const D_M&                         Wmat_transformed      (unsigned int experimentId) const;
  const D_M&                         Wmat_transformed_y    () const;
  const D_M&                         Wmat_transformed_y_inv() const;

  const BaseEnvironment&      env                   () const;
        void                         print                 (std::ostream& os) const;

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
        std::vector<const D_M* >     m_Wmats_transformed;
        VectorSpace<D_V,D_M>* m_y_space;
        D_V*                         m_yVec_transformed;
        D_M*                         m_Wmat_transformed_y;
        D_M*                         m_Wmat_transformed_y_inv;
};

template<class S_V,class S_M,class D_V,class D_M>
std::ostream& operator<<(std::ostream& os, const ExperimentStorage<S_V,S_M,D_V,D_M>& obj);

template<class S_V,class S_M,class D_V,class D_M>
ExperimentStorage<S_V,S_M,D_V,D_M>::ExperimentStorage(
  const VectorSpace<S_V,S_M>& scenarioSpace,
  unsigned int                       numExperiments)
  :
  m_env                   (scenarioSpace.env()),
  m_scenarioSpace         (scenarioSpace),
  m_paper_n               (numExperiments),
  m_paper_n_ys_transformed(m_paper_n,0),
  m_paper_n_y             (0),
  m_addId                 (0),
  m_scenarioVecs_standard (m_paper_n, (S_V*) NULL),
  m_dataVecs_transformed  (m_paper_n, (D_V*) NULL),
  m_Wmats_transformed     (m_paper_n, (D_M*) NULL),
  m_y_space               (NULL),
  m_yVec_transformed      (NULL),
  m_Wmat_transformed_y    (NULL),
  m_Wmat_transformed_y_inv(NULL)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering ExperimentStorage<S_V,S_M,D_V,D_M>::constructor()"
                            << "\n  m_paper_n = " << m_paper_n
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving ExperimentStorage<S_V,S_M,D_V,D_M>::constructor()"
                            << std::endl;
  }
}

template<class S_V,class S_M,class D_V,class D_M>
ExperimentStorage<S_V,S_M,D_V,D_M>::~ExperimentStorage()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering ExperimentStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::destructor()..."
                            << std::endl;
  }

  delete m_Wmat_transformed_y_inv;
  delete m_Wmat_transformed_y;
  delete m_yVec_transformed;
  delete m_y_space;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving ExperimentStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::destructor()"
                            << std::endl;
  }
}

template<class S_V,class S_M,class D_V,class D_M>
void
ExperimentStorage<S_V,S_M,D_V,D_M>::addExperiment(
  const S_V& scenarioVec_standard,
  const D_V& dataVec_transformed,
  const D_M& Wmat_transformed)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering ExperimentStorage<S_V,S_M,D_V,D_M>::addExperiment()"
                            << ": m_addId = " << m_addId
                            << "\n  scenarioVec_standard = " << scenarioVec_standard
                            << "\n  dataVec_transformed = "  << dataVec_transformed
                            << "\n  WmatVec_transformed = "  << Wmat_transformed
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(m_addId >= m_paper_n,
                      m_env.worldRank(),
                      "ExperimentStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::addExperiment()",
                      "too many adds...");

  m_scenarioVecs_standard [m_addId] = &scenarioVec_standard;
  m_dataVecs_transformed  [m_addId] = &dataVec_transformed;
  m_Wmats_transformed     [m_addId] = &Wmat_transformed;
  m_paper_n_ys_transformed[m_addId] = dataVec_transformed.sizeLocal();
  m_paper_n_y += m_paper_n_ys_transformed[m_addId];
  m_addId++;

  if (m_addId == m_paper_n) {
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
      *m_env.subDisplayFile() << "KEY In ExperimentStorage<S_V,S_M,D_V,D_M>::addExperiment()"
                              << ": m_addId = " << m_addId
                              << ", m_paper_n_y = " << m_paper_n_y
                              << std::endl;
    }

    //***********************************************************************
    // Form 'yVec_transformed', 'Wmat_transformed_y' matrix, and compute its inverse
    //***********************************************************************
    m_y_space = new VectorSpace<D_V,D_M>(m_env, "m_y_exp_storage", m_paper_n_y, NULL),
    m_yVec_transformed = new D_V(m_y_space->zeroVector());
    m_yVec_transformed->cwSetConcatenated(m_dataVecs_transformed);

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In ExperimentStorage<S_V,S_M,D_V,D_M>::addExperiment()"
                              << ": populated yVec_transformed of size = " << m_yVec_transformed->sizeLocal()
                              << "\n *m_yVec_transformed = " << *m_yVec_transformed
                              << std::endl;
    }

    m_Wmat_transformed_y = new D_M(m_y_space->zeroVector());
    m_Wmat_transformed_y->fillWithBlocksDiagonally(0,0,m_Wmats_transformed,true,true);
    m_Wmat_transformed_y_inv = new D_M(m_Wmat_transformed_y->inverse()); // inversion savings

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In ExperimentStorage<S_V,S_M,D_V,D_M>::addExperiment()"
                              << ": m_Wmat_transformed_y->lnDeterminant() = "     << m_Wmat_transformed_y->lnDeterminant()
                              << ", m_Wmat_transformed_y_inv->lnDeterminant() = " << m_Wmat_transformed_y_inv->lnDeterminant()
                              << std::endl;
    }
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving ExperimentStorage<S_V,S_M,D_V,D_M>::addExperiment()"
                            << ": m_addId = " << m_addId
                            << std::endl;
  }

  return;
}

template<class S_V,class S_M,class D_V,class D_M>
const std::vector<const S_V* >&
ExperimentStorage<S_V,S_M,D_V,D_M>::xs_standard() const
{
  return m_scenarioVecs_standard;
}

template<class S_V,class S_M,class D_V,class D_M>
unsigned int
ExperimentStorage<S_V,S_M,D_V,D_M>::numExperiments() const
{
  return m_paper_n;
}

template<class S_V,class S_M,class D_V,class D_M>
const VectorSpace<S_V,S_M>&
ExperimentStorage<S_V,S_M,D_V,D_M>::scenarioSpace() const
{
  return m_scenarioSpace;
}

template<class S_V,class S_M,class D_V,class D_M>
const std::vector<unsigned int>&
ExperimentStorage<S_V,S_M,D_V,D_M>::n_ys_transformed() const
{
  return m_paper_n_ys_transformed;
}

template<class S_V,class S_M,class D_V,class D_M>
unsigned int
ExperimentStorage<S_V,S_M,D_V,D_M>::n_y() const
{
  return m_paper_n_y;
}

template<class S_V,class S_M,class D_V,class D_M>
const S_V&
ExperimentStorage<S_V,S_M,D_V,D_M>::scenarioVec_standard(unsigned int experimentId) const
{
  UQ_FATAL_TEST_MACRO(experimentId >= m_scenarioVecs_standard.size(),
                      m_env.worldRank(),
                      "ExperimentStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::scenarioVec_standard()",
                      "experimentId is too large");

  UQ_FATAL_TEST_MACRO(m_scenarioVecs_standard[experimentId] == NULL,
                      m_env.worldRank(),
                      "ExperimentStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::scenarioVec_standard()",
                      "vector is NULL");

  return *(m_scenarioVecs_standard[experimentId]);
}

template<class S_V,class S_M,class D_V,class D_M>
const D_V&
ExperimentStorage<S_V,S_M,D_V,D_M>::dataVec_transformed(unsigned int experimentId) const
{
  UQ_FATAL_TEST_MACRO(experimentId >= m_dataVecs_transformed.size(),
                      m_env.worldRank(),
                      "ExperimentStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::dataVec_transformed()",
                      "experimentId is too large");

  UQ_FATAL_TEST_MACRO(m_dataVecs_transformed[experimentId] == NULL,
                      m_env.worldRank(),
                      "ExperimentStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::dataVec_transformed()",
                      "vector is NULL");

  return *(m_dataVecs_transformed[experimentId]);
}

template<class S_V,class S_M,class D_V,class D_M>
const D_V&
ExperimentStorage<S_V,S_M,D_V,D_M>::yVec_transformed() const
{
  UQ_FATAL_TEST_MACRO(m_yVec_transformed == NULL,
                      m_env.worldRank(),
                      "ExprimentStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::yVec_transformed()",
                      "'m_yVec_transformed' is NULL");

  return *m_yVec_transformed;
}

template<class S_V,class S_M,class D_V,class D_M>
const D_M&
ExperimentStorage<S_V,S_M,D_V,D_M>::Wmat_transformed(unsigned int experimentId) const
{
  UQ_FATAL_TEST_MACRO(experimentId >= m_Wmats_transformed.size(),
                      m_env.worldRank(),
                      "ExperimentStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::Wmat_transformed()",
                      "experimentId is too large");

  UQ_FATAL_TEST_MACRO(m_Wmats_transformed[experimentId] == NULL,
                      m_env.worldRank(),
                      "ExperimentStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::Wmat_transformed()",
                      "matrix is NULL");

  return *(m_Wmats_transformed[experimentId]);
}

template<class S_V,class S_M,class D_V,class D_M>
const D_M&
ExperimentStorage<S_V,S_M,D_V,D_M>::Wmat_transformed_y() const
{
  UQ_FATAL_TEST_MACRO(m_Wmat_transformed_y == NULL,
                      m_env.worldRank(),
                      "ExprimentStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::Wmat_transformed_y()",
                      "'m_Wmat_transformed_y' is NULL");

  return *m_Wmat_transformed_y;
}

template<class S_V,class S_M,class D_V,class D_M>
const D_M&
ExperimentStorage<S_V,S_M,D_V,D_M>::Wmat_transformed_y_inv() const
{
  UQ_FATAL_TEST_MACRO(m_Wmat_transformed_y_inv == NULL,
                      m_env.worldRank(),
                      "ExprimentStorage<S_V,S_M,P_V,P_M,Q_V,Q_M>::Wmat_transformed_y_inv()",
                      "'m_Wmat_transformed_y_inv' is NULL");

  return *m_Wmat_transformed_y_inv;
}

template<class S_V,class S_M,class D_V,class D_M>
const BaseEnvironment&
ExperimentStorage<S_V,S_M,D_V,D_M>::env() const
{
  return m_env;
}

template<class S_V,class S_M,class D_V,class D_M>
void
ExperimentStorage<S_V,S_M,D_V,D_M>::print(std::ostream& os) const
{
  return;
}

template<class S_V,class S_M,class D_V,class D_M>
std::ostream& operator<<(std::ostream& os, const ExperimentStorage<S_V,S_M,D_V,D_M>& obj)
{
  obj.print(os);

  return os;
}

}  // End namespace QUESO

#endif // __UQ_EXPERIMENT_STORAGE_H__
