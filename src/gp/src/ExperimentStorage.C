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

#include <queso/ExperimentStorage.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

template<class S_V,class S_M,class D_V,class D_M>
ExperimentStorage<S_V,S_M,D_V,D_M>::ExperimentStorage(
  const VectorSpace<S_V,S_M>& scenarioSpace,
  unsigned int                numExperiments)
  :
  m_env                    (scenarioSpace.env()),
  m_scenarioSpace          (scenarioSpace),
  m_paper_n                (numExperiments),
  m_paper_n_ys_transformed (m_paper_n,0),
  m_paper_n_y              (0),
  m_addId                  (0),
  m_scenarioVecs_standard  (m_paper_n, (S_V*) NULL),
  m_dataVecs_transformed   (m_paper_n, (D_V*) NULL),
  m_covMats_transformed_inv(m_paper_n, (D_M*) NULL),
  m_y_space                (NULL),
  m_yVec_transformed       (NULL),
  m_Wy                     (NULL)
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

  delete m_Wy;
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
  const D_M& covMat_transformed_inv)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering ExperimentStorage<S_V,S_M,D_V,D_M>::addExperiment()"
                            << ": m_addId = " << m_addId
                            << "\n  scenarioVec_standard = "   << scenarioVec_standard
                            << "\n  dataVec_transformed = "    << dataVec_transformed
                            << "\n  covMat_transformed_inv = " << covMat_transformed_inv
                            << std::endl;
  }

  queso_require_less_msg(m_addId, m_paper_n, "too many adds...");

  m_scenarioVecs_standard  [m_addId] = &scenarioVec_standard;
  m_dataVecs_transformed   [m_addId] = &dataVec_transformed;
  m_covMats_transformed_inv[m_addId] = &covMat_transformed_inv;
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
    // Form 'yVec_transformed', 'Wy' matrix, and compute its inverse
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

    m_Wy = new D_M(m_y_space->zeroVector());
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
      *m_env.subDisplayFile() << "In ExperimentStorage<S_V,S_M,D_V,D_M>::addExperiment()"
                              << ": key-debug"
                              << ", m_Wy just created (not yet populated)"
                              << ", numRowsLocal = " << m_Wy->numRowsLocal()
                              << ", numCols = "      << m_Wy->numCols()
                              << std::endl;
    }

    m_Wy->fillWithBlocksDiagonally(0,0,m_covMats_transformed_inv,true,true);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
      *m_env.subDisplayFile() << "In ExperimentStorage<S_V,S_M,D_V,D_M>::addExperiment()"
                              << ": key-debug"
                              << ", m_Wy just populated"
                              << std::endl;
    }

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
      *m_env.subDisplayFile() << "In ExperimentStorage<S_V,S_M,D_V,D_M>::addExperiment()"
                              << ": m_Wy->lnDeterminant() = " << m_Wy->lnDeterminant()
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
  queso_require_less_msg(experimentId, m_scenarioVecs_standard.size(), "experimentId is too large");

  queso_require_msg(m_scenarioVecs_standard[experimentId], "vector is NULL");

  return *(m_scenarioVecs_standard[experimentId]);
}

template<class S_V,class S_M,class D_V,class D_M>
const D_V&
ExperimentStorage<S_V,S_M,D_V,D_M>::dataVec_transformed(unsigned int experimentId) const
{
  queso_require_less_msg(experimentId, m_dataVecs_transformed.size(), "experimentId is too large");

  queso_require_msg(m_dataVecs_transformed[experimentId], "vector is NULL");

  return *(m_dataVecs_transformed[experimentId]);
}

template<class S_V,class S_M,class D_V,class D_M>
const D_V&
ExperimentStorage<S_V,S_M,D_V,D_M>::yVec_transformed() const
{
  queso_require_msg(m_yVec_transformed, "'m_yVec_transformed' is NULL");

  return *m_yVec_transformed;
}

template<class S_V,class S_M,class D_V,class D_M>
const D_M&
ExperimentStorage<S_V,S_M,D_V,D_M>::Wy() const
{
  queso_require_msg(m_Wy, "'m_Wy' is NULL");

  return *m_Wy;
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

}  // End namespace QUESO

template class QUESO::ExperimentStorage<QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix>;
