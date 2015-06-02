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

#include <queso/ExperimentModel.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

template<class S_V,class S_M,class D_V,class D_M>
ExperimentModel<S_V,S_M,D_V,D_M>::ExperimentModel(
  const char*                                      prefix,
  const EmOptionsValues*                    alternativeOptionsValues, // dakota
  const ExperimentStorage<S_V,S_M,D_V,D_M>& experimentStorage,
  const std::vector<D_M* >&                        Dmats,
  const std::vector<D_M* >&                        Kmats_interp)
  :
  m_env                     (Dmats[0]->env()),
  m_optionsObj              (alternativeOptionsValues),
  m_paper_p_x               (experimentStorage.scenarioSpace().dimLocal()),
  m_paper_n                 (Dmats.size()),
  m_paper_p_delta           (Dmats[0]->numCols()),
  m_paper_n_y               (0),
  m_Dmats                   (Dmats),
  m_Kmats_interp            (Kmats_interp),
  m_n_y_space               (NULL),
  m_Dmat_BlockDiag          (NULL)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering ExperimentModel<S_V,S_M,D_V,D_M>::constructor()"
                            << ": prefix = " << prefix
                            << ", alternativeOptionsValues = " << alternativeOptionsValues
                            << ", m_env.optionsInputFileName() = " << m_env.optionsInputFileName()
                            << "\n  m_paper_p_x = "                << m_paper_p_x
                            << "\n  m_paper_n = "                  << m_paper_n
                            << "\n  m_paper_p_delta = "            << m_paper_p_delta
                            << "\n  m_paper_n_y = "                << "to be set soon"
                            << std::endl;
  }

  // If NULL, we create one
  if (m_optionsObj == NULL) {
    EmOptionsValues * tempOptions = new EmOptionsValues(&m_env, prefix);

    // We did this dance because scanOptionsValues is not a const method, but
    // m_optionsObj is a pointer to const
    m_optionsObj = tempOptions;
  }

  // We'll need to remove this later
  m_experimentModelOptions = new ExperimentModelOptions(m_env, prefix);

  queso_require_greater_equal_msg(m_optionsObj->m_Gvalues.size(), 1, "invalid m_Gs");

  queso_require_equal_to_msg(m_paper_n, experimentStorage.xs_standard().size(), "invalid m_paper_n");

  unsigned int sumGs = 0;
  for (unsigned int i = 0; i < m_optionsObj->m_Gvalues.size(); ++i) {
    sumGs += m_optionsObj->m_Gvalues[i];
  }
  queso_require_equal_to_msg(m_paper_p_delta, sumGs, "inconsistent input");

  //***********************************************************************
  // Form 'Dmat_BlockDiag' matrix
  //***********************************************************************
  for (unsigned int i = 0; i < m_Dmats.size(); ++i) {
    queso_require_equal_to_msg(m_Dmats[i]->numCols(), m_paper_p_delta, "inconsistent m_Dmats");
    m_paper_n_y += m_Dmats[i]->numRowsLocal();
  }
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In ExperimentModel<S_V,S_M,D_V,D_M>::constructor()"
                            << ": prefix = " << prefix
                            << "\n  m_paper_n_y = " << m_paper_n_y
                            << std::endl;
  }
  m_n_y_space = new VectorSpace<D_V,D_M>(m_env, "n_y_", m_paper_n_y, NULL),
  m_Dmat_BlockDiag = new D_M(m_env,m_n_y_space->map(),m_paper_n*m_paper_p_delta);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "In ExperimentModel<S_V,S_M,D_V,D_M>::constructor()"
                            << ": key-debug"
                            << ", m_Dmat_BlockDiag just created (not yet populated)"
                            << ", numRowsLocal = " << m_Dmat_BlockDiag->numRowsLocal()
                            << ", numCols = "      << m_Dmat_BlockDiag->numCols()
                            << std::endl;
  }
  m_Dmat_BlockDiag->fillWithBlocksDiagonally(0,0,m_Dmats,true,true);

  std::set<unsigned int> tmpSet;
  tmpSet.insert(m_env.subId());
  if (false) { //gcmOptionsObj.m_dataOutputAllowedSet.find(m_env.subId()) != gcmOptionsObj.m_dataOutputAllowedSet.end()) {
    m_Dmat_BlockDiag->subWriteContents("Dmat_BlockDiag",
                                       "Dmat_BlockDiag",
                                       "m",
                                       tmpSet);
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving ExperimentModel<S_V,S_M,D_V,D_M>::constructor()"
                            << ": prefix = " << m_optionsObj->m_prefix
                            << std::endl;
  }
}

template<class S_V,class S_M,class D_V,class D_M>
ExperimentModel<S_V,S_M,D_V,D_M>::~ExperimentModel()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering ExperimentModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::destructor()..."
                            << std::endl;
  }

  if (m_Dmat_BlockDiag) delete m_Dmat_BlockDiag; // to be deleted on destructor
  if (m_n_y_space     ) delete m_n_y_space;      // to be deleted on destructor
  if (m_optionsObj    ) delete m_optionsObj;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving ExperimentModel<S_V,S_M,P_V,P_M,Q_V,Q_M>::destructor()"
                            << std::endl;
  }

  if (m_experimentModelOptions) {
    delete m_experimentModelOptions;
  }
}

template<class S_V,class S_M,class D_V,class D_M>
unsigned int
ExperimentModel<S_V,S_M,D_V,D_M>::numBasis() const
{
  return m_paper_p_delta;
}

template<class S_V,class S_M,class D_V,class D_M>
unsigned int
ExperimentModel<S_V,S_M,D_V,D_M>::numBasisGroups() const
{
  return m_optionsObj->m_Gvalues.size();
}

template<class S_V,class S_M,class D_V,class D_M>
const std::vector<unsigned int>&
ExperimentModel<S_V,S_M,D_V,D_M>::Gs() const
{
  return m_optionsObj->m_Gvalues;
}

template<class S_V,class S_M,class D_V,class D_M>
const std::vector<D_M* >&
ExperimentModel<S_V,S_M,D_V,D_M>::Kmats_interp() const
{
  return m_Kmats_interp;
}

template<class S_V,class S_M,class D_V,class D_M>
const D_M&
ExperimentModel<S_V,S_M,D_V,D_M>::Dmat(unsigned int basisId) const
{
  return *(m_Dmats[basisId]);
}

template<class S_V,class S_M,class D_V,class D_M>
const D_M&
ExperimentModel<S_V,S_M,D_V,D_M>::Dmat_BlockDiag() const
{
  return *m_Dmat_BlockDiag;
}


template<class S_V,class S_M,class D_V,class D_M>
const ExperimentModelOptions&
ExperimentModel<S_V,S_M,D_V,D_M>::optionsObj() const
{
  // This methods returns the old ExperimentModelOptions object.  We're using
  // the new EmOptionsValues so this method is deprecated.
  queso_deprecated();
  return *m_experimentModelOptions;
}

template<class S_V,class S_M,class D_V,class D_M>
void
ExperimentModel<S_V,S_M,D_V,D_M>::print(std::ostream& os) const
{
  return;
}

}  // End namespace QUESO

template class QUESO::ExperimentModel<QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix>;
