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

#include <queso/asserts.h>
#include <queso/DistArray.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/ScalarSequence.h>
#include <queso/SampledScalarCdf.h>
#include <queso/OneDGrid.h>

#ifdef QUESO_HAS_TRILINOS
#include <Epetra_MpiComm.h>
#include <Epetra_SerialComm.h>
#endif

namespace QUESO {

// Constructor for a given inputMap and inputRowSize.
template<typename T>
DistArray<T>::DistArray(const Map& inputMap, const int inputRowSize)
  : m_Map(inputMap),
#ifdef QUESO_HAS_TRILINOS
  m_epetraDistArray(new EpetraExt::DistArray<T>(inputMap.epetraMap(),inputRowSize))
#else
  m_rowSize(inputRowSize)
#endif
{
#ifdef QUESO_HAS_TRILINOS
#else
  m_elements.resize(m_Map.NumGlobalElements());
  for (int i = 0; i < m_Map.NumGlobalElements(); ++i) {
    m_elements[i].resize(m_rowSize);
  }
#endif
}

// Destructor
template<typename T>
DistArray<T>::~DistArray()
{
#ifdef QUESO_HAS_TRILINOS
  delete m_epetraDistArray;
  m_epetraDistArray = NULL;
#else
  for (int i = 0; i < m_Map.NumGlobalElements(); ++i) {
    m_elements[i].clear();
  }
  m_elements.clear();
#endif
}

// Set methods
template<typename T>
DistArray<T>&
DistArray<T>::operator=(const DistArray<T>& rhs)
{
  return *this;
}

// Query methods
template<typename T>
T&
DistArray<T>::operator()(int localElementId, int colId)
{
#ifdef QUESO_HAS_TRILINOS
  return (*m_epetraDistArray)(localElementId,colId);
#else
  return m_elements[localElementId][colId];
#endif
}

template<typename T>
const T&
DistArray<T>::operator()(int localElementId, int colId) const
{
#ifdef QUESO_HAS_TRILINOS
  return (*m_epetraDistArray)(localElementId,colId);
#else
  return m_elements[localElementId][colId];
#endif
}

template<typename T>
int
DistArray<T>::GlobalLength() const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraDistArray->GlobalLength();
#else
  return m_Map.NumGlobalElements();
#endif
}

template<typename T>
int
DistArray<T>::MyLength() const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraDistArray->MyLength();
#else
  return m_Map.NumMyElements();
#endif
}

template<typename T>
int
DistArray<T>::RowSize() const
{
#ifdef QUESO_HAS_TRILINOS
  return m_epetraDistArray->RowSize();
#else
  return m_rowSize;
#endif
}

// I/O methods
template<typename T>
void
DistArray<T>::print(std::ostream& os) const
{
#ifdef QUESO_HAS_TRILINOS
  os << *m_epetraDistArray;
#else
  os << "m_rowSize = "           << m_rowSize
     << ", m_elements.size() = " << m_elements.size()
     << std::endl;
#endif

  return;
}

}  // End namespace QUESO

// Class prototypes with all the types that QUESO needs
template class QUESO::DistArray<QUESO::GslVector*>;
template class QUESO::DistArray<QUESO::GslMatrix*>;
template class QUESO::DistArray<QUESO::ScalarSequence<double>*>;
template class QUESO::DistArray<QUESO::SampledScalarCdf<double>*>;
template class QUESO::DistArray<QUESO::BaseOneDGrid<double>*>;
template class QUESO::DistArray<std::string>;
template class QUESO::DistArray<std::vector<double>*>;
