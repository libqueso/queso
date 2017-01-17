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

#include <queso/Defines.h>
#include <queso/2dArrayOfStuff.h>

namespace QUESO {

// Default constructor -------------------------------------------------
template <class T>
TwoDArray<T>::TwoDArray(unsigned int numRows, unsigned int numCols)
  : m_numRows(numRows),
    m_numCols(numCols),
    m_data   (m_numRows,NULL)
{
  for (unsigned int i = 0; i < m_numRows; ++i) {
    m_data[i] = new std::vector<T*>(m_numCols,NULL);
  }
}

// Destructor ----------------------------------------------------------
template <class T>
TwoDArray<T>::~TwoDArray()
{
  for (unsigned int i = 0; i < m_numRows; ++i) {
    for (unsigned int j = 0; j < m_numCols; ++j) {
      if ((*(m_data[i]))[j] != NULL) delete (*(m_data[i]))[j];
    }
    delete m_data[i];
  }
}

// Property methods-----------------------------------------------------
template <class T>
unsigned int TwoDArray<T>::numRows() const
{
  return m_numRows;
}

template <class T>
unsigned int TwoDArray<T>::numCols() const
{
  return m_numCols;
}

template <class T>
void TwoDArray<T>::setLocation(unsigned int i, unsigned int j, T* info)
{
  queso_require_msg(!((i >= m_numRows) || (j >= m_numCols) || (m_data[i] == NULL)), "invalid situation");
  (*(m_data[i]))[j] = info;
  return;
}

// Accessor methods-----------------------------------------------------
template <class T>
T& TwoDArray<T>::operator()(unsigned int i, unsigned int j)
{
  // FIXME - Should this really be active in optimized codes? - RHS
  queso_require_less(i, m_numRows);
  queso_require_less(j, m_numCols);
  queso_require(m_data[i]);
  queso_require((*m_data[i])[j]);

  return *(*(m_data[i]))[j];
}

template <class T>
const T& TwoDArray<T>::operator()(unsigned int i, unsigned int j) const
{
  // FIXME - Should this really be active in optimized codes? - RHS
  queso_require_less(i, m_numRows);
  queso_require_less(j, m_numCols);
  queso_require(m_data[i]);
  queso_require((*m_data[i])[j]);

  return *(*(m_data[i]))[j];
}

}  // End namespace QUESO
