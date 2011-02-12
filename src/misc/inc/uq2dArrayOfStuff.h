//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011 The PECOS Development Team
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

#ifndef __UQ_2D_ARRAY_OF_STUFF_H__
#define __UQ_2D_ARRAY_OF_STUFF_H__

template <class T>
class uq2dArrayOfStuff
{
public:
  uq2dArrayOfStuff(unsigned int numRows, unsigned int numCols);
 ~uq2dArrayOfStuff();

        unsigned int numRows    ()                                 const;
        unsigned int numCols    ()                                 const;
        void         setLocation(unsigned int i, unsigned int j, T* info);
        T&           operator   ()(unsigned int i, unsigned int j);
  const T&           operator   ()(unsigned int i, unsigned int j) const;

private:
  unsigned int m_numRows;
  unsigned int m_numCols;

  std::vector<std::vector<T*>* > m_data;
};

template <class T>
uq2dArrayOfStuff<T>::uq2dArrayOfStuff(
  unsigned int numRows,
  unsigned int numCols)
  :
  m_numRows(numRows),
  m_numCols(numCols),
  m_data   (m_numRows,NULL)
{
  for (unsigned int i = 0; i < m_numRows; ++i) {
    m_data[i] = new std::vector<T*>(m_numCols,NULL);
  }
}

template <class T>
uq2dArrayOfStuff<T>::~uq2dArrayOfStuff()
{
  for (unsigned int i = 0; i < m_numRows; ++i) {
    for (unsigned int j = 0; j < m_numCols; ++j) {
      if ((*(m_data[i]))[j] != NULL) delete (*(m_data[i]))[j];
    }
    delete m_data[i];
  }
}

template <class T>
unsigned int
uq2dArrayOfStuff<T>::numRows() const
{
  return m_numRows;
}

template <class T>
unsigned int
uq2dArrayOfStuff<T>::numCols() const
{
  return m_numCols;
}

template <class T>
void
uq2dArrayOfStuff<T>::setLocation(unsigned int i, unsigned int j, T* info)
{
  UQ_FATAL_TEST_MACRO((i >= m_numRows) || (j >= m_numCols) || (m_data[i] == NULL),
                      UQ_UNAVAILABLE_RANK,
                      "uq2dArrayOfStuff<T>::setLocation()",
                      "invalid situation");
  (*(m_data[i]))[j] = info;
  return;
}

template <class T>
T&
uq2dArrayOfStuff<T>::operator()(unsigned int i, unsigned int j)
{
  UQ_FATAL_TEST_MACRO((i >= m_numRows) || (j >= m_numCols) || (m_data[i] == NULL) || ((*m_data[i])[j] == NULL),
                      UQ_UNAVAILABLE_RANK,
                      "uq2dArrayOfStuff<T>::operator(1)",
                      "invalid situation");
  return *(*(m_data[i]))[j];
}

template <class T>
const T&
uq2dArrayOfStuff<T>::operator()(unsigned int i, unsigned int j) const
{
  UQ_FATAL_TEST_MACRO((i >= m_numRows) || (j >= m_numCols) || (m_data[i] == NULL) || ((*m_data[i])[j] == NULL),
                      UQ_UNAVAILABLE_RANK,
                      "uq2dArrayOfStuff<T>::operator(2)",
                      "invalid situation");
  return *(*(m_data[i]))[j];
}
#endif // __UQ_2D_ARRAY_OF_STUFF_H__

