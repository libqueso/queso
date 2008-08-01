/* uq/libs/basic/inc/uqMatrixOfStuff.h
 *
 * Copyright (C) 2008 The PECOS Team, http://queso.ices.utexas.edu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_MATRIX_OF_STUFF_H__
#define __UQ_MATRIX_OF_STUFF_H__

template <class T>
class uqMatrixOfStuff
{
public:
  uqMatrixOfStuff(unsigned int numRows, unsigned int numCols);
 ~uqMatrixOfStuff();

  const unsigned int numRows    ()                                 const;
  const unsigned int numCols    ()                                 const;
        void         setLocation(unsigned int i, unsigned int j, T* info);
        T&           operator   ()(unsigned int i, unsigned int j);
  const T&           operator   ()(unsigned int i, unsigned int j) const;

private:
  unsigned int m_numRows;
  unsigned int m_numCols;

  std::vector<std::vector<T*>* > m_data;
};

template <class T>
uqMatrixOfStuff<T>::uqMatrixOfStuff(
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
uqMatrixOfStuff<T>::~uqMatrixOfStuff()
{
  for (unsigned int i = 0; i < m_numRows; ++i) {
    for (unsigned int j = 0; j < m_numCols; ++j) {
      if ((*(m_data[i]))[j] != NULL) delete (*(m_data[i]))[j];
    }
    delete m_data[i];
  }
}

template <class T>
const unsigned int
uqMatrixOfStuff<T>::numRows() const
{
  return m_numRows;
}

template <class T>
const unsigned int
uqMatrixOfStuff<T>::numCols() const
{
  return m_numCols;
}

template <class T>
void
uqMatrixOfStuff<T>::setLocation(unsigned int i, unsigned int j, T* info)
{
  (*(m_data[i]))[j] = info;
  return;
}

template <class T>
T&
uqMatrixOfStuff<T>::operator()(unsigned int i, unsigned int j)
{
  return *(*(m_data[i]))[j];
}

template <class T>
const T&
uqMatrixOfStuff<T>::operator()(unsigned int i, unsigned int j) const
{
  return *(*(m_data[i]))[j];
}
#endif // __UQ_MATRIX_OF_STUFF_H__

