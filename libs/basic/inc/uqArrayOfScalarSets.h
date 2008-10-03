/* uq/libs/basic/inc/uqArrayOfScalarSets.h
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

#ifndef __UQ_ARRAY_OF_SCALAR_SETS_H__
#define __UQ_ARRAY_OF_SCALAR_SETS_H__

#include <EpetraExt_DistArray.h>

template <class V, class M>
class uqArrayOfScalarSetsClass
{
public:
  uqArrayOfScalarSetsClass(const uqVectorSpaceClass<V,M>& rowSpace);
 ~uqArrayOfScalarSetsClass();

  void                       setScalarSet(unsigned int rowId, const std::vector<double>& values);
  const std::vector<double>& scalarSet   (unsigned int rowId) const;

private:
  const uqEnvironmentClass&                  m_env;
  const uqVectorSpaceClass<V,M>&             m_rowSpace;
  EpetraExt::DistArray<std::vector<double>*> m_scalarSets;
};

template <class V, class M>
uqArrayOfScalarSetsClass<V,M>::uqArrayOfScalarSetsClass(const uqVectorSpaceClass<V,M>& rowSpace)
  :
  m_env       (rowSpace.env()),
  m_rowSpace  (rowSpace      ),
  m_scalarSets(m_rowSpace.map(),1)
{
  for (unsigned int i = 0; i < (unsigned int) m_scalarSets.MyLength(); ++i) {
    m_scalarSets(i,0) = NULL;
  }
}

template <class V, class M>
uqArrayOfScalarSetsClass<V,M>::~uqArrayOfScalarSetsClass()
{
  for (unsigned int i = 0; i < (unsigned int) m_scalarSets.MyLength(); ++i) {
    if (m_scalarSets(i,0)) delete m_scalarSets(i,0);
  }
}

template <class V, class M>
void
uqArrayOfScalarSetsClass<V,M>::setScalarSet(unsigned int rowId, const std::vector<double>& values)
{
  UQ_FATAL_TEST_MACRO(rowId >= (unsigned int) m_scalarSets.MyLength(),
                      m_env.rank(),
                      "uqArrayOfScalarSetClass<T>::setScalarSet()",
                      "rowId is out of range");

  if (m_scalarSets(rowId,0) == NULL) m_scalarSets(rowId,0) = new std::vector<double>(0);
  else                               m_scalarSets(rowId,0)->clear();

  std::vector<double>& vec = *(m_scalarSets(rowId,0));
  vec.resize(values.size(),0.);
  for (unsigned int j = 0; j < values.size(); ++j) {
    vec[j] = values[j];
  }

  return;
}

template <class V, class M>
const std::vector<double>&
uqArrayOfScalarSetsClass<V,M>::scalarSet(unsigned int rowId) const
{
  UQ_FATAL_TEST_MACRO(rowId >= (unsigned int) m_scalarSets.MyLength(),
                      m_env.rank(),
                      "uqArrayOfScalarSetClass<T>::scalarSet()",
                      "rowId is out of range");

  uqArrayOfScalarSetsClass<V,M>* tmp = const_cast<uqArrayOfScalarSetsClass<V,M>*>(this);

  UQ_FATAL_TEST_MACRO(tmp->m_scalarSets(rowId,0) == NULL,
                      m_env.rank(),
                      "uqArrayOfScalarSetClass<T>::scalarSet()",
                      "requested row is still NULL");

  return *(tmp->m_scalarSets(rowId,0));
}

#endif // __UQ_ARRAY_OF_SCALAR_SETS_H__

