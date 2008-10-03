/* uq/libs/basic/inc/uqArrayOfOneDUniformGrids.h
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

#ifndef __UQ_ARRAY_OF_ONE_D_UNIFORM_GRIDS_H__
#define __UQ_ARRAY_OF_ONE_D_UNIFORM_GRIDS_H__

#include <EpetraExt_DistArray.h>

template <class V, class M>
class uqArrayOfOneDUniformGridsClass
{
public:
  uqArrayOfOneDUniformGridsClass(const uqVectorSpaceClass<V,M>& rowSpace);
 ~uqArrayOfOneDUniformGridsClass();

  const uqVectorSpaceClass<V,M>& rowSpace           () const;
  const V&                       numEvaluationPoints() const;
  const V&                       minValues          () const;
  const V&                       maxValues          () const;
        void                     setGrids           (const V& numEvaluationPointsVec,
                                                     const V& minValuesVec,
                                                     const V& maxValuesVec);
        void                     getAGrid           (unsigned int         rowId,
                                                     std::vector<double>& oneDGrid) const;

private:
  const uqEnvironmentClass&      m_env;
  const uqVectorSpaceClass<V,M>& m_rowSpace;
  V*                             m_numEvaluationPoints;
  V*                             m_minValues;
  V*                             m_maxValues;
};

template <class V, class M>
uqArrayOfOneDUniformGridsClass<V,M>::uqArrayOfOneDUniformGridsClass(const uqVectorSpaceClass<V,M>& rowSpace)
  :
  m_env                (rowSpace.env()),
  m_rowSpace           (rowSpace      ),
  m_numEvaluationPoints(NULL),
  m_minValues          (NULL),
  m_maxValues          (NULL)
{
}

template <class V, class M>
uqArrayOfOneDUniformGridsClass<V,M>::~uqArrayOfOneDUniformGridsClass()
{
  if (m_maxValues          ) delete m_maxValues;
  if (m_minValues          ) delete m_minValues;
  if (m_numEvaluationPoints) delete m_numEvaluationPoints;
}

template <class V, class M>
const uqVectorSpaceClass<V,M>&
uqArrayOfOneDUniformGridsClass<V,M>::rowSpace() const
{
  return m_rowSpace;
}

template <class V, class M>
const V&
uqArrayOfOneDUniformGridsClass<V,M>::numEvaluationPoints() const
{
  UQ_FATAL_TEST_MACRO(m_numEvaluationPoints == NULL,
                      m_env.rank(),
                      "uqArrayOfScalarSetClass<T>::numEvaluationPoints()",
                      "numEvaluationPoints is still NULL");

  return *m_numEvaluationPoints;
}

template <class V, class M>
const V&
uqArrayOfOneDUniformGridsClass<V,M>::minValues() const
{
  UQ_FATAL_TEST_MACRO(m_minValues == NULL,
                      m_env.rank(),
                      "uqArrayOfScalarSetClass<T>::minValues()",
                      "minValues is still NULL");

  return *m_minValues;
}

template <class V, class M>
const V&
uqArrayOfOneDUniformGridsClass<V,M>::maxValues() const
{
  UQ_FATAL_TEST_MACRO(m_maxValues == NULL,
                      m_env.rank(),
                      "uqArrayOfScalarSetClass<T>::maxValues()",
                      "maxValues is still NULL");

  return *m_maxValues;
}

template <class V, class M>
void
uqArrayOfOneDUniformGridsClass<V,M>::setGrids(
  const V& numEvaluationPointsVec,
  const V& minValuesVec,
  const V& maxValuesVec)
{
  if (m_numEvaluationPoints == NULL) m_numEvaluationPoints = new V(numEvaluationPointsVec);
  else                              *m_numEvaluationPoints = numEvaluationPointsVec;

  if (m_minValues           == NULL) m_minValues           = new V(minValuesVec);
  else                              *m_minValues           = minValuesVec;

  if (m_maxValues           == NULL) m_maxValues           = new V(maxValuesVec);
  else                              *m_maxValues           = maxValuesVec;

  return;
}

template <class V, class M>
void
uqArrayOfOneDUniformGridsClass<V,M>::getAGrid(
  unsigned int         rowId,
  std::vector<double>& oneDGrid) const
{
  UQ_FATAL_TEST_MACRO(rowId >= m_rowSpace.dim(),
                      m_env.rank(),
                      "uqArrayOfOneDUnformGridsClass<T>::getAGrid()",
                      "rowId is out of range");

  double numEvaluationPoints = (*m_numEvaluationPoints)[rowId];
  double minValue            = (*m_minValues          )[rowId];
  double maxValue            = (*m_maxValues          )[rowId];

  oneDGrid.clear();
  oneDGrid.resize((unsigned int)numEvaluationPoints);
  for (double j = 0.; j < numEvaluationPoints; ++j) {
    double ratio = j/(numEvaluationPoints-1.); // IMPORTANT: Yes, '-1.'
    oneDGrid[(unsigned int)j] = (1.-ratio)*minValue + ratio*maxValue;
  }

  return;
}

#endif // __UQ_ARRAY_OF_ONE_D_UNIFORM_GRIDS_H__
