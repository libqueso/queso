//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010 The PECOS Development Team
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

#ifndef __UQ_ARRAY_OF_ONE_D_GRIDS_H__
#define __UQ_ARRAY_OF_ONE_D_GRIDS_H__

#include <uqOneDGrid.h>

template <class V, class M>
class uqArrayOfOneDGridsClass
{
public:
  uqArrayOfOneDGridsClass(const char* prefix, const uqVectorSpaceClass<V,M>& rowSpace);
 ~uqArrayOfOneDGridsClass();

  const uqVectorSpaceClass<V,M>&     rowSpace       () const;
  const V&                           sizes          () const;
  const V&                           minPositions   () const;
  const V&                           maxPositions   () const;
//      void                         setGrid        (unsigned int                 rowId,
//                                                   uqBaseOneDGridClass<double>& oneDGrid);
        void                         setUniformGrids(const V& sizesVec,
                                                     const V& minPositionsVec,
                                                     const V& maxPositionsVec);
  const uqBaseOneDGridClass<double>& grid           (unsigned int rowId) const;
        void                         print          (std::ostream& os) const;

private:
  const uqBaseEnvironmentClass&                  m_env;
        std::string                              m_prefix;
  const uqVectorSpaceClass<V,M>&                 m_rowSpace;
  uqDistArrayClass<uqBaseOneDGridClass<double>*> m_oneDGrids;

  V* m_sizes;
  V* m_minPositions;
  V* m_maxPositions;

};

template <class V, class M>
uqArrayOfOneDGridsClass<V,M>::uqArrayOfOneDGridsClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& rowSpace)
  :
  m_env         (rowSpace.env()    ),
  m_prefix      ((std::string)(prefix)+""),
  m_rowSpace    (rowSpace          ),
  m_oneDGrids   (m_rowSpace.map(),1),
  m_sizes       (NULL),
  m_minPositions(NULL),
  m_maxPositions(NULL)
{
  for (unsigned int i = 0; i < (unsigned int) m_oneDGrids.MyLength(); ++i) {
    m_oneDGrids(i,0) = NULL;
  }
}

template <class V, class M>
uqArrayOfOneDGridsClass<V,M>::~uqArrayOfOneDGridsClass()
{
  if (m_maxPositions) delete m_maxPositions;
  if (m_minPositions) delete m_minPositions;
  if (m_sizes       ) delete m_sizes;

  for (unsigned int i = 0; i < (unsigned int) m_oneDGrids.MyLength(); ++i) {
    if (m_oneDGrids(i,0)) delete m_oneDGrids(i,0);
  }
}

template <class V, class M>
const uqVectorSpaceClass<V,M>&
uqArrayOfOneDGridsClass<V,M>::rowSpace() const
{
  return m_rowSpace;
}

template <class V, class M>
const V&
uqArrayOfOneDGridsClass<V,M>::sizes() const
{
  UQ_FATAL_TEST_MACRO(m_sizes == NULL,
                      m_env.worldRank(),
                      "uqArrayOfOneDGridsClass<T>::sizes()",
                      "sizes is still NULL");

  return *m_sizes;
}

template <class V, class M>
const V&
uqArrayOfOneDGridsClass<V,M>::minPositions() const
{
  UQ_FATAL_TEST_MACRO(m_minPositions == NULL,
                      m_env.worldRank(),
                      "uqArrayOfOneDGridsClass<T>::minPositions()",
                      "minPositions is still NULL");

  return *m_minPositions;
}

template <class V, class M>
const V&
uqArrayOfOneDGridsClass<V,M>::maxPositions() const
{
  UQ_FATAL_TEST_MACRO(m_maxPositions == NULL,
                      m_env.worldRank(),
                      "uqArrayOfOneDGridsClass<T>::maxPositions()",
                      "maxPositions is still NULL");

  return *m_maxPositions;
}

template <class V, class M>
void
uqArrayOfOneDGridsClass<V,M>::setUniformGrids(
  const V& sizesVec,
  const V& minPositionsVec,
  const V& maxPositionsVec)
{
  if (m_sizes        == NULL) m_sizes        = new V(sizesVec);
  else                       *m_sizes        = sizesVec;

  if (m_minPositions == NULL) m_minPositions = new V(minPositionsVec);
  else                       *m_minPositions = minPositionsVec;

  if (m_maxPositions == NULL) m_maxPositions = new V(maxPositionsVec);
  else                       *m_maxPositions = maxPositionsVec;

  char strI[65];
  for (unsigned int i = 0; i < (unsigned int) m_oneDGrids.MyLength(); ++i) {
    sprintf(strI,"%u_",i);
    m_oneDGrids(i,0) = new uqUniformOneDGridClass<double>(m_env,
                                                          (m_prefix+strI).c_str(),
                                                          (unsigned int) sizesVec[i],
                                                          minPositionsVec[i],
                                                          maxPositionsVec[i]);
  }

  return;
}

template <class V, class M>
const uqBaseOneDGridClass<double>&
uqArrayOfOneDGridsClass<V,M>::grid(unsigned int rowId) const
{
  UQ_FATAL_TEST_MACRO(rowId >= m_rowSpace.dimLocal(),
                      m_env.worldRank(),
                      "uqArrayOfOneDUnformGridsClass<T>::grid()",
                      "rowId is out of range");

  uqArrayOfOneDGridsClass<V,M>* tmp = const_cast<uqArrayOfOneDGridsClass<V,M>*>(this);
  return *(tmp->m_oneDGrids(rowId,0));
}

template <class V, class M>
void
uqArrayOfOneDGridsClass<V,M>::print(std::ostream& os) const
{
  uqArrayOfOneDGridsClass<V,M>* tmp = const_cast<uqArrayOfOneDGridsClass<V,M>*>(this);
  for (unsigned int i = 0; i < (unsigned int) m_oneDGrids.MyLength(); ++i) {
    os << *(tmp->m_oneDGrids(i,0))
       << std::endl;
  }

  return;
}

template <class V, class M>
std::ostream& operator<< (std::ostream& os, const uqArrayOfOneDGridsClass<V,M>& obj)
{
  obj.print(os);
  return os;
}
#endif // __UQ_ARRAY_OF_ONE_D_GRIDS_H__
