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

#ifndef __UQ_ARRAY_OF_ONE_D_GRIDS_H__
#define __UQ_ARRAY_OF_ONE_D_GRIDS_H__

#include <queso/OneDGrid.h>

namespace QUESO {

/*!\file uqArrayOfOneDGrids.h
 * \brief Class to accommodate arrays of one-dimensional grid.
 * 
 * \class ArrayOfOneDGrids
 * \brief Class to accommodate arrays of one-dimensional grid.
 * 
 * Arrays of one-dimensional grids are necessary in the calculation, for instance, of CDFs
 * and MDF of vector functions (refer to BaseVectorCdf, BaseVectorMdf, and 
 * derived classes).
 */
template <class V, class M>
class ArrayOfOneDGrids
{
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  ArrayOfOneDGrids(const char* prefix, const VectorSpace<V,M>& rowSpace);
  
  //! Destructor.
  ~ArrayOfOneDGrids();
  //@}
  
  //! @name Property methods
  //@{
  //! Returns the (vector) space to which the row belongs to.
  const VectorSpace<V,M>&     rowSpace       () const;
  
  //! Returns an array with the sizes of the grids.
  const V&  sizes          () const;
  
  //! Returns an array with the minimum position of each grid.
  const V&  minPositions   () const;
  
  //! Returns an array with the maximum position of each grid.
  const V&  maxPositions   () const;
  //@}
  
  //! @name Math methods
  //@{
//  void      setGrid        (unsigned int                 rowId,
//			    BaseOneDGrid<double>& oneDGrid);
  
  //! Sets an array of uniform grids.
  void      setUniformGrids(const V& sizesVec,
			    const V& minPositionsVec,
			    const V& maxPositionsVec);
  //@}

  //! @name Accessor methods
  //@{
  //! Returns the grid stored in the <c>rowId</c>-th position of the array of grids. 
  const BaseOneDGrid<double>& grid           (unsigned int rowId) const;
  //@}
  
  //! @name I/O methods
  //@{
  //! Prints the values of the array of grids (points). 
  void      print          (std::ostream& os) const;
  //@}
  
private:
  const BaseEnvironment&                  m_env;
        std::string                              m_prefix;
  const VectorSpace<V,M>&                 m_rowSpace;
  DistArray<BaseOneDGrid<double>*> m_oneDGrids;

  V* m_sizes;
  V* m_minPositions;
  V* m_maxPositions;

};
// Default constructor -------------------------------------------------
template <class V, class M>
ArrayOfOneDGrids<V,M>::ArrayOfOneDGrids(
  const char*                    prefix,
  const VectorSpace<V,M>& rowSpace)
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
// Destructor ----------------------------------------------------------
template <class V, class M>
ArrayOfOneDGrids<V,M>::~ArrayOfOneDGrids()
{
  if (m_maxPositions) delete m_maxPositions;
  if (m_minPositions) delete m_minPositions;
  if (m_sizes       ) delete m_sizes;

  for (unsigned int i = 0; i < (unsigned int) m_oneDGrids.MyLength(); ++i) {
    if (m_oneDGrids(i,0)) delete m_oneDGrids(i,0);
  }
}
// Property methods-----------------------------------------------------
template <class V, class M>
const VectorSpace<V,M>&
ArrayOfOneDGrids<V,M>::rowSpace() const
{
  return m_rowSpace;
}
//----------------------------------------------------------------------
template <class V, class M>
const V&
ArrayOfOneDGrids<V,M>::sizes() const
{
  UQ_FATAL_TEST_MACRO(m_sizes == NULL,
                      m_env.worldRank(),
                      "ArrayOfOneDGrids<T>::sizes()",
                      "sizes is still NULL");

  return *m_sizes;
}
//----------------------------------------------------------------------
template <class V, class M>
const V&
ArrayOfOneDGrids<V,M>::minPositions() const
{
  UQ_FATAL_TEST_MACRO(m_minPositions == NULL,
                      m_env.worldRank(),
                      "ArrayOfOneDGrids<T>::minPositions()",
                      "minPositions is still NULL");

  return *m_minPositions;
}
//----------------------------------------------------------------------
template <class V, class M>
const V&
ArrayOfOneDGrids<V,M>::maxPositions() const
{
  UQ_FATAL_TEST_MACRO(m_maxPositions == NULL,
                      m_env.worldRank(),
                      "ArrayOfOneDGrids<T>::maxPositions()",
                      "maxPositions is still NULL");

  return *m_maxPositions;
}
//----------------------------------------------------------------------
template <class V, class M>
void
ArrayOfOneDGrids<V,M>::setUniformGrids(
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
    m_oneDGrids(i,0) = new UniformOneDGrid<double>(m_env,
                                                          (m_prefix+strI).c_str(),
                                                          (unsigned int) sizesVec[i],
                                                          minPositionsVec[i],
                                                          maxPositionsVec[i]);
  }

  return;
}
//----------------------------------------------------------------------
template <class V, class M>
const BaseOneDGrid<double>&
ArrayOfOneDGrids<V,M>::grid(unsigned int rowId) const
{
  UQ_FATAL_TEST_MACRO(rowId >= m_rowSpace.dimLocal(),
                      m_env.worldRank(),
                      "ArrayOfOneDUnformGrids<T>::grid()",
                      "rowId is out of range");

  ArrayOfOneDGrids<V,M>* tmp = const_cast<ArrayOfOneDGrids<V,M>*>(this);
  return *(tmp->m_oneDGrids(rowId,0));
}
// I/O methods----------------------------------------------------------
template <class V, class M>
void
ArrayOfOneDGrids<V,M>::print(std::ostream& os) const
{
  ArrayOfOneDGrids<V,M>* tmp = const_cast<ArrayOfOneDGrids<V,M>*>(this);
  for (unsigned int i = 0; i < (unsigned int) m_oneDGrids.MyLength(); ++i) {
    os << *(tmp->m_oneDGrids(i,0))
       << std::endl;
  }

  return;
}
//----------------------------------------------------------------------
template <class V, class M>
std::ostream& operator<< (std::ostream& os, const ArrayOfOneDGrids<V,M>& obj)
{
  obj.print(os);
  return os;
}

}  // End namespace QUESO

#endif // __UQ_ARRAY_OF_ONE_D_GRIDS_H__
