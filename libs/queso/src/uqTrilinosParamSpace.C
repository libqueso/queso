/* uq/libs/queso/src/uqTrilinosParamSpace.C
 *
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <uqTrilinosParamSpace.h>

#if 0
uqTrilinosParamSpaceClass::uqTrilinosParamSpaceClass(const Epetra_MpiComm& comm, unsigned int dimension)
  :
  uqParamSpaceClass<uqTrilinosVectorClass,uqTrilinosMatrixClass>(dimension),
  //m_rank(comm.MyPID()), // compiler complains FIXME
  m_map(new Epetra_Map(dimension,0,comm))
{
  //std::cout << "Entering uqTrilinosParamSpaceClass::constructor()"
  //          << std::endl;

  m_rank = comm.MyPID();

  //std::cout << "Leaving uqTrilinosParamSpaceClass::constructor()"
  //          << std::endl;
}
  
uqTrilinosParamSpaceClass::~uqTrilinosParamSpaceClass()
{
  std::cout << "Entering uqTrilinosParamSpaceClass::destructor()"
            << std::endl;

  if (m_map) delete m_map;

  std::cout << "Leaving uqTrilinosParamSpaceClass::destructor()"
            << std::endl;
}

uqTrilinosVectorClass*
uqTrilinosParamSpaceClass::newVector() const
{
  return new uqTrilinosVectorClass(*m_map);
}

uqTrilinosMatrixClass*
uqTrilinosParamSpaceClass::newMatrix() const
{
  return new uqTrilinosMatrixClass(*m_map,this->dim());
}
#endif

template<>
void
uqParamSpaceClass<uqTrilinosVectorClass,uqTrilinosMatrixClass>::createInitialValues() const
{
  return;
}

template<>
void
uqParamSpaceClass<uqTrilinosVectorClass,uqTrilinosMatrixClass>::createMinValues() const
{
  return;
}

template<>
void
uqParamSpaceClass<uqTrilinosVectorClass,uqTrilinosMatrixClass>::createMaxValues() const
{
  return;
}

template<>
void
uqParamSpaceClass<uqTrilinosVectorClass,uqTrilinosMatrixClass>::createPriorMuValues() const
{
  return;
}

template<>
void
uqParamSpaceClass<uqTrilinosVectorClass,uqTrilinosMatrixClass>::createPriorSigmaValues() const
{
  return;
}

template<>
void
uqParamSpaceClass<uqTrilinosVectorClass,uqTrilinosMatrixClass>::print(std::ostream& os) const
{
  return;
}

#if 0
template<>
const Epetra_Map&
uqParamSpaceClass<uqTrilinosVectorClass,uqTrilinosMatrixClass>::map() const
{
 return *m_map;
}

std::ostream&
operator<<(std::ostream& os, const uqParamSpaceClass& space)
{
  space.uqParamSpaceClass<uqTrilinosVectorClass,uqTrilinosMatrixClass>::print(os);
  space.print(os);

  return os;
}
#endif
