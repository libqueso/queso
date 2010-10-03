//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
// 
// $Id$
//
//--------------------------------------------------------------------------

#include <uqVectorSpace.h>
#include <uqTrilinosMatrix.h>

template <>
Epetra_Map*
uqVectorSpaceClass<uqTrilinosVectorClass, uqTrilinosMatrixClass>::newMap()
{
  return new Epetra_Map(m_dimGlobal,0,m_env.subComm());
}

template<>
uqTrilinosVectorClass*
uqVectorSpaceClass<uqTrilinosVectorClass,uqTrilinosMatrixClass>::newVector() const
{
  return new uqTrilinosVectorClass(m_env,*m_map);
}

template<>
uqTrilinosVectorClass*
uqVectorSpaceClass<uqTrilinosVectorClass,uqTrilinosMatrixClass>::newVector(double value) const
{
  return new uqTrilinosVectorClass(m_env,*m_map,value);
}

template<>
uqTrilinosMatrixClass*
uqVectorSpaceClass<uqTrilinosVectorClass,uqTrilinosMatrixClass>::newMatrix() const
{
  return new uqTrilinosMatrixClass(m_env,*m_map,this->dimGlobal());
}

template<>
uqTrilinosMatrixClass*
uqVectorSpaceClass<uqTrilinosVectorClass,uqTrilinosMatrixClass>::newDiagMatrix(double diagValue) const
{
  return new uqTrilinosMatrixClass(m_env,*m_map,diagValue);
}
