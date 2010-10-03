//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
// 
// $Id:$
//
//--------------------------------------------------------------------------

#include <uqVectorSpace.h>
#include <uqGslMatrix.h>

template <>
Epetra_Map*
uqVectorSpaceClass<uqGslVectorClass, uqGslMatrixClass>::newMap()
{
  return new Epetra_Map(m_dimGlobal,0,m_env.selfComm());
}

template<>
uqGslVectorClass*
uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>::newVector() const
{
  return new uqGslVectorClass(m_env,*m_map);
}

template<>
uqGslVectorClass*
uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>::newVector(double value) const
{
  return new uqGslVectorClass(m_env,*m_map,value);
}

template<>
uqGslMatrixClass*
uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>::newMatrix() const
{
  return new uqGslMatrixClass(m_env,*m_map,this->dimGlobal());
}

template<>
uqGslMatrixClass*
uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>::newDiagMatrix(double diagValue) const
{
  return new uqGslMatrixClass(m_env,*m_map,diagValue);
}
