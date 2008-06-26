#include <uqFinDimLinearSpace.h>
#include <uqGslMatrix.h>

template<>
uqGslVectorClass*
uqFinDimLinearSpaceClass<uqGslVectorClass,uqGslMatrixClass>::newVector() const
{
  return new uqGslVectorClass(m_env,this->dim());
}

template<>
uqGslMatrixClass*
uqFinDimLinearSpaceClass<uqGslVectorClass,uqGslMatrixClass>::newMatrix() const
{
  return new uqGslMatrixClass(m_env,this->dim(),this->dim());
}

template<>
uqGslMatrixClass*
uqFinDimLinearSpaceClass<uqGslVectorClass,uqGslMatrixClass>::newDiagMatrix(double diagValue) const
{
  return new uqGslMatrixClass(m_env,this->dim(),diagValue);
}
