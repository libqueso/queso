/*For the uqGslVectorSpace.C, it is linked into que inside src/Makefile.am:  
 * libqueso_la_SOURCES += $(top_srcdir)/src/basic/src/uqGslVectorSpace.C
 * For now, I will add it in the makefile
*/
#include <uqVectorSpace.h>
#include <uqTeuchosMatrix.h>

template <>
uqMapClass*
uqVectorSpaceClass<uqTeuchosVectorClass, uqTeuchosMatrixClass>::newMap()
{
  return new uqMapClass(m_dimGlobal,0,m_env.selfComm());
}

template<>
uqTeuchosVectorClass*
uqVectorSpaceClass<uqTeuchosVectorClass,uqTeuchosMatrixClass>::newVector() const
{
  return new uqTeuchosVectorClass(m_env,*m_map);
}

template<>
uqTeuchosVectorClass*
uqVectorSpaceClass<uqTeuchosVectorClass,uqTeuchosMatrixClass>::newVector(double value) const
{
  return new uqTeuchosVectorClass(m_env,*m_map,value);
}

template<>
uqTeuchosMatrixClass*
uqVectorSpaceClass<uqTeuchosVectorClass,uqTeuchosMatrixClass>::newMatrix() const
{
  return new uqTeuchosMatrixClass(m_env,*m_map,this->dimGlobal());
}

template<>
uqTeuchosMatrixClass*
uqVectorSpaceClass<uqTeuchosVectorClass,uqTeuchosMatrixClass>::newDiagMatrix(double diagValue) const
{
  return new uqTeuchosMatrixClass(m_env,*m_map,diagValue);
}
