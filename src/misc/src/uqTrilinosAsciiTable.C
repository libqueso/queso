//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
// 
// $Id$
//
//--------------------------------------------------------------------------

#include <uqAsciiTable.h>

template <>
Epetra_Map*
uqAsciiTableClass<class uqTrilinosVectorClass, class uqTrilinosMatrixClass>::newMap()
{
  return new Epetra_Map(m_numRows,0,m_env.subComm());
}
