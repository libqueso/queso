//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
// 
// $Id$
//
//--------------------------------------------------------------------------

#include <uqAsciiTable.h>
#ifdef QUESO_HAS_TRILINOS

template <>
Epetra_Map*
uqAsciiTableClass<class uqTrilinosVectorClass, class uqTrilinosMatrixClass>::newMap()
{
  return new Epetra_Map(m_numRows,0,m_env.subComm());
}

#endif // #ifdef QUESO_HAS_TRILINOS
