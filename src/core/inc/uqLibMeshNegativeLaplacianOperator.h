//-----------------------------------------------------------------------bl-
//
//-----------------------------------------------------------------------el-

#ifndef __QUESO_LIBMESHNEGATIVELAPLACIANOPERATOR_H__
#define __QUESO_LIBMESHNEGATIVELAPLACIANOPERATOR_H__

#include <string>
#include <uqLibMeshOperatorBase.h>

namespace libMesh {
  class EquationSystems;
}

class uqLibMeshNegativeLaplacianOperator : public uqLibMeshOperatorBase {
public:
  // TODO: Document these
  uqLibMeshNegativeLaplacianOperator();
  uqLibMeshNegativeLaplacianOperator(const std::string& filename);
  ~uqLibMeshNegativeLaplacianOperator();

  // from system::assembly
  virtual void assemble();

  //! Print libmesh related information
  virtual void print_info() const;

};

#endif // __QUESO_LIBMESHNEGATIVELAPLACIANOPERATOR_H__
