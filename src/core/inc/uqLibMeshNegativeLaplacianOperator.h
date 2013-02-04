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
  virtual void get_boundary_dofs(
      std::set<unsigned int>& global_boundary_dofs_set);

  // from system::assembly
  virtual void assemble();
};

#endif // __QUESO_LIBMESHNEGATIVELAPLACIANOPERATOR_H__
