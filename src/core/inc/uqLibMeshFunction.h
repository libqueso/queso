//-----------------------------------------------------------------------bl-
//
//-----------------------------------------------------------------------el-

#ifndef __QUESO_LIBMESHFUNCTION__
#define __QUESO_LIBMESHFUNCTION__

#include <uqFunctionBase.h>

namespace libMesh {
  class Mesh;
  class EquationSystems;
}

/*!
 * \file uqLibMeshFunction.h
 * \brief Function objects using lib mesh for the backend
 */

class uqLibMeshFunction : public uqFunctionBase {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor. Zero everywhere.
  uqLibMeshFunction();

  //! Create a function that is equal to \c value everywhere
  uqLibMeshFunction(double value);

  //! Destructor
  ~uqLibMeshFunction();

  //! Will print mesh-related libMesh foo to the screen
  void print_info();

  //@}
private:
  libMesh::Mesh *mesh;
  libMesh::EquationSystems *equations_systems;
};

#endif // __QUESO_LIBMESHFUNCTION__
