//-----------------------------------------------------------------------bl-
//
//-----------------------------------------------------------------------el-

#ifndef __QUESO_LIBMESHOPERATOR_BASE__
#define __QUESO_LIBMESHOPERATOR_BASE__

/*!
 * \file uqLibMeshOperatorBase.h
 * \brief Abstract base class for operator objects using libmesh in the
 *        backend
 */

#include <string>
#include <set>
#include <uqOperatorBase.h>
#include <libmesh/system.h>

namespace libMesh {
  class Mesh;
  class EquationSystems;
  class EigenSystem;
}

class uqLibMeshOperatorBase : public uqOperatorBase, public libMesh::System::Assembly {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor
  uqLibMeshOperatorBase();

  //! Construct an operator on the meshed domain in \c filename
  uqLibMeshOperatorBase(const std::string& filename);

  //! Destructor
  ~uqLibMeshOperatorBase();
  //@}

  //! Must implement this for the solve to work
  virtual void assemble()=0;
  
protected:
  libMesh::Mesh *mesh;
  libMesh::EquationSystems *equation_systems;
};

#endif // __QUESO_LIBMESHOPERATOR_BASE__
