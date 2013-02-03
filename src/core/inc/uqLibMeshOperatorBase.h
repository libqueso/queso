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
  
  // virtual void assemble_matrices(libMesh::EquationSystems &es,
  //     const std::string& name) = 0;

  // virtual void get_boundary_dofs(libMesh::EquationSystems& es,
  //     const std::string& system_name,
  //     std::set<unsigned int>& global_boundary_dofs_set) = 0;
protected:
  libMesh::Mesh *mesh;
  libMesh::EquationSystems *equation_systems;
  // libMesh::EigenSystem *eigen_system;
};

#endif // __QUESO_LIBMESHOPERATOR_BASE__
