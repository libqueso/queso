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
  virtual void assemble() = 0;

  //! Print libmesh related information
  virtual void print_info() const = 0;

  //! Save the eigenvalues to file \c filename
  virtual void save_converged_evals(const std::string &filename) const;

  //! Save converged eigenfunction \c i to \c filename
  virtual void save_converged_evec(const std::string &filename, unsigned int i) const;

protected:
  libMesh::Mesh *mesh;
  libMesh::EquationSystems *equation_systems;

  //! The number of converged eigenvalue/eigenvector pairs
  unsigned int nconv;

private:
  //! Common initialisation function for all the constructors
  void init();
};

#endif // __QUESO_LIBMESHOPERATOR_BASE__
