//-----------------------------------------------------------------------bl-
//
//-----------------------------------------------------------------------el-

#ifndef __QUESO_OPERATOR_BASE__
#define __QUESO_OPERATOR_BASE__

/*!
 * \file uqOperatorBase.h
 * \brief Abstract base class for operator objects
 */

#include <string>

class uqOperatorBase {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor
  uqOperatorBase();

  //! Construct using the mesh from file \c filename
  uqOperatorBase(const std::string& filename);

  //! Destructor
  ~uqOperatorBase();
  //@}
};

#endif // __QUESO_OPERATOR_BASE__
