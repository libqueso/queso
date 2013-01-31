//-----------------------------------------------------------------------bl-
//
//-----------------------------------------------------------------------el-

#ifndef __QUESO_FUNCTION_BASE__
#define __QUESO_FUNCTION_BASE__

/*!
 * \file uqFunctionBase.h
 * \brief Abstract base class for function objects
 */

class uqFunctionBase {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor. Zero everywhere.
  uqFunctionBase();

  //   libmesh object
  //   read from a file

  //! Destructor
  ~uqFunctionBase();
  //@}

protected:
  // Number of degrees of freedom
  unsigned int ndofs;
};

#endif // __QUESO_FUNCTION_BASE__
