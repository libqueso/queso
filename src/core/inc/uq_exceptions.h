//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#ifndef QUESO_EXCEPTIONS_H
#define QUESO_EXCEPTIONS_H

// C++
#include <stdexcept>
#include <string>

namespace Queso
{
  /*!
   * A class to represent the internal "this should never happen"
   * errors, to be thrown by "queso_error();"
   */
  class LogicError : public std::logic_error
  {
  public:
    LogicError() : std::logic_error( "Error in Queso internal logic" ) {}
  };

  /*!
   * A class to stub for features that should be in Queso, but
   * haven't been written yet, to be thrown by
   * "queso_not_implemented();"
   */
  class NotImplemented : public std::logic_error
  {
    public:
    NotImplemented() : std::logic_error( "Error: not implemented!" ) {}
  };
  
  /*!
   * A class representing a failed attempt by the library to open a
   * file (or construct an fstream, etc), to be thrown by
   * "queso_file_error(filename);" For ease of debugging, "filename"
   * should include any (absolute or relative or implicit) pathname
   * that was part of the failed open.
   */
  class FileError : public std::runtime_error
  {
  public:
    FileError(const std::string& filename) : std::runtime_error( "Error accessing file: " + filename ) {}
  };

  /*!
   * A class representing failed Newton convergence in T/Tv inversion.
   */
  class FailedNewtonTTvInversion : public std::runtime_error
  {
  public:
    FailedNewtonTTvInversion (const std::string &description) 
      : std::runtime_error (description)
    {}
  };

} // end namespace Queso

#define QUESO_THROW(e) do { throw e; } while (0)

#endif // QUESO_EXCEPTIONS_H
