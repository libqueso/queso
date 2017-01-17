//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2017 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

//
// nick (nicholas.malaya) providing a direct port of the antioch assertions
// to queso (see: https://github.com/libantioch/antioch/blob/master/src/utilities/include/antioch/antioch_asserts.h)
//
// adding MPI-friendly exits
//

#ifndef QUESO_ASSERTS_H
#define QUESO_ASSERTS_H

#include "exceptions.h"
#include "queso/config_queso.h" // for QUESO_HAVE_CXX11

// C++
#include <iostream>
#include <iomanip>

#define queso_here()     do { std::cerr << __FILE__ << ", line " << __LINE__ << ", compiled " << __DATE__ << " at " << __TIME__ << std::endl; } while (0)

// queso_error and kin throw exceptions to indicate various possible
// errors

#define queso_error_msg(msg) do { queso_here(); std::cerr << msg << std::endl; QUESO_THROW(QUESO::LogicError()); } while(0)

#define queso_not_implemented_msg(msg) do { queso_here(); std::cerr << msg << std::endl; QUESO_THROW(QUESO::NotImplemented()); } while(0)

#define queso_file_error_msg(filename, msg) do { queso_here(); std::cerr << msg << std::endl; QUESO_THROW(QUESO::FileError(filename)); } while(0)

#define queso_error() \
  queso_error_msg("")

#define queso_not_implemented() \
  queso_not_implemented_msg("")

#define queso_file_error(filename) \
  queso_file_error_msg(filename, "")

// The queso_assert() macro acts like C's assert(), but throws a
// queso_error() (enabling exception handling, stack trace, etc)
// instead of just exiting.

// The queso_require() macro does the same, but remains active even
// when NDEBUG is not defined

#define queso_require_msg(asserted, msg)  do { if (!(asserted)) { std::cerr << "Assertion `" #asserted "' failed.\n" << msg << std::endl; queso_error(); } } while(0)

// We can get more descriptive error messages if we're testing a
// simple standard binary operation.
#define queso_require_equal_to_msg(expr1,expr2,msg)  do { if (!(expr1 == expr2)) { std::cerr << "Assertion `" #expr1 " == " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; queso_error(); } } while(0)
#define queso_require_not_equal_to_msg(expr1,expr2,msg)  do { if (!(expr1 != expr2)) { std::cerr << "Assertion `" #expr1 " != " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; queso_error(); } } while(0)
#define queso_require_less_msg(expr1,expr2,msg)  do { if (!(expr1 < expr2)) { std::cerr << "Assertion `" #expr1 " < " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; queso_error(); } } while(0)
#define queso_require_greater_msg(expr1,expr2,msg)  do { if (!(expr1 > expr2)) { std::cerr << "Assertion `" #expr1 " > " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; queso_error(); } } while(0)
#define queso_require_less_equal_msg(expr1,expr2,msg)  do { if (!(expr1 <= expr2)) { std::cerr << "Assertion `" #expr1 " <= " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; queso_error(); } } while(0)
#define queso_require_greater_equal_msg(expr1,expr2,msg)  do { if (!(expr1 >= expr2)) { std::cerr << "Assertion `" #expr1 " >= " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; queso_error(); } } while(0)

// When not debugging, we don't test any asserts
#ifdef NDEBUG

#define queso_assert_msg(asserted,msg)  ((void) 0)
#define queso_assert_equal_to_msg(expr1,expr2,msg)  ((void) 0)
#define queso_assert_not_equal_to_msg(expr1,expr2,msg)  ((void) 0)
#define queso_assert_less_msg(expr1,expr2,msg)  ((void) 0)
#define queso_assert_greater_msg(expr1,expr2,msg)  ((void) 0)
#define queso_assert_less_equal_msg(expr1,expr2,msg)  ((void) 0)
#define queso_assert_greater_equal_msg(expr1,expr2,msg)  ((void) 0)

#else

#define queso_assert_msg(asserted,msg) \
  queso_require_msg(asserted,msg)
#define queso_assert_equal_to_msg(expr1,expr2,msg) \
  queso_require_equal_to_msg(expr1,expr2,msg)
#define queso_assert_not_equal_to_msg(expr1,expr2,msg) \
  queso_require_not_equal_to_msg(expr1,expr2,msg)
#define queso_assert_less_msg(expr1,expr2,msg) \
  queso_require_less_msg(expr1,expr2,msg)
#define queso_assert_greater_msg(expr1,expr2,msg) \
  queso_require_greater_msg(expr1,expr2,msg)
#define queso_assert_less_equal_msg(expr1,expr2,msg) \
  queso_require_less_equal_msg(expr1,expr2,msg)
#define queso_assert_greater_equal_msg(expr1,expr2,msg) \
  queso_require_greater_equal_msg(expr1,expr2,msg)

#endif // NDEBUG

#define queso_require(asserted) \
  queso_require_msg(asserted, "")

#define queso_require_equal_to(expr1,expr2) \
  queso_require_equal_to_msg(expr1,expr2,"")

#define queso_require_not_equal_to(expr1,expr2) \
  queso_require_not_equal_to_msg(expr1,expr2,"") \

#define queso_require_less(expr1,expr2) \
  queso_require_less_msg(expr1,expr2,"")

#define queso_require_greater(expr1,expr2) \
  queso_require_greater_msg(expr1,expr2,"")

#define queso_require_less_equal(expr1,expr2) \
  queso_require_less_equal_msg(expr1,expr2,"")

#define queso_require_greater_equal(expr1,expr2) \
  queso_require_greater_equal_msg(expr1,expr2,"")


#define queso_assert(asserted) \
  queso_assert_msg(asserted, "")

#define queso_assert_equal_to(expr1,expr2) \
  queso_assert_equal_to_msg(expr1,expr2,"")

#define queso_assert_not_equal_to(expr1,expr2) \
  queso_assert_not_equal_to_msg(expr1,expr2,"") \

#define queso_assert_less(expr1,expr2) \
  queso_assert_less_msg(expr1,expr2,"")

#define queso_assert_greater(expr1,expr2) \
  queso_assert_greater_msg(expr1,expr2,"")

#define queso_assert_less_equal(expr1,expr2) \
  queso_assert_less_equal_msg(expr1,expr2,"")

#define queso_assert_greater_equal(expr1,expr2) \
  queso_assert_greater_equal_msg(expr1,expr2,"")



#endif // QUESO_ASSERTS_H
