//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

//
// nick (nicholas.malaya) providing a direct port of the antioch assertions
// to queso (see: https://github.com/libantioch/antioch/blob/master/src/utilities/include/antioch/antioch_asserts.h)
//
// adding MPI-friendly exits
//
#ifndef QUESO_ASSERTS_H
#define QUESO_ASSERTS_H

#include "exceptions.h"
#include "config_queso.h" // for QUESO_HAVE_CXX11

// C++
#include <iostream>
#include <iomanip>

#define queso_here()     do { std::cerr << __FILE__ << ", line " << __LINE__ << ", compiled " << __DATE__ << " at " << __TIME__ << std::endl; } while (0)

// The queso_assert() macro acts like C's assert(), but throws a
// queso_error() (including stack trace, etc) instead of just exiting

// When not debugging, we don't test any asserts
#ifdef NDEBUG
#define queso_assert(asserted)  ((void) 0)
#define queso_assert_msg(asserted, msg)  ((void) 0)
#define queso_assert_equal_to(expr1,expr2)  ((void) 0)
#define queso_assert_not_equal_to(expr1,expr2)  ((void) 0)
#define queso_assert_less(expr1,expr2)  ((void) 0)
#define queso_assert_greater(expr1,expr2)  ((void) 0)
#define queso_assert_less_equal(expr1,expr2)  ((void) 0)
#define queso_assert_greater_equal(expr1,expr2)  ((void) 0)
#else

#define queso_assert(asserted)  do { if (!(asserted)) { std::cerr << "Assertion `" #asserted "' failed." << std::endl; queso_error(); } } while(0)

// When using C++11, we can test asserts comparing two different types
// robustly
// #if __cplusplus > 199711L // http://gcc.gnu.org/bugzilla/show_bug.cgi?id=1773
#ifdef QUESO_HAVE_CXX11
#define queso_assert_equal_to(expr1,expr2)  do { typedef decltype(expr1) type1; typedef decltype(expr2) type2; if (!((expr1 == static_cast<type1>(expr2)) && static_cast<type2>(expr1) == expr2)) { std::cerr << "Assertion `" #expr1 " == " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; queso_error(); } } while(0)
#define queso_assert_not_equal_to(expr1,expr2)  do { typedef decltype(expr1) type1; typedef decltype(expr2) type2; if (!((expr1 != static_cast<type1>(expr2)) && (static_cast<type2>(expr1) != expr2))) { std::cerr << "Assertion `" #expr1 " != " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; queso_error(); } } while(0)
#define queso_assert_less(expr1,expr2)  do { typedef decltype(expr1) type1; typedef decltype(expr2) type2; if (!((static_cast<type2>(expr1) < expr2) && (expr1 < static_cast<type1>(expr2)))) { std::cerr << "Assertion `" #expr1 " < " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; queso_error(); } } while(0)
#define queso_assert_greater(expr1,expr2)  do { typedef decltype(expr1) type1; typedef decltype(expr2) type2; if (!((static_cast<type2>(expr1) > expr2) && (expr1 > static_cast<type1>(expr2)))) { std::cerr << "Assertion `" #expr1 " > " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; queso_error(); } } while(0)
#define queso_assert_less_equal(expr1,expr2)  do { typedef decltype(expr1) type1; typedef decltype(expr2) type2; if (!((static_cast<type2>(expr1) <= expr2) && (expr1 <= static_cast<type1>(expr2)))) { std::cerr << "Assertion `" #expr1 " <= " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; queso_error(); } } while(0)
#define queso_assert_greater_equal(expr1,expr2)  do { typedef decltype(expr1) type1; typedef decltype(expr2) type2; if (!((static_cast<type2>(expr1) >= expr2) && (expr1 >= static_cast<type1>(expr2)))) { std::cerr << "Assertion `" #expr1 " >= " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; queso_error(); } } while(0)

// When using C++98, we let the compiler pick the type conversion and
// hope for the best.
#else
#define queso_assert_equal_to(expr1,expr2)  do { if (!(expr1 == expr2)) { std::cerr << "Assertion `" #expr1 " == " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; queso_error(); } } while(0)
#define queso_assert_not_equal_to(expr1,expr2)  do { if (!(expr1 != expr2)) { std::cerr << "Assertion `" #expr1 " != " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; queso_error(); } } while(0)
#define queso_assert_less(expr1,expr2)  do { if (!(expr1 < expr2)) { std::cerr << "Assertion `" #expr1 " < " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; queso_error(); } } while(0)
#define queso_assert_greater(expr1,expr2)  do { if (!(expr1 > expr2)) { std::cerr << "Assertion `" #expr1 " > " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; queso_error(); } } while(0)
#define queso_assert_less_equal(expr1,expr2)  do { if (!(expr1 <= expr2)) { std::cerr << "Assertion `" #expr1 " <= " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; queso_error(); } } while(0)
#define queso_assert_greater_equal(expr1,expr2)  do { if (!(expr1 >= expr2)) { std::cerr << "Assertion `" #expr1 " >= " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; queso_error(); } } while(0)

#endif // C++11

#endif // NDEBUG


#define queso_error()    do { queso_here(); QUESO_THROW(QUESO::LogicError()); } while(0)
#define queso_not_implemented()    do { queso_here(); QUESO_THROW(QUESO::NotImplemented()); } while(0)
#define queso_file_error(filename)    do { queso_here(); QUESO_THROW(QUESO::FileError(filename)); } while(0)

#endif // QUESO_ASSERTS_H
