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

#ifndef UQ_SCOPED_PTR_H
#define UQ_SCOPED_PTR_H

#include <queso/config_queso.h>

#ifdef QUESO_HAVE_CXX11_UNIQUE_PTR
#include <memory>
#elif QUESO_HAVE_BOOST_SCOPED_PTR_HPP
#include <boost/scoped_ptr.hpp>
#elif QUESO_HAVE_CXX_AUTO_PTR
#include <memory>
#endif

/*!
 * \file ScopedPtr.h
 * \brief This file declares and defines a scoped pointer
 */

/*!
 * \class ScopedPtr
 * \brief Definition of a scoped pointer.
 *
 * Scoped pointers are smart pointers that clean up after themselves when they
 * go out of scope, with no shared ownership, and with no transfer of
 * ownership regardless of the underlying implementation.
 *
 * If QUESO detects C++11 functionality at configure-time, then the underlying
 * implementation of ScopedPtr is that of std::unique_ptr.  Even when this is
 * the case, please note that shared ownership and transfer of ownership are
 * not supported.
 *
 * If QUESO does not detect C++11 functionality then QUESO will check if boost
 * is present and, if so, will use boost::scoped_ptr as the underlying
 * implementation.  It should now be clear why shared ownership and transfer
 * of ownership are not permitted.
 *
 * If QUESO fails to detect both C++11 functionality and boost then QUESO's
 * last resort is to fall back on the C++03's std::auto_ptr.  If std::auto_ptr
 * isn't detected, QUESO will fail to configure.
 *
 * ScopedPtr is implemented as a template typedef pattern.  That is, one
 * declares a ScopedPtr like so:
 *
 * QUESO::ScopedPtr< name_of_type >::Type name_of_variable;
 */
namespace QUESO
{
#ifdef QUESO_HAVE_CXX11_UNIQUE_PTR
  template<typename T>
  struct ScopedPtr
  {
    typedef std::unique_ptr<T> Type;
  };
#elif QUESO_HAVE_BOOST_SCOPED_PTR_HPP
  template<typename T>
  struct ScopedPtr
  {
    typedef boost::scoped_ptr<T> Type;
  };
#elif QUESO_HAVE_CXX_AUTO_PTR
  template<typename T>
  struct ScopedPtr
  {
    typedef std::auto_ptr<T> Type;
  };
#else
#     error "No valid definition for ScopedPtr found!"
#endif

} // end namespace QUESO

#endif // UQ_SCOPED_PTR_H
