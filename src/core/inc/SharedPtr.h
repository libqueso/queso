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

#ifndef UQ_SHARED_PTR_H
#define UQ_SHARED_PTR_H

#include <queso/config_queso.h>

#ifdef QUESO_HAVE_CXX11_SHARED_PTR
#include <memory>
#elif QUESO_HAVE_BOOST_SHARED_PTR_HPP
#include <boost/shared_ptr.hpp>
#endif

/*!
 * \file SharedPtr.h
 * \brief This file declares and defines a shared pointer
 */

/*!
 * \class SharedPtr
 * \brief Definition of a shared pointer.
 *
 * Shared pointers may share ownership of some dynamically allocated object.
 * The object pointed to is guaranteed to be deleted when the last underlying
 * implementation of SharedPtr is destroyed or reset.
 *
 * If QUESO detects C++11 functionality ot configure-time, then the underlying
 * implementation of SharedPtr is that of std::shared_ptr.
 *
 * If QUESO does not detect C++11 functionality, boost is required.  In that
 * case, the underlying implementation is that of boost::shared_ptr.
 *
 * If QUESO does not detect C++11 functionality and does not detect boost then
 * QUESO will fail to configure.
 *
 * SharedPtr is implemented as a template typdef pattern.  That is, one
 * declares a SharedPtr like so:
 *
 * QUESO::SharedPtr< name_of_type >::Type name_of_variable;
 */
namespace QUESO
{
#ifdef QUESO_HAVE_CXX11_SHARED_PTR
  template<typename T>
  struct SharedPtr
  {
    typedef std::shared_ptr<T> Type;
  };
#elif QUESO_HAVE_BOOST_SHARED_PTR_HPP
  template<typename T>
  struct SharedPtr
  {
    typedef boost::shared_ptr<T> Type;
  };
#else
#     error "No valid definition for SharedPtr found!"
#endif

} // end namespace QUESO

#endif // UQ_SHARED_PTR_H
