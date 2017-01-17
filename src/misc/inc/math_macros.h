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

#ifndef UQ_MATH_MACROS_H
#define UQ_MATH_MACROS_H

#include <queso/config_queso.h>

#ifdef QUESO_HAVE_CXX11_ISNAN
#include <cmath>
#elif QUESO_HAVE_BOOST_MATH_SPECIAL_FUNCTIONS_HPP
#include <boost/math/special_functions.hpp>
#endif

namespace QUESO
{
  template<typename T>
  bool queso_isnan( T arg )
  {
#ifdef QUESO_HAVE_CXX11_ISNAN
    return std::isnan(arg);
#elif QUESO_HAVE_BOOST_MATH_SPECIAL_FUNCTIONS_HPP
    return (boost::math::isnan)(arg);
#else
#     error "No valid definition for is_nan found!"
#endif
  }

  template<typename T>
  bool queso_isfinite( T arg )
  {
#ifdef QUESO_HAVE_CXX11_ISFINITE
    return std::isfinite(arg);
#elif QUESO_HAVE_BOOST_MATH_SPECIAL_FUNCTIONS_HPP
    return (boost::math::isfinite)(arg);
#else
#     error "No valid definition for isfinite found!"
#endif
  }

} // end namespace QUESO

#endif // UQ_MATH_MACROS_H
