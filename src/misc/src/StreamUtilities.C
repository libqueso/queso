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

// This class
#include <queso/asserts.h>
#include <queso/StreamUtilities.h>

// C++
#include <istream>
#include <locale> // isspace

namespace QUESO
{
  void StreamUtilities::skip_comment_lines( std::istream &in,
                                            const char comment_start)
  {
    char c, line[256];

    in.get(c);
    // std::isspace used to be std::isblank, but the latter requires C++11
    // There are slight differences, but for our purposes, this should
    // be OK.
    while(std::isspace(c))in.get(c);

    in.putback(c);
    queso_require( !in.fail() );

    while (in.get(c), c==comment_start)
       in.getline (line, 255);

    // put back first character of
    // first non-comment line
    in.putback (c);
    queso_require( !in.fail() );

  }

} // end namespace QUESO
