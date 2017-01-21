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

#include "config_queso.h"

#ifdef QUESO_HAVE_CPPUNIT

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

#include <queso/EnvironmentOptions.h>
#include <queso/BoxSubset.h>
#include <queso/MultiDQuadratureBase.h>

namespace QUESOTesting
{
  class LegendreQuadratureTestingHelper
  {
  public:

    void testing_orders( std::vector<unsigned int> & orders )
    {
      // These are the valid listed orders for Legendre quadrature in QUESO;
      // TODO: With C++11, we can initialize this with array syntax
      orders.resize(11);
      orders[0] = 1;
      orders[1] = 2;
      orders[2] = 3;
      orders[3] = 4;
      orders[4] = 5;
      orders[5] = 6;
      orders[6] = 7;
      orders[7] = 10;
      orders[8] = 11;
      orders[9] = 12;
      orders[10] = 16;
    }

  };

} // end namespace QUESOTesting

#endif // QUESO_HAVE_CPPUNIT
