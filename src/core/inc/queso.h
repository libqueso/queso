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
// grvy.h: Basic API Definitions
//
//--------------------------------------------------------------------------

#ifndef QUESO_H_
#define QUESO_H_

#warning Header queso.h is deprecated; use config_queso.h instead.

#include "config_queso.h"

// Library version/build info

// Deprecated backward-compatible duplicate definitions, now derived
// from config_queso.h

#ifndef QUESO_LIB_VERSION
#define QUESO_LIB_VERSION QUESO_VERSION
#warning QUESO_LIB_VERSION is deprecated; use QUESO_VERSION instead.
#endif

#ifndef QUESO_LIB_RELEASE
#define QUESO_LIB_RELEASE QUESO_BUILD_DEVSTATUS
#warning QUESO_LIB_RELEASE is deprecated; use QUESO_BUILD_DEVSTATUS instead.
#endif

#endif
