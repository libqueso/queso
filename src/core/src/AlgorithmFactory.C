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

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/AlgorithmFactory.h>

namespace QUESO
{

template <>
std::map<std::string, Factory<Algorithm<GslVector, GslMatrix> > *> &
Factory<Algorithm<GslVector, GslMatrix> >::factory_map()
{
  static std::map<std::string, Factory<Algorithm<GslVector, GslMatrix> > *> _factory_map;

  return _factory_map;
}

const BaseEnvironment * AlgorithmFactory::m_env = NULL;
const BaseTKGroup<GslVector, GslMatrix> * AlgorithmFactory::m_tk = NULL;

} // namespace QUESO
