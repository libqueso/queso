//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
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

#include <queso/Environment.h>
#include <queso/BaseInputOptions.h>

namespace QUESO {

BaseInputOptions::BaseInputOptions(const BaseEnvironment * env)
  :
    m_env(env),
    m_optionsDescription(new po::options_description("Input options")),
    m_optionsMap(new po::variables_map())
{
}

BaseInputOptions::BaseInputOptions()
  :
    m_env(NULL),
    m_optionsDescription(new po::options_description("Input options")),
    m_optionsMap(new po::variables_map())
{
}

BaseInputOptions::~BaseInputOptions()
{
  if (m_optionsDescription) {
    delete m_optionsDescription;
  }

  if (m_optionsMap) {
    delete m_optionsMap;
  }
}

void
BaseInputOptions::scanOptionsValues()
{
  queso_require_msg(m_optionsDescription, "m_optionsDescription variable is NULL");

  // If it's NULL then the defaults are used
  if (m_env != NULL) {
    defineOptions();
    (*m_env).scanInputFileForMyOptions(*m_optionsDescription);
    getOptionValues();
  }

  // if (m_env.subDisplayFile() != NULL) {
  //   *m_env.subDisplayFile() << "In BaseInputOptions::scanOptionsValues()"
  //                           << ": after reading values of options with prefix '" << m_prefix
  //                           << "', state of  object is:"
  //                           << "\n" << *this
  //                           << std::endl;
  // }
}

}  // End namespace QUESO
