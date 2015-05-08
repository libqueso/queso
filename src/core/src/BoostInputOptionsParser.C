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

#include <boost/program_options.hpp>

#include <queso/Environment.h>
#include <queso/BoostInputOptionsParser.h>

namespace QUESO {

BoostInputOptionsParser::BoostInputOptionsParser(const BaseEnvironment * env)
  :
    m_env(env),
    m_optionsDescription(new boost::program_options::options_description("Input options")),
    m_optionsMap(new boost::program_options::variables_map())
{
}

BoostInputOptionsParser::BoostInputOptionsParser()
  :
    m_env(NULL),
    m_optionsDescription(new boost::program_options::options_description("Input options")),
    m_optionsMap(new boost::program_options::variables_map())
{
}

BoostInputOptionsParser::~BoostInputOptionsParser()
{
  if (m_optionsDescription) {
    delete m_optionsDescription;
  }

  if (m_optionsMap) {
    delete m_optionsMap;
  }
}

void
BoostInputOptionsParser::scanOptionsValues()
{
  queso_require_msg(m_optionsDescription, "m_optionsDescription variable is NULL");

  // If it's NULL then the defaults are used
  if (m_env != NULL) {
    defineOptions();

    queso_require_not_equal_to_msg(m_env->optionsInputFileName(), "", "m_optionsInputFileName is 'nothing'");

    std::ifstream ifs;
    ifs.open(m_env->optionsInputFileName().c_str());

    queso_require_msg(m_optionsMap, "m_allOptionsMap variable is NULL");
    boost::program_options::store(boost::program_options::parse_config_file(ifs, *m_optionsDescription, true), *m_optionsMap);
    boost::program_options::notify(*m_optionsMap);

    ifs.close();

    getOptionValues();
  }

  // if (m_env.subDisplayFile() != NULL) {
  //   *m_env.subDisplayFile() << "In BoostInputOptionsParser::scanOptionsValues()"
  //                           << ": after reading values of options with prefix '" << m_prefix
  //                           << "', state of  object is:"
  //                           << "\n" << *this
  //                           << std::endl;
  // }
}

}  // End namespace QUESO
