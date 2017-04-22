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

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>

#include <queso/BoostInputOptionsParser.h>
#include <queso/Miscellaneous.h>
#include <queso/Defines.h>

namespace QUESO {

BoostInputOptionsParser::BoostInputOptionsParser(const std::string & filename)
  :
    m_filename(filename),
    m_optionsDescription(new boost::program_options::options_description("Input options")),
    m_optionsMap(new boost::program_options::variables_map()),
    m_scannedInputFile(false)
{
  queso_deprecated();
}

BoostInputOptionsParser::BoostInputOptionsParser()
  :
    m_filename(""),
    m_optionsDescription(new boost::program_options::options_description("Input options")),
    m_optionsMap(new boost::program_options::variables_map()),
    m_scannedInputFile(false)
{
  queso_deprecated();
}

BoostInputOptionsParser::~BoostInputOptionsParser()
{
  queso_deprecated();
  // Do nothing
}

void
BoostInputOptionsParser::scanInputFile()
{
  queso_deprecated();
  queso_require_msg(m_optionsDescription, "m_optionsDescription variable is NULL");

  // If it's the empty string then the defaults are used
  if (m_filename != "") {
    std::ifstream ifs;
    ifs.open(m_filename.c_str());

    queso_require_msg(m_optionsMap, "m_allOptionsMap variable is NULL");
    boost::program_options::store(
        boost::program_options::parse_config_file(
          ifs, *m_optionsDescription, true), *m_optionsMap);
    boost::program_options::notify(*m_optionsMap);

    ifs.close();

    m_scannedInputFile = true;
  }
}

template <typename T>
void
BoostInputOptionsParser::registerOption(const std::string & name,
                                        const T & defaultValue,
                                        const std::string & description)
{
  queso_deprecated();
  m_optionsDescription->add_options()
    (name.c_str(),
     boost::program_options::value<T>()->default_value(defaultValue),
     description.c_str());
}

void
BoostInputOptionsParser::registerOption(const std::string & name,
                                        const std::string & description)
{
  queso_deprecated();
  m_optionsDescription->add_options()
    (name.c_str(),
     description.c_str());
}

template <typename T>
void
BoostInputOptionsParser::getOption(const std::string & name, T & value) const
{
  queso_deprecated();
  if (m_scannedInputFile) {
    value = (*m_optionsMap)[name].as<T>();
  }
}

template <>
void
BoostInputOptionsParser::getOption(const std::string & name, std::set<unsigned int, std::less<unsigned int>, std::allocator<unsigned int> > & value) const
{
  queso_deprecated();
  if (m_scannedInputFile) {
    // Clear before putting things in it
    value.clear();

    // Get the option as a string
    // DM:  Why do it this way?  Doesn't boost support vectors as input
    //      options?
    std::vector<double> tmpVec(0, 0.0);
    std::string optionValue;
    this->getOption<std::string>(name, optionValue);
    MiscReadDoublesFromString(optionValue, tmpVec);

    for (unsigned int i = 0; i < tmpVec.size(); i++) {
      value.insert((unsigned int) tmpVec[i]);  // Why cast?!
    }
  }
}

template <>
void
BoostInputOptionsParser::getOption<std::vector<double, std::allocator<double> > >(const std::string & name, std::vector<double, std::allocator<double> > & value) const
{
  queso_deprecated();
  if (m_scannedInputFile) {
    // Need to reset value?

    // Get the option as a string
    // DM:  Why do it this way?  Doesn't boost support vectors as input options?
    std::string optionValue;
    this->getOption<std::string>(name, optionValue);
    MiscReadDoublesFromString(optionValue, value);
  }
}

std::ostream &
operator<<(std::ostream & os, const BoostInputOptionsParser & parser)
{
  queso_deprecated();
  os << *(parser.m_optionsDescription);
  return os;
}

template void BoostInputOptionsParser::registerOption<int>(const std::string &, const int &, const std::string &);
template void BoostInputOptionsParser::registerOption<unsigned int>(const std::string &, const unsigned int &, const std::string &);
template void BoostInputOptionsParser::registerOption<bool>(const std::string &, const bool &, const std::string &);
template void BoostInputOptionsParser::registerOption<double>(const std::string &, const double &, const std::string &);
template void BoostInputOptionsParser::registerOption<std::string>(const std::string &, const std::string &, const std::string &);

template void BoostInputOptionsParser::getOption<int>(const std::string&, int&) const;
template void BoostInputOptionsParser::getOption<unsigned int>(const std::string&, unsigned int&) const;
template void BoostInputOptionsParser::getOption<double>(const std::string&, double&) const;
template void BoostInputOptionsParser::getOption<bool>(const std::string&, bool&) const;
template void BoostInputOptionsParser::getOption<std::string>(const std::string&, std::string&) const;
template void BoostInputOptionsParser::getOption<std::set<unsigned int, std::less<unsigned int>, std::allocator<unsigned int> > >(const std::string&, std::set<unsigned int, std::less<unsigned int>, std::allocator<unsigned int> >&) const;
template void BoostInputOptionsParser::getOption<std::vector<unsigned int, std::allocator<unsigned int> > >(const std::string&, std::vector<unsigned int, std::allocator<unsigned int> >&) const;
template void BoostInputOptionsParser::getOption<std::vector<double, std::allocator<double> > >(const std::string&, std::vector<double, std::allocator<double> >&) const;

}  // End namespace QUESO
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
