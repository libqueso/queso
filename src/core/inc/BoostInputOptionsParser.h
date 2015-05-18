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

#ifndef UQ_BOOST_INPUT_OPTIONS_H
#define UQ_BOOST_INPUT_OPTIONS_H

#include <string>
#include <queso/BaseInputOptionsParser.h>

namespace boost
{
  namespace program_options
  {
    class options_description;
    class variables_map;
  }
}

namespace QUESO {

class BoostInputOptionsParser : public BaseInputOptionsParser
{
public:
  //! Constructor that parses file \c filename
  BoostInputOptionsParser(const std::string & filename);

  //! Default constructor that sets m_filename to ""
  /*!
   * The use-case for m_env being NULL is that options are *not* read in from
   * a file but are instead wholly provided by the user.
   */
  BoostInputOptionsParser();

  //! Destructor
  /*!
   * Deletes m_optionsDescription
   */
  virtual ~BoostInputOptionsParser();

  //! This is the method that parses the input file
  /*!
   * It calls defineOptions, which sets m_optionsDescription, to define the
   * boost options, then it parses the input file and sets the m_optionsMap
   * member.  After both of those are done, it *sets* the options by calling
   * getOptionValues.
   */
  void scanInputFile();

  template <class T>
  void registerOption(std::string name, T defaultValue, std::string description);

  //! For flags *without* values
  void registerOption(std::string name, std::string description);

  template <class T>
  void getOption(std::string & name, T & value);

protected:
  const std::string & m_filename;
  boost::program_options::options_description * m_optionsDescription;
  boost::program_options::variables_map * m_optionsMap;

private:
  bool m_scannedInputFile;
};

}  // End namespace QUESO

#endif  // UQ_BOOST_INPUT_OPTIONS_H
