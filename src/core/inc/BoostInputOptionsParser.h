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

#ifndef UQ_BOOST_INPUT_OPTIONS_H
#define UQ_BOOST_INPUT_OPTIONS_H

#include <string>
#include <queso/ScopedPtr.h>
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
   * The use-case that options are *not* read in from a file but are instead
   * wholly provided by the user.
   */
  BoostInputOptionsParser();

  //! Destructor
  /*!
   * Deletes m_optionsDescription
   */
  virtual ~BoostInputOptionsParser();

  //! This is the method that parses the input file
  /*!
   * It parses the input file and sets the m_optionsMap member.
   */
  void scanInputFile();

  //! Call this to register an option with the parser.
  /*!
   * The name of the option to look for in the input file is given by \c name
   * is given by.  If no option is present in the input file, then use
   * \c defaultValue is used for the default value for the option.  Describe
   * the option with a helpful message using \c description.
   */
  template <class T>
  void registerOption(const std::string & name,
                      const T & defaultValue,
                      const std::string & description);

  //! For flags *without* values.  Like a help message, for example.
  void registerOption(const std::string & name, const std::string & description);

  //! Get option \c name from the parser and set \c value to the parsed value.
  template <class T>
  void getOption(const std::string & name, T & value) const;

  //! Helpful stream operator for printing the parser state
  friend std::ostream & operator<<(std::ostream & os,
      const BoostInputOptionsParser & parser);

protected:
  // Needs to be a copy, not a reference; we don't want to force users
  // to keep their strings around and non-temporary, and we don't have
  // a non-temporary "" string for internal use.
  const std::string m_filename;

  ScopedPtr<boost::program_options::options_description>::Type m_optionsDescription;
  ScopedPtr<boost::program_options::variables_map>::Type m_optionsMap;

private:
  bool m_scannedInputFile;
};

}  // End namespace QUESO

#endif  // UQ_BOOST_INPUT_OPTIONS_H

#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
