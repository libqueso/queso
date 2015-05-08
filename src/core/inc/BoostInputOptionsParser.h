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

namespace boost
{
  namespace program_options
  {
    class options_description;
    class variables_map;
  }
}

namespace QUESO {

class BaseEnvironment;

class BoostInputOptionsParser
{
public:
  //! Constructor that sets the internal environment
  BoostInputOptionsParser(const BaseEnvironment * env);

  //! Default constructor that sets m_env to NULL
  /*!
   * The use-case for m_env being NULL is that options are *not* read in from
   * a file but are instead wholly provided by the user.
   */
  BoostInputOptionsParser();

  //! Destructor
  /*!
   * Deletes m_optionsDescription, but not m_env
   */
  virtual ~BoostInputOptionsParser();

  //! Calls the relevant QUESO BaseEnvironment methods to scan the input file
  void scanOptionsValues();

private:
  //! Subclasses implement this to *define* input options
  /*!
   * This will act on the boost-specific m_optionsDescription
   */
  virtual void defineOptions() = 0;

  //! Subclasses implement this to *set* input options from an input string
  /*!
   * This will act on the boost-specific m_optionsDescription
   */
  virtual void getOptionValues() = 0;

protected:
  const BaseEnvironment * m_env;
  boost::program_options::options_description * m_optionsDescription;
  boost::program_options::variables_map * m_optionsMap;
};

}  // End namespace QUESO

#endif  // UQ_BOOST_INPUT_OPTIONS_H
