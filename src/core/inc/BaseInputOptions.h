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

#ifndef UQ_BASE_INPUT_OPTIONS_H
#define UQ_BASE_INPUT_OPTIONS_H

namespace boost
{
  namespace program_options
  {
    class options_description;
  }
}

namespace QUESO {

class BaseEnvironment;

class BaseInputOptions
{
public:
  BaseInputOptions(const BaseEnvironment * env);
  BaseInputOptions();
  virtual ~BaseInputOptions();

  void scanOptionsValues();

private:
  virtual void defineOptions() = 0;
  virtual void getOptionValues() = 0;

protected:
  const BaseEnvironment * m_env;
  boost::program_options::options_description * m_optionsDescription;
};

}  // End namespace QUESO

#endif // UQ_BASE_INPUT_OPTIONS_H
