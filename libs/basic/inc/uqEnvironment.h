/* libs/basic/inc/uqEnvironment.h
 * 
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_ENVIRONMENT_H__
#define __UQ_ENVIRONMENT_H__

#undef UQ_USES_COMMAND_LINE_OPTIONS

#define UQ_ENV_VERBOSITY_ODV        0
#define UQ_ENV_SEED_ODV             0
#define UQ_ENV_NUM_DEBUG_PARAMS_ODV 0
#define UQ_ENV_DEBUG_PARAM_ODV      0.

#include <Epetra_MpiComm.h>
#include <gsl/gsl_rng.h>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <iostream>
#include <fstream>

extern unsigned long int gsl_rng_default_seed;

struct uqEnvOptionsStruct {
  uqEnvOptionsStruct(unsigned int verbosity,
                     int          seed);
 ~uqEnvOptionsStruct();

  unsigned int        m_verbosity;
  int                 m_seed;
  unsigned int        m_numDebugParams;
  std::vector<double> m_debugParams;
};

//*****************************************************
// Base class
//*****************************************************
class uqBaseEnvironmentClass {
public:
  uqBaseEnvironmentClass();
  uqBaseEnvironmentClass(int& argc, char** &argv);
  uqBaseEnvironmentClass(const uqEnvOptionsStruct& options);
  uqBaseEnvironmentClass(const uqBaseEnvironmentClass& obj);
  virtual ~uqBaseEnvironmentClass();

          uqBaseEnvironmentClass& operator=                (const uqBaseEnvironmentClass& rhs);
          int                     rank                     () const;
          void                    barrier                  () const;
          const Epetra_MpiComm&   comm                     () const; 
#ifdef UQ_USES_COMMAND_LINE_OPTIONS
  const po::options_description&  allOptionsDesc           () const;
#endif
          po::variables_map&      allOptionsMap            () const;
          void                    scanInputFileForMyOptions(const po::options_description& optionsDesc) const;
          unsigned int            verbosity                () const;
          const gsl_rng*          rng                      () const;
          bool                    isThereInputFile         () const;
  virtual void                    print                    (std::ostream& os) const = 0;

protected:
  int                      m_argc;
  char**                   m_argv;
  Epetra_MpiComm*          m_comm;
  int                      m_rank;
  int                      m_commSize;
  bool                     m_argsWereProvided;
  bool                     m_thereIsInputFile;
  std::string              m_inputFileName;
  po::options_description* m_allOptionsDesc;
  po::options_description* m_envOptionsDesc;
  po::variables_map*       m_allOptionsMap;
  unsigned int             m_verbosity;
  int                      m_seed;
  unsigned int             m_numDebugParams;
  std::vector<double>      m_debugParams;
  gsl_rng*                 m_rng;
  struct timeval           m_timevalBegin;
};

//*****************************************************
// Empty Environment
//*****************************************************
class uqEmptyEnvironmentClass : public uqBaseEnvironmentClass {
public:
  uqEmptyEnvironmentClass();
 ~uqEmptyEnvironmentClass();

        void                     print                    (std::ostream& os) const;
};

//*****************************************************
// Full Environment
//*****************************************************
class uqEnvironmentClass : public uqBaseEnvironmentClass {
public:
  uqEnvironmentClass();
  uqEnvironmentClass(int& argc, char** &argv);
  uqEnvironmentClass(const uqEnvOptionsStruct& options);
 ~uqEnvironmentClass();

        void                     print                    (std::ostream& os) const;

private:
        void                     commonConstructor        ();
        void                     readEventualInputFile    ();
        void                     defineMyOptions          (po::options_description& optionsDesc) const;
        void                     getMyOptionValues        (po::options_description& optionsDesc);
};

std::ostream& operator<<(std::ostream& os, const uqBaseEnvironmentClass& obj);

#endif // __UQ_ENVIRONMENT_H__
