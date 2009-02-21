/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * $Id$
 *
 * Brief description of this file: 
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __UQ_ENVIRONMENT_H__
#define __UQ_ENVIRONMENT_H__

#undef UQ_USES_COMMAND_LINE_OPTIONS

#define UQ_ENV_NUM_APPL_INSTANCES_ODV 1
#define UQ_ENV_VERBOSITY_ODV           0
#define UQ_ENV_SEED_ODV                0
#define UQ_ENV_RUN_NAME_ODV            ""
#define UQ_ENV_NUM_DEBUG_PARAMS_ODV    0
#define UQ_ENV_DEBUG_PARAM_ODV         0.

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

  unsigned int        m_numApplInstances;
  unsigned int        m_verbosity;
  int                 m_seed;
  std::string         m_runName;
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
          int                     myApplRank               () const;
          const Epetra_MpiComm&   worldComm                () const; 
          const Epetra_MpiComm&   myApplComm               () const; 

#ifdef UQ_USES_COMMAND_LINE_OPTIONS
          const po::options_description& allOptionsDesc    () const;
#endif
          po::variables_map&      allOptionsMap            () const;
          void                    scanInputFileForMyOptions(const po::options_description& optionsDesc) const;
          unsigned int            verbosity                () const;
          const std::string&      runName                  () const;
          const gsl_rng*          rng                      () const;
          bool                    isThereInputFile         () const;
  virtual void                    print                    (std::ostream& os) const = 0;

protected:
  int                      m_argc;
  char**                   m_argv;

  Epetra_MpiComm*          m_worldComm;
  int                      m_worldRank;
  int                      m_worldCommSize;
  MPI_Group                m_worldGroup;

  bool                     m_argsWereProvided;
  bool                     m_thereIsInputFile;
  std::string              m_inputFileName;
  po::options_description* m_allOptionsDesc;
  po::options_description* m_envOptionsDesc;
  po::variables_map*       m_allOptionsMap;
  unsigned int             m_numApplInstances;
  unsigned int             m_verbosity;
  int                      m_seed;
  std::string              m_runName;
  unsigned int             m_numDebugParams;
  std::vector<double>      m_debugParams;

  unsigned int             m_myApplInstanceId;
  MPI_Group                m_myApplGroup;
  MPI_Comm                 m_myApplRawComm;
  Epetra_MpiComm*          m_myApplComm;
  int                      m_myApplRank;
  int                      m_myApplCommSize;

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
class uqFullEnvironmentClass : public uqBaseEnvironmentClass {
public:
  uqFullEnvironmentClass();
  uqFullEnvironmentClass(int& argc, char** &argv);
  uqFullEnvironmentClass(const uqEnvOptionsStruct& options);
 ~uqFullEnvironmentClass();

        void                     print                    (std::ostream& os) const;

private:
        void                     commonConstructor        ();
        void                     readEventualInputFile    ();
        void                     defineMyOptions          (po::options_description& optionsDesc) const;
        void                     getMyOptionValues        (po::options_description& optionsDesc);
};

std::ostream& operator<<(std::ostream& os, const uqBaseEnvironmentClass& obj);

#endif // __UQ_ENVIRONMENT_H__
