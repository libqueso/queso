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

#include <uqDefines.h>

#undef UQ_USES_COMMAND_LINE_OPTIONS

#define UQ_ENV_FILENAME_FOR_NO_OUTPUT_FILE "."

#define UQ_ENV_NUM_SUB_ENVIRONMENTS_ODV        1
#define UQ_ENV_SUB_SCREEN_WRITE_ODV            0
#define UQ_ENV_SUB_SCREEN_OUTPUT_FILE_NAME_ODV UQ_ENV_FILENAME_FOR_NO_OUTPUT_FILE
#define UQ_ENV_VERBOSITY_ODV                   0
#define UQ_ENV_SEED_ODV                        0
#define UQ_ENV_NUM_DEBUG_PARAMS_ODV            0
#define UQ_ENV_DEBUG_PARAM_ODV                 0.

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

  unsigned int        m_numSubEnvironments;
  bool                m_subScreenWrite;
  std::string         m_subScreenOutputFileName;
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
  uqBaseEnvironmentClass(MPI_Comm inputComm, const char* prefix);
  uqBaseEnvironmentClass(int& argc, char** &argv, MPI_Comm inputComm, const char* prefix);
  uqBaseEnvironmentClass(const uqEnvOptionsStruct& options, MPI_Comm inputComm, const char* prefix);
  uqBaseEnvironmentClass(const uqBaseEnvironmentClass& obj);
  virtual ~uqBaseEnvironmentClass();

          uqBaseEnvironmentClass& operator=                (const uqBaseEnvironmentClass& rhs);

          int                     rank                     () const;
          const Epetra_MpiComm&   fullComm                 () const; 

          int                     subRank                  () const;
          const Epetra_MpiComm&   subComm                  () const; 

          const Epetra_MpiComm&   selfComm                 () const; 

                std::ofstream*    subScreenFile            () const;

          unsigned int            numSubEnvironments       () const;
          unsigned int            subId                    () const;
          const std::string&      subIdString              () const;

#ifdef UQ_USES_COMMAND_LINE_OPTIONS
          const po::options_description& allOptionsDesc    () const;
#endif
          po::variables_map&      allOptionsMap            () const;
          void                    scanInputFileForMyOptions(const po::options_description& optionsDesc) const;
          unsigned int            verbosity                () const;
          const gsl_rng*          rng                      () const;
          bool                    isThereInputFile         () const;
          void                    syncPrintDebugMsg        (const char* msg, unsigned int numUSecs, const Epetra_MpiComm& commObj) const;
  virtual void                    print                    (std::ostream& os) const = 0;

protected:
  int                      m_argc;
  char**                   m_argv;
  std::string              m_prefix;

  Epetra_MpiComm*          m_fullComm;
  int                      m_fullRank;
  int                      m_fullCommSize;
  MPI_Group                m_fullGroup;

  bool                     m_argsWereProvided;
  bool                     m_thereIsInputFile;
  std::string              m_inputFileName;
  po::options_description* m_allOptionsDesc;
  po::options_description* m_envOptionsDesc;
  po::variables_map*       m_allOptionsMap;

  std::string              m_option_help;
  std::string              m_option_numSubEnvironments;
  std::string              m_option_subScreenWrite;
  std::string              m_option_subScreenOutputFileName;
  std::string              m_option_verbosity;
  std::string              m_option_seed;

  unsigned int             m_numSubEnvironments;
  bool                     m_subScreenWrite;
  std::string              m_subScreenOutputFileName;
  unsigned int             m_verbosity;
  int                      m_seed;

  unsigned int             m_numDebugParams;
  std::vector<double>      m_debugParams;

  unsigned int             m_subId;
  std::string              m_subIdString;
  MPI_Group                m_subGroup;
  MPI_Comm                 m_subRawComm;
  Epetra_MpiComm*          m_subComm;
  int                      m_subRank;
  int                      m_subCommSize;

  Epetra_MpiComm*          m_selfComm;

  mutable std::ofstream*   m_subScreenFile;
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

        void                     print                (std::ostream& os) const;
};

//*****************************************************
// Full Environment
//*****************************************************
class uqFullEnvironmentClass : public uqBaseEnvironmentClass {
public:
  uqFullEnvironmentClass(MPI_Comm inputComm, const char* prefix);
  uqFullEnvironmentClass(int& argc, char** &argv, MPI_Comm inputComm, const char* prefix);
  uqFullEnvironmentClass(const uqEnvOptionsStruct& options, MPI_Comm inputComm, const char* prefix);
 ~uqFullEnvironmentClass();

        void                     print                (std::ostream& os) const;

private:
        void                     commonConstructor    (MPI_Comm inputComm);
        void                     readEventualInputFile();
        void                     defineMyOptions      (po::options_description& optionsDesc) const;
        void                     getMyOptionValues    (po::options_description& optionsDesc);
};

std::ostream& operator<<(std::ostream& os, const uqBaseEnvironmentClass& obj);

#endif // __UQ_ENVIRONMENT_H__
