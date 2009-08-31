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
class uqEnvironmentOptionsClass;

#undef UQ_USES_COMMAND_LINE_OPTIONS

#include <Epetra_MpiComm.h>
#include <gsl/gsl_rng.h>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <iostream>
#include <fstream>

extern unsigned long int gsl_rng_default_seed;

//*****************************************************
// Base class
//*****************************************************

/*! A class that sets up the environment underlying a QUESO run:
    - assigns rank numbers, other than the world rank, to nodes participating in a parallel job;
    - open output files for messages that would otherwise be written to the screen; one output file per allowed rank is opened; allowed ranks are specified in the options input file; */
/*! */
/*! This class is virtual. It is inherited by 'uqEmptyEnvironmentClass' and 'uqFullEnvironmentClass'.
    The environment options can be set in an input file. The name of this input file is one of the parameters in the constructor. */
/*! The uqEnvironmentClass internally delegates the task of reading input options to uqEnvironmentOptionsClass.                   */
class uqBaseEnvironmentClass {
public:
  uqBaseEnvironmentClass(MPI_Comm inputComm, const char* inputFileName);
  uqBaseEnvironmentClass(const uqBaseEnvironmentClass& obj);
  virtual ~uqBaseEnvironmentClass();

          uqBaseEnvironmentClass& operator=                (const uqBaseEnvironmentClass& rhs);

          int                     worldRank                () const;

          int                     fullRank                 () const;
          const Epetra_MpiComm&   fullComm                 () const; 

          int                     subRank                  () const;
          const Epetra_MpiComm&   subComm                  () const; 

          const Epetra_MpiComm&   selfComm                 () const; 

          int                     inter0Rank               () const;
          const Epetra_MpiComm&   inter0Comm               () const; 

                std::ofstream*    subDisplayFile           () const;

          unsigned int            numSubEnvironments       () const;
          unsigned int            subId                    () const;
          const std::string&      subIdString              () const;

#ifdef UQ_USES_COMMAND_LINE_OPTIONS
          const po::options_description& allOptionsDesc    () const;
#endif
          po::variables_map&      allOptionsMap            () const;
          void                    scanInputFileForMyOptions(const po::options_description& optionsDesc) const;
          unsigned int            displayVerbosity         () const;
          unsigned int            syncVerbosity            () const;
          const gsl_rng*          rng                      () const;
          bool                    isThereInputFile         () const;
          void                    syncPrintDebugMsg        (const char* msg, unsigned int msgVerbosity, unsigned int numUSecs, const Epetra_MpiComm& commObj) const;

          void                    openOutputFile           (const std::string&            fileName,
                                                            const std::string&            fileType,
                                                            const std::set<unsigned int>& allowedSubEnvIds,
                                                                  bool                    writeOver,
                                                                  std::ofstream*&         ofsvar) const;

          void                    openUnifiedOutputFile    (const std::string&            fileName,
                                                            const std::string&            fileType,
                                                                  bool                    writeOver,
                                                                  std::ofstream*&         ofsvar) const;

  virtual void                    print                    (std::ostream& os) const = 0;

protected:
  int                        m_worldRank;

  MPI_Comm                   m_fullRawComm;
  Epetra_MpiComm*            m_fullComm;
  int                        m_fullRank;
  int                        m_fullCommSize;
  MPI_Group                  m_fullGroup;

  std::string                m_optionsInputFileName;
  po::options_description*   m_allOptionsDesc;
  po::variables_map*         m_allOptionsMap;

  unsigned int               m_subId;
  std::string                m_subIdString;
  MPI_Group                  m_subGroup;
  MPI_Comm                   m_subRawComm;
  Epetra_MpiComm*            m_subComm;
  int                        m_subRank;
  int                        m_subCommSize;

  Epetra_MpiComm*            m_selfComm;

  MPI_Group                  m_inter0Group;
  MPI_Comm                   m_inter0RawComm;
  Epetra_MpiComm*            m_inter0Comm;
  int                        m_inter0Rank;
  int                        m_inter0CommSize;

  mutable std::ofstream*     m_subDisplayFile;
  gsl_rng*                   m_rng;
  struct timeval             m_timevalBegin;

  uqEnvironmentOptionsClass* m_options;
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
  uqFullEnvironmentClass(MPI_Comm inputComm, const char* optionsInputFileName, const char* prefix);
 ~uqFullEnvironmentClass();

        void                     print                (std::ostream& os) const;

private:
        void                     readOptionsInputFile ();
};

std::ostream& operator<<(std::ostream& os, const uqBaseEnvironmentClass& obj);

#endif // __UQ_ENVIRONMENT_H__
