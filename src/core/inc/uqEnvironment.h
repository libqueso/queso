//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010 The PECOS Development Team
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
// 
// $Id$
//
//--------------------------------------------------------------------------

#ifndef __UQ_ENVIRONMENT_H__
#define __UQ_ENVIRONMENT_H__

#include <uqDefines.h>
class uqEnvironmentOptionsClass;

#undef UQ_USES_COMMAND_LINE_OPTIONS

#include <hdf5.h>
#include <Epetra_MpiComm.h>
#include <gsl/gsl_rng.h>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <iostream>
#include <fstream>

extern unsigned long int gsl_rng_default_seed;

struct uqFilePtrSetStruct {
  uqFilePtrSetStruct();
 ~uqFilePtrSetStruct();

  std::ofstream* ofsVar;
  std::ifstream* ifsVar;
  hid_t          h5Var;
};

//------------------------------------------------------------------------
// Library versioning routines: we include them in a QUESO namespace
// here so that mutiple classes can use them as required (and so we
// can have a standalone versioning binary which does not require a full
// QUESO MPI environment.
//------------------------------------------------------------------------

namespace QUESO {
  void QUESO_version_print       (std::ostream &os);
  int  QUESO_get_numeric_version ();
}

//*****************************************************
// Base class
//*****************************************************

/*! This class sets up the environment underlying the use of the QUESO library by an executable. It:
<list type=number>
<item> assigns rank numbers, other than the world rank, to nodes participating in a parallel job,
<item> provides communicators for generating a sequence of vectors in a distributed way,
<item> provides functionality to read options from the 'options input file' (whose name is passed in the constructor of this environment class),
<item> opens output files for messages that would otherwise be written to the screen (one output file per allowed rank is opened and allowed ranks can be specified through the 'options input file').
</list>
*/
/*! -------------------------------------------------------------
*/
/*! This class is virtual. It is inherited by 'uqEmptyEnvironmentClass' and 'uqFullEnvironmentClass'.
    The QUESO environment class is instantiated at the application level, right after 'MPI_Init(&argc,&argv)'. 
    The QUESO environment is required by reference by many constructors in the QUESO library, and is available by reference from many classes as well.
*/
/*! -------------------------------------------------------------
*/
/*! Throughout QUESO, there are five classes whose constructors check options in the 'options input file':
<list type=number>
<item> uqBaseEnvironmentClass
<item> uqStatisticalInverseProblemClass
<item> uqStatisticalForwardProblemClass
<item> uqMetropolisHastingsSGClass ('SG' stands for 'sequence generator')
<item> uqMonteCarloSGClass
</list>
*/
/*! These classes rely on 'options classes' to read their options from the input file.
    The options classes are, respectively:
<list type=number>
<item> uqEnvironmentOptionsClass
<item> uqStatisticalInverseProblemOptionsClass
<item> uqStatisticalForwardProblemOptionsClass
<item> uqMetropolisHastingsSGOptionsClass
<item> uqMonteCarloSGOptionsClass
</list>
    The last two classes also rely on uqSequenceStatisticalOptionsClass for reading the
    options specifying which statistics have to be computed on the sequences of vectors
    involved.
*/
/*! -------------------------------------------------------------
*/
/*! The QUESO environment class manages five types of communicators. Let:
<list type=number>
<item> 'W >= 1' be the size of whole world communicator involved in a parallel run;
<item> 'N >= 1' be the size of the communicator passed to the QUESO environment constructor;
<item> 'S >= 1' be the number of statistical problems a QUESO environment will be handling at the same time, in parallel.
</list>
    Usually 'W'='N', but such equality is not necessary.
    The number 'S' is equal to the QUESO environment option 'm_numSubEnvironments', and is equal to 1 by default.
    The number 'N' must be a multiple of 'S', otherwise the QUESO class prints a fatal error message and MPI aborts.
    The five types of communicators that QUESO manages are referred to as:
<list type=number>
<item> world = MPI_WORLD_COMM, of size W;
<item> full = communicator passed to the constructor of uqBaseEnvironmentClass, of size N and usually equal to the world communicator;
<item> sub = communicator of size N/S that contains the number of MPI nodes necessary to solve a statistical inverse problem or a statistical forward problem.
<item> self = MPI_SELF_COMM, of size 1;
<item> inter0 = communicator of size S formed by all MPI nodes that have 'sub' rank 0 in their respective 'sub' communicators.
</list>
    So, any given node has potentially five different ranks. Of course, if the user is solving just one statistical problem with just one MPI node, then all ranks are equal to zero.
*/
/*! -------------------------------------------------------------
*/
/*! In the QUESO library terminology, one might refer to a QUESO "full" environment composed of 'S' QUESO "sub" environments.
    Each sub environment is assigned a "sub" id varying from 0 (zero) to S-1.
    Each sub environment is able to generate a statistical inverse problem and/or a statistical forward problem.
    That is, each sub environment is able to handle a "sub" Markov chain (a sequence) of vectors and/or
    a "sub" Monte Carlo sequence of output vectors.
    The "sub" sequences can be seen as forming a "unified" sequence in a distributed way.
    Indeed, the virtual class 'uqVectorSequenceClass' provides "sub" and "unified" statistical operations.

    A QUESO "sub" environment eventually prints messages to its own output file. In order for that to happen, the requirements are:
<list type=number>
<item> option 'm_subDisplayFileName', a string, must be different than the default value ".";
<item> option 'm_subDisplayAllowedSet', a set of sub ids, must contain the id of the sub environment wanting to write a message to the output file;
<item> the previous requirement is automatically satisfied if the option 'm_subDisplayAllowAll', a boolean, is set to 1 (the default value is 0);
<item> the processor wanting to write a message to the output file must have sub rank 0 (zero).
</list>
    If all requirements are satisfied, then QUESO will generate a file with name '\<m_subDisplayFileName\>_sub\<sub id\>.txt'.
    For instance, if 'm_subDisplayFileName' is 'pROblem_775_' then a node of sub rank 0 in sub environment 17
    will write a message to the file 'pROblem_775_sub17.txt'.
*/
class uqBaseEnvironmentClass {
public:
  uqBaseEnvironmentClass(MPI_Comm inputComm, const char* passedOptionsInputFileName, const uqEnvOptionsValuesClass* alternativeOptionsValues);
  uqBaseEnvironmentClass(const uqBaseEnvironmentClass& obj);
  virtual ~uqBaseEnvironmentClass();

          uqBaseEnvironmentClass& operator=                  (const uqBaseEnvironmentClass& rhs);

          int                     worldRank                  () const;

          int                     fullRank                   () const;
          const Epetra_MpiComm&   fullComm                   () const; 

          int                     subRank                    () const;
          const Epetra_MpiComm&   subComm                    () const; 

          const Epetra_MpiComm&   selfComm                   () const; 

          int                     inter0Rank                 () const;
          const Epetra_MpiComm&   inter0Comm                 () const;

                std::ofstream*    subDisplayFile             () const;

          unsigned int            numSubEnvironments         () const;
          unsigned int            subId                      () const;
          const std::string&      subIdString                () const;
          void                    checkTheParallelEnvironment() const;

          std::string             optionsInputFileName       () const;
#ifdef UQ_USES_COMMAND_LINE_OPTIONS
          const po::options_description& allOptionsDesc      () const;
#endif
          po::variables_map&      allOptionsMap              () const;
          void                    scanInputFileForMyOptions  (const po::options_description& optionsDesc) const;
          unsigned int            displayVerbosity           () const;
          unsigned int            syncVerbosity              () const;
          const gsl_rng*          rng                        () const;
          int                     seed                       () const;
          void                    resetGslSeed               (int newSeedOption);
	  std::string             identifyingString          () const;
          void                    resetIdentifyingString     (const std::string& newString) const; // Yes, const
          bool                    isThereInputFile           () const;
          void                    syncPrintDebugMsg          (const char* msg, unsigned int msgVerbosity, unsigned int numUSecs, const Epetra_MpiComm& commObj) const;

          bool                    openOutputFile             (const std::string&            fileName,
                                                              const std::string&            fileType,
                                                              const std::set<unsigned int>& allowedSubEnvIds,
                                                                    bool                    writeOver,
                                                                    uqFilePtrSetStruct&     filePtrSet) const;

          bool                    openUnifiedOutputFile      (const std::string&            fileName,
                                                              const std::string&            fileType,
                                                                    bool                    writeOver,
                                                                    uqFilePtrSetStruct&     filePtrSet) const;
          bool                    openInputFile              (const std::string&            fileName,
                                                              const std::string&            fileType,
                                                              const std::set<unsigned int>& allowedSubEnvIds,
                                                                    uqFilePtrSetStruct&     filePtrSet) const;
          bool                    openUnifiedInputFile       (const std::string&            fileName,
                                                              const std::string&            fileType,
                                                                    uqFilePtrSetStruct&     filePtrSet) const;
          void                    closeFile                  (      uqFilePtrSetStruct&     filePtrSet,
                                                              const std::string&            fileType) const; 
          void                    setExceptionalCircunstance (bool value) const;
          bool                    exceptionalCircunstance    () const;
  virtual void                    print                      (std::ostream& os) const = 0;

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
  mutable bool               m_exceptionalCircunstance;

  uqEnvOptionsValuesClass    m_alternativeOptionsValues;
  uqEnvironmentOptionsClass* m_optionsObj;
};

//*****************************************************
// Empty Environment
//*****************************************************
class uqEmptyEnvironmentClass : public uqBaseEnvironmentClass {
public:
  uqEmptyEnvironmentClass();
 ~uqEmptyEnvironmentClass();

        void print(std::ostream& os) const;
};

//*****************************************************
// Full Environment
//*****************************************************
class uqFullEnvironmentClass : public uqBaseEnvironmentClass {
public:
  uqFullEnvironmentClass(MPI_Comm inputComm, const char* passedOptionsInputFileName, const char* prefix, const uqEnvOptionsValuesClass* alternativeOptionsValues);
 ~uqFullEnvironmentClass();

        void        print               (std::ostream& os) const;

private:
        void        readOptionsInputFile();
};

std::ostream& operator<<(std::ostream& os, const uqBaseEnvironmentClass& obj);

#endif // __UQ_ENVIRONMENT_H__
