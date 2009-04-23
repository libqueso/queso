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

#include <uqEnvironment.h>
#include <uqMiscellaneous.h>
#include <sys/time.h>
#include <gsl/gsl_randist.h>

// Version "0.1"   on "Aug/11/2008"
// Version "0.11"  on "Aug/15/2008"
// Version "0.2"   on "Oct/01/2008"
// Version "0.21"  on "Oct/08/2008"
// Version "0.3.0" on "Feb/13/2009"
// Version "0.3.1" on "Feb/19/2009"
// Version "0.4.0" on "MMM/DD/2009"
#define QUESO_LIBRARY_CURRENT_VERSION "0.4.0"
#define QUESO_LIBRARY_RELEASE_DATE    "MMM/DD/2009"

uqEnvOptionsStruct::uqEnvOptionsStruct(
  unsigned int verbosity,
  int          seed)
  :
  m_numSubEnvironments     (UQ_ENV_NUM_SUB_ENVIRONMENTS_ODV),
  m_subScreenOutputFileName(UQ_ENV_SUB_SCREEN_OUTPUT_FILE_NAME_ODV),
//m_subScreenOutputFilter  (),
  m_verbosity              (verbosity),
  m_syncVerbosity          (UQ_ENV_SYNC_VERBOSITY_ODV),
  m_seed                   (seed),
  m_numDebugParams         (0),
  m_debugParams            (0,0.)
{
}

uqEnvOptionsStruct::~uqEnvOptionsStruct()
{
}

//*****************************************************
// Base class
//*****************************************************
uqBaseEnvironmentClass::uqBaseEnvironmentClass(
  MPI_Comm    inputComm,
  const char* prefix)
  :
  m_argc                          (0),
  m_argv                          (NULL),
  m_prefix                        ((std::string)(prefix) + "env_"),
  m_worldRank                     (-1),
  m_fullRawComm                   (inputComm),
  m_fullComm                      (NULL),
  m_fullRank                      (-1),
  m_fullCommSize                  (1),
  m_argsWereProvided              (false),
  m_thereIsInputFile              (false),
  m_inputFileName                 (""),
  m_allOptionsDesc                (NULL),
  m_envOptionsDesc                (NULL),
  m_allOptionsMap                 (NULL),
  m_option_help                   (m_prefix + "help"                   ),
  m_option_numSubEnvironments     (m_prefix + "numSubEnvironments"     ),
  m_option_subScreenOutputFileName(m_prefix + "subScreenOutputFileName"),
  m_option_subScreenOutputAllowAll(m_prefix + "subScreenOutputAllowAll"),
  m_option_subScreenOutputAllow   (m_prefix + "subScreenOutputAllow"   ),
  m_option_verbosity              (m_prefix + "verbosity"              ),
  m_option_syncVerbosity          (m_prefix + "syncVerbosity"          ),
  m_option_seed                   (m_prefix + "seed"                   ),
  m_numSubEnvironments            (UQ_ENV_NUM_SUB_ENVIRONMENTS_ODV),
  m_subScreenOutputFileName       (UQ_ENV_SUB_SCREEN_OUTPUT_FILE_NAME_ODV),
  m_subScreenOutputAllowAll       (UQ_ENV_SUB_SCREEN_OUTPUT_ALLOW_ALL_ODV),
//m_subScreenOutputAllowSet       (),
  m_verbosity                     (UQ_ENV_VERBOSITY_ODV),
  m_syncVerbosity                 (UQ_ENV_SYNC_VERBOSITY_ODV),
  m_seed                          (UQ_ENV_SEED_ODV),
  m_numDebugParams                (UQ_ENV_NUM_DEBUG_PARAMS_ODV),
  m_debugParams                   (m_numDebugParams,0.),
  m_subComm                       (NULL),
  m_subRank                       (-1),
  m_subCommSize                   (1),
  m_selfComm                      (NULL),
  m_inter0Comm                    (NULL),
  m_inter0Rank                    (-1),
  m_inter0CommSize                (1),
  m_subScreenFile                 (NULL),
  m_rng                           (NULL)
{
}

uqBaseEnvironmentClass::uqBaseEnvironmentClass(
  int&        argc,
  char**      &argv,
  MPI_Comm    inputComm,
  const char* prefix)
  :
  m_argc                          (argc),
  m_argv                          (argv),
  m_prefix                        ((std::string)(prefix) + "env_"),
  m_worldRank                     (-1),
  m_fullRawComm                   (inputComm),
  m_fullComm                      (NULL),
  m_fullRank                      (-1),
  m_fullCommSize                  (0),
  m_argsWereProvided              (true),
  m_thereIsInputFile              (false),
  m_inputFileName                 (""),
  m_allOptionsDesc                (NULL),
  m_envOptionsDesc                (NULL),
  m_allOptionsMap                 (NULL),
  m_option_help                   (m_prefix + "help"                   ),
  m_option_numSubEnvironments     (m_prefix + "numSubEnvironments"     ),
  m_option_subScreenOutputFileName(m_prefix + "subScreenOutputFileName"),
  m_option_subScreenOutputAllowAll(m_prefix + "subScreenOutputAllowAll"),
  m_option_subScreenOutputAllow   (m_prefix + "subScreenOutputAllow"   ),
  m_option_verbosity              (m_prefix + "verbosity"              ),
  m_option_syncVerbosity          (m_prefix + "syncVerbosity"          ),
  m_option_seed                   (m_prefix + "seed"                   ),
  m_numSubEnvironments            (UQ_ENV_NUM_SUB_ENVIRONMENTS_ODV),
  m_subScreenOutputFileName       (UQ_ENV_SUB_SCREEN_OUTPUT_FILE_NAME_ODV),
  m_subScreenOutputAllowAll       (UQ_ENV_SUB_SCREEN_OUTPUT_ALLOW_ALL_ODV),
//m_subScreenOutputAllowSet       (),
  m_verbosity                     (UQ_ENV_VERBOSITY_ODV),
  m_syncVerbosity                 (UQ_ENV_SYNC_VERBOSITY_ODV),
  m_seed                          (UQ_ENV_SEED_ODV),
  m_numDebugParams                (UQ_ENV_NUM_DEBUG_PARAMS_ODV),
  m_debugParams                   (m_numDebugParams,0.),
  m_subComm                       (NULL),
  m_subRank                       (-1),
  m_subCommSize                   (0),
  m_selfComm                      (NULL),
  m_inter0Comm                    (NULL),
  m_inter0Rank                    (-1),
  m_inter0CommSize                (0),
  m_subScreenFile                 (NULL),
  m_rng                           (NULL)
{
}

uqBaseEnvironmentClass::uqBaseEnvironmentClass(
  const uqEnvOptionsStruct& options,
  MPI_Comm                  inputComm,
  const char*               prefix)
  :
  m_argc                          (0),
  m_argv                          (NULL),
  m_prefix                        ((std::string)(prefix) + "env_"),
  m_worldRank                     (-1),
  m_fullRawComm                   (inputComm),
  m_fullComm                      (NULL),
  m_fullRank                      (-1),
  m_fullCommSize                  (0),
  m_argsWereProvided              (false),
  m_thereIsInputFile              (false),
  m_inputFileName                 (""),
  m_allOptionsDesc                (NULL),
  m_envOptionsDesc                (NULL),
  m_allOptionsMap                 (NULL),
  m_option_help                   (m_prefix + "help"                   ),
  m_option_numSubEnvironments     (m_prefix + "numSubEnvironments"     ),
  m_option_subScreenOutputFileName(m_prefix + "subScreenOutputFileName"),
  m_option_subScreenOutputAllowAll(m_prefix + "subScreenOutputAllowAll"),
  m_option_subScreenOutputAllow   (m_prefix + "subScreenOutputAllow"   ),
  m_option_verbosity              (m_prefix + "verbosity"              ),
  m_option_syncVerbosity          (m_prefix + "syncVerbosity"          ),
  m_option_seed                   (m_prefix + "seed"                   ),
  m_numSubEnvironments            (UQ_ENV_NUM_SUB_ENVIRONMENTS_ODV),
  m_subScreenOutputFileName       (UQ_ENV_SUB_SCREEN_OUTPUT_FILE_NAME_ODV),
  m_subScreenOutputAllowAll       (UQ_ENV_SUB_SCREEN_OUTPUT_ALLOW_ALL_ODV),
//m_subScreenOutputAllowSet       (),
  m_verbosity                     (options.m_verbosity),
  m_syncVerbosity                 (UQ_ENV_SYNC_VERBOSITY_ODV),
  m_seed                          (options.m_seed),
  m_numDebugParams                (options.m_numDebugParams),
  m_debugParams                   (options.m_debugParams),
  m_subComm                       (NULL),
  m_subRank                       (-1),
  m_subCommSize                   (0),
  m_selfComm                      (NULL),
  m_inter0Comm                    (NULL),
  m_inter0Rank                    (-1),
  m_inter0CommSize                (0),
  m_subScreenFile                 (NULL),
  m_rng                           (NULL)
{
}

uqBaseEnvironmentClass::uqBaseEnvironmentClass(const uqBaseEnvironmentClass& obj)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      obj.rank(),
                      "uqBaseEnvironmentClass::constructor(), copy",
                      "code should not execute through here");
}

uqBaseEnvironmentClass::~uqBaseEnvironmentClass()
{
  //if (m_subScreenFile) {
  //  *m_subScreenFile << "Entering uqBaseEnvironmentClass::destructor()"
  //                   << std::endl;
  //}

  if (m_allOptionsMap) {
    delete m_allOptionsMap;
    delete m_envOptionsDesc;
    delete m_allOptionsDesc;
  }

  if (m_rng) gsl_rng_free(m_rng);

  struct timeval timevalNow;
  /*int iRC = 0;*/
  /*iRC = */gettimeofday(&timevalNow, NULL);

  if (m_subScreenFile) {
    *m_subScreenFile << "Ending run at " << ctime(&timevalNow.tv_sec)
                     << "Total run time = " << timevalNow.tv_sec - m_timevalBegin.tv_sec
                     << " seconds"
                     << std::endl;
  }

  if (m_fullRank == 0) {
    std::cout << "Ending run at " << ctime(&timevalNow.tv_sec)
              << "Total run time = " << timevalNow.tv_sec - m_timevalBegin.tv_sec
              << " seconds"
              << std::endl;
  }

  //if (m_subScreenFile) {
  //  *m_subScreenFile << "Leaving uqBaseEnvironmentClass::destructor()"
  //                   << std::endl;
  //}

  if (m_subScreenFile) delete m_subScreenFile;
  if (m_inter0Comm   ) delete m_inter0Comm;
  if (m_selfComm     ) delete m_selfComm;
  if (m_subComm      ) delete m_subComm;
  if (m_fullComm     ) delete m_fullComm;
}

uqBaseEnvironmentClass&
uqBaseEnvironmentClass::operator= (const uqBaseEnvironmentClass& rhs)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      rhs.rank(),
                      "uqBaseEnvironmentClass::operator=()",
                      "code should not execute through here");
  return *this;
}

int
uqBaseEnvironmentClass::rank() const
{
  return m_fullRank;
}

const Epetra_MpiComm&
uqBaseEnvironmentClass::fullComm() const
{
  return *m_fullComm;
}

int
uqBaseEnvironmentClass::subRank() const
{
  return m_subRank;
}

const Epetra_MpiComm&
uqBaseEnvironmentClass::subComm() const
{
  return *m_subComm;
}

const Epetra_MpiComm&
uqBaseEnvironmentClass::selfComm() const
{
  return *m_selfComm;
}

int
uqBaseEnvironmentClass::inter0Rank() const
{
  return m_inter0Rank;
}

const Epetra_MpiComm&
uqBaseEnvironmentClass::inter0Comm() const
{
  return *m_inter0Comm;
}

std::ofstream*
uqBaseEnvironmentClass::subScreenFile() const
{
  return m_subScreenFile;
}

unsigned int
uqBaseEnvironmentClass::numSubEnvironments() const
{
  return m_numSubEnvironments;
}

unsigned int
uqBaseEnvironmentClass::subId() const
{
  return m_subId;
}

const std::string&
uqBaseEnvironmentClass::subIdString() const
{
  return m_subIdString;
}

void
uqBaseEnvironmentClass::scanInputFileForMyOptions(const po::options_description& optionsDesc) const
{
#ifdef UQ_USES_COMMAND_LINE_OPTIONS
  // If you want to use command line options, the following line does *not* work outside 'main.C',
  // e.g., in the constructor of uqFullEnvironmentClass:
  // Line: po::store(po::parse_command_line(argc, argv, *m_allOptionsDesc), *m_allOptionsMap);
  //
  // Instead, put the following three lines *immediately after* instantianting the UQ environment
  // variable "uqFullEnvironmentClass* env":
  // Line 1: po::store(po::parse_command_line(argc, argv, env->allOptionsDesc()), env->allOptionsMap());
  // Line 2: po::notify(env->allOptionsMap());
  // Line 3: env->getMyOptionValues();
#endif

  m_allOptionsDesc->add(optionsDesc);
  //if (m_subScreenFile) {
  //  *m_subScreenFile << *m_allOptionsDesc
  //                   << std::endl;
  //}
  if (m_thereIsInputFile) {
    std::ifstream ifs(m_inputFileName.c_str());
    po::store(po::parse_config_file(ifs, *m_allOptionsDesc, true), *m_allOptionsMap);
    po::notify(*m_allOptionsMap);
    ifs.close();
  }

  return;
}

#ifdef UQ_USES_COMMAND_LINE_OPTIONS
const po::options_description&
uqBaseEnvironmentClass::allOptionsDesc() const
{
  return *m_allOptionsDesc;
}
#endif

po::variables_map&
uqBaseEnvironmentClass::allOptionsMap() const
{
  return *m_allOptionsMap;
}

unsigned int
uqBaseEnvironmentClass::verbosity() const
{
  return m_verbosity;
}

unsigned int
uqBaseEnvironmentClass::syncVerbosity() const
{
  return m_syncVerbosity;
}

const gsl_rng*
uqBaseEnvironmentClass::rng() const
{
  return m_rng;
}

bool
uqBaseEnvironmentClass::isThereInputFile() const
{
  return m_thereIsInputFile;
}

void
uqBaseEnvironmentClass::syncPrintDebugMsg(const char* msg, unsigned int msgVerbosity, unsigned int numUSecs, const Epetra_MpiComm& commObj) const
{
  commObj.Barrier();
  if (this->syncVerbosity() >= msgVerbosity) {
    for (int i = 0; i < commObj.NumProc(); ++i) {
      if (i == commObj.MyPID()) {
        std::cout << msg
                  << ": fullRank "       << this->rank()
                  << ", subEnvironment " << this->subId()
                  << ", subRank "        << this->subRank()
                  << ", inter0Rank "     << this->inter0Rank()
                  << std::endl;
      }
      commObj.Barrier();
    }
    if (this->rank() == 0) std::cout << "Sleeping " << numUSecs << " microseconds..."
                                     << std::endl;
    usleep(numUSecs);
  }
  commObj.Barrier();

  return;
}

void
uqBaseEnvironmentClass::openOutputFile(
  const std::string&            baseFileName,
  const std::string&            fileType,
  const std::set<unsigned int>& allowedSubEnvIds,
        bool                    writeOver,
        std::ofstream*&         ofsvar) const
{
  ofsvar = NULL;
  if ((baseFileName                         == UQ_ENV_FILENAME_FOR_NO_OUTPUT_FILE) ||
      (allowedSubEnvIds.find(this->subId()) == allowedSubEnvIds.end()            )) {
    if (this->subScreenFile()) {
      *this->subScreenFile() << "In openOutputFile()"
                             << ": no output file opened with base name '" << baseFileName
                             << "'"
                             << std::endl;
    }
  }
  else {
    if (this->subScreenFile()) {
      *this->subScreenFile() << "In openOutputFile()"
                             << ": opening output file with base name '" << baseFileName
                             << "'"
                             << std::endl;
    }

    if (this->subRank() == 0) {
      // Open file
      if (writeOver) {
        // Write over an eventual pre-existing file
        ofsvar = new std::ofstream((baseFileName+"_sub"+this->subIdString()+"."+fileType).c_str(), std::ofstream::out | std::ofstream::trunc);
      }
      else {
        // Write at the end of an eventual pre-existing file
        ofsvar = new std::ofstream((baseFileName+"_sub"+this->subIdString()+"."+fileType).c_str(), std::ofstream::out | std::ofstream::in | std::ofstream::ate);
        if ((ofsvar            == NULL ) ||
            (ofsvar->is_open() == false)) {
          delete ofsvar;
          ofsvar = new std::ofstream((baseFileName+"_sub"+this->subIdString()+"."+fileType).c_str(), std::ofstream::out | std::ofstream::trunc);
        }
      }
      if (ofsvar == NULL) {
        std::cerr << "In openOutputFile()"
                  << ": failed to open output file with base name '" << baseFileName
                  << "'"
                  << std::endl;
      }
      UQ_FATAL_TEST_MACRO((ofsvar && ofsvar->is_open()) == false,
                          this->rank(),
                          "openOutputFile()",
                          "failed to open output file");
    }
  }

  return;
}

void
uqBaseEnvironmentClass::openUnifiedOutputFile(
  const std::string&            baseFileName,
  const std::string&            fileType,
        bool                    writeOver,
        std::ofstream*&         ofsvar) const
{
  ofsvar = NULL;
  if (baseFileName == ".") {
    if (this->subScreenFile()) {
      *this->subScreenFile() << "In openUnifiedOutputFile()"
                             << ": no unified output file opened with base name '" << baseFileName
                             << "'"
                             << std::endl;
    }
  }
  else {
    if (this->subScreenFile()) {
      *this->subScreenFile() << "In openUnifiedOutputFile()"
                             << ": opening unified output file with base name '" << baseFileName
                             << "'"
                             << std::endl;
    }

    //if ((this->subRank   () == 0) &&
    //    (this->inter0Rank() == 0)) {
      // Open file
      if (writeOver) {
        // Write over an eventual pre-existing file
        ofsvar = new std::ofstream((baseFileName+"."+fileType).c_str(), std::ofstream::out | std::ofstream::trunc);
      }
      else {
        // Write at the end of an eventual pre-existing file
        ofsvar = new std::ofstream((baseFileName+"."+fileType).c_str(), std::ofstream::out | std::ofstream::in | std::ofstream::ate);
        if ((ofsvar            == NULL ) ||
            (ofsvar->is_open() == false)) {
          delete ofsvar;
          ofsvar = new std::ofstream((baseFileName+"."+fileType).c_str(), std::ofstream::out | std::ofstream::trunc);
        }
      }
      if (ofsvar == NULL) {
        std::cerr << "In openUnifiedOutputFile()"
                  << ": failed to open unified output file with base name '" << baseFileName
                  << "'"
                  << std::endl;
      }
      UQ_FATAL_TEST_MACRO((ofsvar && ofsvar->is_open()) == false,
                          this->rank(),
                          "openUnifiedOutputFile()",
                          "failed to open output file");
    //}
  }

  return;
}

//*****************************************************
// Empty Environment
//*****************************************************
uqEmptyEnvironmentClass::uqEmptyEnvironmentClass()
  :
  uqBaseEnvironmentClass(MPI_COMM_WORLD,"")
{
}

uqEmptyEnvironmentClass::~uqEmptyEnvironmentClass()
{
}

void
uqEmptyEnvironmentClass::print(std::ostream& os) const
{
  os.flush(); // just to avoid icpc warnings
  return;
}

//*****************************************************
// Full Environment
//*****************************************************
uqFullEnvironmentClass::uqFullEnvironmentClass(MPI_Comm inputComm, const char* prefix)
  :
  uqBaseEnvironmentClass(inputComm, prefix)
{
  commonConstructor(inputComm);
}

uqFullEnvironmentClass::uqFullEnvironmentClass(
  int&        argc,
  char**      &argv,
  MPI_Comm    inputComm,
  const char* prefix)
  :
  uqBaseEnvironmentClass(argc,argv,inputComm,prefix)
{
  commonConstructor(inputComm);
}

uqFullEnvironmentClass::uqFullEnvironmentClass(
  const uqEnvOptionsStruct& options,
  MPI_Comm                  inputComm,
  const char*               prefix)
  :
  uqBaseEnvironmentClass(options,inputComm,prefix)
{
  commonConstructor(inputComm);
}

uqFullEnvironmentClass::~uqFullEnvironmentClass()
{
}

void
uqFullEnvironmentClass::commonConstructor(MPI_Comm inputComm)
{
  //////////////////////////////////////////////////
  // Initialize "full" communicator
  //////////////////////////////////////////////////
  int mpiRC = MPI_Comm_rank(inputComm,&m_worldRank);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      UQ_UNAVAILABLE_RANK,
                      "uqFullEnvironmentClass::commonConstructor()",
                      "failed to get world rank()");

  m_fullComm = new Epetra_MpiComm(inputComm);
  m_fullRank     = m_fullComm->MyPID();
  m_fullCommSize = m_fullComm->NumProc();
  mpiRC = MPI_Comm_group(m_fullComm->Comm(), &m_fullGroup);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_fullRank,
                      "uqFullEnvironmentClass::commonConstructor()",
                      "failed MPI_Comm_group()");

  //////////////////////////////////////////////////
  // Display main initial messages
  // 'std::cout' is for: main trace messages + synchronized trace messages + error messages prior to 'exit()' or 'abort()'
  //////////////////////////////////////////////////
  /*int iRC = 0;*/
  /*iRC = */gettimeofday(&m_timevalBegin, NULL);

  if (m_fullRank == 0) {
    std::cout << "\n======================================================="
              << "\n QUESO library, version " << QUESO_LIBRARY_CURRENT_VERSION
              << ", released on "             << QUESO_LIBRARY_RELEASE_DATE
              << "\n======================================================="
              << "\n"
              << std::endl;
  }

  if (m_fullRank == 0) {
    std::cout << "Beginning run at " << ctime(&m_timevalBegin.tv_sec)
              << std::endl;
  }

  //////////////////////////////////////////////////
  // Read options
  //////////////////////////////////////////////////
  m_allOptionsMap  = new po::variables_map();
  m_allOptionsDesc = new po::options_description("Allowed options");
  m_envOptionsDesc = new po::options_description("Environment options");

  if (m_argsWereProvided) readEventualInputFile();

  defineMyOptions          (*m_envOptionsDesc);
  scanInputFileForMyOptions(*m_envOptionsDesc);
  getMyOptionValues        (*m_envOptionsDesc);

  // 'm_subScreenFile' is still not available at this moment. Use 'std::cout'
  if ((m_fullRank == 0) && (this->verbosity() >= 3)) {
    std::cout << "After getting option values, state of uqFullEnvironmentClass object is:"
              << "\n" << *this
              << std::endl;
  }

  //////////////////////////////////////////////////
  // Deal with multiple subEnvironments: create the sub communicators, one for each subEnvironment
  //////////////////////////////////////////////////
  unsigned int numRanksPerSubEnvironment = m_fullCommSize/m_numSubEnvironments;

  m_subId = m_fullRank/numRanksPerSubEnvironment;
  char tmpSubId[16];
  sprintf(tmpSubId,"%u",m_subId);
  m_subIdString = tmpSubId;

  std::vector<int> fullRanksOfMySubEnvironment(numRanksPerSubEnvironment,0);
  for (unsigned int i = 0; i < numRanksPerSubEnvironment; ++i) {
    fullRanksOfMySubEnvironment[i] = m_subId * numRanksPerSubEnvironment + i;
  }
  mpiRC = MPI_Group_incl(m_fullGroup, (int) numRanksPerSubEnvironment, &fullRanksOfMySubEnvironment[0], &m_subGroup);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_fullRank,
                      "uqFullEnvironmentClass::commonConstructor()",
                      "failed MPI_Group_incl() for a subEnvironment");
  mpiRC = MPI_Comm_create(m_fullComm->Comm(), m_subGroup, &m_subRawComm);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_fullRank,
                      "uqFullEnvironmentClass::commonConstructor()",
                      "failed MPI_Comm_group() for a subEnvironment");
  m_subComm = new Epetra_MpiComm(m_subRawComm);
  m_subRank     = m_subComm->MyPID();
  m_subCommSize = m_subComm->NumProc();

  //////////////////////////////////////////////////
  // Deal with multiple subEnvironments: create the self communicator
  //////////////////////////////////////////////////
  m_selfComm = new Epetra_MpiComm(MPI_COMM_SELF);

  //////////////////////////////////////////////////
  // Deal with multiple subEnvironments: create the inter0 communicator
  //////////////////////////////////////////////////
  std::vector<int> fullRanksOfInter0(m_numSubEnvironments,0);
  for (unsigned int i = 0; i < m_numSubEnvironments; ++i) {
    fullRanksOfInter0[i] = i * numRanksPerSubEnvironment;
  }
  mpiRC = MPI_Group_incl(m_fullGroup, (int) m_numSubEnvironments, &fullRanksOfInter0[0], &m_inter0Group);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_fullRank,
                      "uqFullEnvironmentClass::commonConstructor()",
                      "failed MPI_Group_incl() for inter0");
  mpiRC = MPI_Comm_create(m_fullComm->Comm(), m_inter0Group, &m_inter0RawComm);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_fullRank,
                      "uqFullEnvironmentClass::commonConstructor()",
                      "failed MPI_Comm_group() for inter0");
  if (m_fullRank%numRanksPerSubEnvironment == 0) {
    m_inter0Comm = new Epetra_MpiComm(m_inter0RawComm);
    m_inter0Rank     = m_inter0Comm->MyPID();
    m_inter0CommSize = m_inter0Comm->NumProc();
  }

  //////////////////////////////////////////////////
  // Open "screen" file
  //////////////////////////////////////////////////
  bool openFile = false;
  if ((m_subRank                               == 0                                 ) &&
      (m_subScreenOutputFileName               != UQ_ENV_FILENAME_FOR_NO_OUTPUT_FILE) &&
      (m_subScreenOutputAllowSet.find(m_subId) != m_subScreenOutputAllowSet.end()   )) {
    openFile = true;
  }

  if (openFile) {
    // Always write over an eventual pre-existing file
    m_subScreenFile = new std::ofstream((m_subScreenOutputFileName+"_sub"+m_subIdString+".txt").c_str(), std::ofstream::out | std::ofstream::trunc);
    UQ_FATAL_TEST_MACRO((m_subScreenFile && m_subScreenFile->is_open()) == false,
                        m_fullRank,
                        "uqEnvironment::constructor()",
                        "failed to open sub screen file");

    *m_subScreenFile << "\n================================="
                     << "\n QUESO library, version " << QUESO_LIBRARY_CURRENT_VERSION
                     << ", released on "             << QUESO_LIBRARY_RELEASE_DATE
                     << "\n================================="
                     << "\n"
                     << std::endl;

    *m_subScreenFile << "Beginning run at " << ctime(&m_timevalBegin.tv_sec)
                     << std::endl;
  }

  //////////////////////////////////////////////////
  // Debug message related to subEnvironments
  //////////////////////////////////////////////////
  if (this->verbosity() >= 2) {
    for (int i = 0; i < m_fullCommSize; ++i) {
      if (i == m_fullRank) {
        //std::cout << "In FullEnvironmentClass::commonConstructor()"
        std::cout << "MPI node of worldRank "             << m_worldRank
                  << " has fullRank "                     << m_fullRank
                  << ", belongs to subEnvironment of id " << m_subId
                  << ", and has subRank "                 << m_subRank
                  << std::endl;

        std::cout << "MPI node of worldRank " << m_worldRank
                  << " belongs to sub communicator with full ranks";
        for (unsigned int j = 0; j < fullRanksOfMySubEnvironment.size(); ++j) {
          std::cout << " " << fullRanksOfMySubEnvironment[j];
        }
	std::cout << "\n";

        if (m_inter0Comm) {
          std::cout << "MPI node of worldRank " << m_worldRank
                    << " also belongs to inter0 communicator with full ranks";
          for (unsigned int j = 0; j < fullRanksOfInter0.size(); ++j) {
            std::cout << " " << fullRanksOfInter0[j];
          }
          std::cout << ", and has inter0Rank " << m_inter0Rank;
        }
	std::cout << "\n";

	std::cout << std::endl;
      }
      m_fullComm->Barrier();
    }
    //if (this->rank() == 0) std::cout << "Sleeping 3 seconds..."
    //                                 << std::endl;
    //sleep(3);
  }

  //////////////////////////////////////////////////
  // Deal with seed
  //////////////////////////////////////////////////
  if (m_seed >= 0) {
    gsl_rng_default_seed = (unsigned long int) m_seed;
  }
  else if (m_seed == -1) {
    gsl_rng_default_seed = (unsigned long int) (1+m_fullRank);
  }
  else {
    struct timeval timevalNow;
    /*iRC = */gettimeofday(&timevalNow, NULL);
    gsl_rng_default_seed = (unsigned long int) timevalNow.tv_usec;
  }

  m_rng = gsl_rng_alloc(gsl_rng_ranlxd2);
  UQ_FATAL_TEST_MACRO((m_rng == NULL),
                      m_fullRank,
                      "uqFullEnvironmentClass::commonConstructor()",
                      "null m_rng");

  //gsl_rng_set(m_rng, gsl_rng_default_seed);

  if ((m_subScreenFile)/* && (this->verbosity() > 0)*/) {
    *m_subScreenFile << "In uqFullEnvironmentClass::commonConstructor():"
                     << "\n  m_seed = "                                              << m_seed
                     << "\n  internal seed = "                                       << gsl_rng_default_seed
                   //<< "\n  first generated sample from uniform distribution = "    << gsl_rng_uniform(m_rng)
                   //<< "\n  first generated sample from std normal distribution = " << gsl_ran_gaussian(m_rng,1.)
                     << std::endl;
  }

  //////////////////////////////////////////////////
  // Leave commonConstructor()
  //////////////////////////////////////////////////
  if ((m_subScreenFile) && (this->verbosity() >= 5)) {
    *m_subScreenFile << "Done with initializations at uqFullEnvironmentClass::commonConstructor()"
                     << std::endl;
  }

  return;
}

void
uqFullEnvironmentClass::readEventualInputFile()
{
  bool displayHelpMessageAndExit = false;
  bool invalidCmdLineParameters  = false;
  int maxArgIndex = m_argc - 1;
  int curArgIndex = 1;

  if (curArgIndex > maxArgIndex) {
    invalidCmdLineParameters = true;
  }
  else while (curArgIndex <= maxArgIndex) {
    if ((strcmp(m_argv[curArgIndex],"-h"    ) == 0) ||
        (strcmp(m_argv[curArgIndex],"--help") == 0)) {
      displayHelpMessageAndExit = true;
    }
    else if ((strcmp(m_argv[curArgIndex],"-i"     ) == 0) ||
             (strcmp(m_argv[curArgIndex],"--input") == 0)) {
      curArgIndex++;
      if (curArgIndex > maxArgIndex) {
        invalidCmdLineParameters = true;
        break;
      }
      m_inputFileName = std::string(m_argv[curArgIndex]);
      std::ifstream* ifs = new std::ifstream(m_inputFileName.c_str());
      if (ifs->is_open()) {
        m_thereIsInputFile = true;
        //if (m_subScreenFile) {
        //  int numLines = std::count(std::istreambuf_iterator<char>(*ifs),
        //                            std::istreambuf_iterator<char>(),'\n');
        //  *m_subScreenFile << "Input file has " << numLines
        //                   << " lines."
        //                   << std::endl;
        //}
        ifs->close();
        delete ifs;
      }
      else {
        invalidCmdLineParameters = true;
      }
    }
    //else {
    //  invalidCmdLineParameters = true;
    //}
    curArgIndex++;
  }

  if (invalidCmdLineParameters) {
    if (m_fullRank == 0) std::cout << "Invalid command line parameters!"
                                    << std::endl;
    displayHelpMessageAndExit = true;
  }

  if (displayHelpMessageAndExit) {
    if (m_fullRank == 0) std::cout << "\nThis is a help message of the QUESO library."
                                   << "\nAn application using the QUESO library shall be executed by typing"
                                   << "\n  '<eventual mpi commands and options> <uqApplication> -i <uqInputFile>'"
                                   << "\nin the command line."
                                   << "\n"
                                   << std::endl;
    /*int mpiRC = 0;*/
    /*mpiRC = */MPI_Abort(m_fullComm->Comm(),-999);
    exit(1);
  }

  return;
}

void
uqFullEnvironmentClass::defineMyOptions(po::options_description& options) const
{
  options.add_options()
    (m_option_help.c_str(),                                                                                                      "produce help message for uq environment"  )
    (m_option_numSubEnvironments.c_str(),      po::value<unsigned int>()->default_value(UQ_ENV_NUM_SUB_ENVIRONMENTS_ODV),        "number of subEnvironments"                )
    (m_option_subScreenOutputFileName.c_str(), po::value<std::string >()->default_value(UQ_ENV_SUB_SCREEN_OUTPUT_FILE_NAME_ODV), "output filename for subscreen writing"    )
    (m_option_subScreenOutputAllowAll.c_str(), po::value<bool        >()->default_value(UQ_ENV_SUB_SCREEN_OUTPUT_ALLOW_ALL_ODV), "Allow all subEnvs to write to output file")
    (m_option_subScreenOutputAllow.c_str(),    po::value<std::string >()->default_value(UQ_ENV_SUB_SCREEN_OUTPUT_ALLOW_ODV),     "subEnvs that will write to output file"   )
    (m_option_verbosity.c_str(),               po::value<unsigned int>()->default_value(UQ_ENV_VERBOSITY_ODV),                   "set verbosity"                            )
    (m_option_syncVerbosity.c_str(),           po::value<unsigned int>()->default_value(UQ_ENV_SYNC_VERBOSITY_ODV),              "set sync verbosity"                       )
    (m_option_seed.c_str(),                    po::value<int         >()->default_value(UQ_ENV_SEED_ODV),                        "set seed"                                 )
  //(m_option_numDebugParams.c_str(),          po::value<unsigned int>()->default_value(UQ_ENV_NUM_DEBUG_PARAMS_ODV),            "set number of debug parameters"           )
  ;

  return;
}

void
uqFullEnvironmentClass::getMyOptionValues(po::options_description& optionsDesc)
{
  if (m_allOptionsMap->count(m_option_help.c_str())) {
    // 'm_subScreenFile' is still not available at this moment. Use 'std::cout'
    if (m_fullRank == 0) std::cout << optionsDesc
                                   << std::endl;
  }

  if (m_allOptionsMap->count(m_option_numSubEnvironments.c_str())) {
    m_numSubEnvironments = (*m_allOptionsMap)[m_option_numSubEnvironments.c_str()].as<unsigned int>();
  }
  if ((m_fullCommSize%m_numSubEnvironments) != 0) {
    std::cerr << "In uqBaseEnvironmentClass::getMyOptionValues()"
              << ": m_fullCommSize = "       << m_fullCommSize
              << ", m_numSubEnvironments = " << m_numSubEnvironments
              << std::endl;
  }
  UQ_FATAL_TEST_MACRO((m_fullCommSize%m_numSubEnvironments) != 0,
                      this->rank(),
                      "uqBaseEnvironmentClass::getMyOptionValues()",
                      "total number of processors in environment must be multiple of the specified number of subEnvironments");

  if (m_allOptionsMap->count(m_option_subScreenOutputFileName.c_str())) {
    m_subScreenOutputFileName = (*m_allOptionsMap)[m_option_subScreenOutputFileName.c_str()].as<std::string>();
  }

  if (m_allOptionsMap->count(m_option_subScreenOutputAllowAll.c_str())) {
    m_subScreenOutputAllowAll = (*m_allOptionsMap)[m_option_subScreenOutputAllowAll.c_str()].as<bool>();
  }

  if (m_subScreenOutputAllowAll) {
    m_subScreenOutputAllowSet.clear();
    m_subScreenOutputAllowSet.insert((unsigned int) m_subId);
  }
  else if (m_allOptionsMap->count(m_option_subScreenOutputAllow.c_str())) {
    m_subScreenOutputAllowSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = (*m_allOptionsMap)[m_option_subScreenOutputAllow.c_str()].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpAllow);
    //if (m_subScreenFile) {
    //  *m_subScreenFile << "In uqFullEnvironmentClass::getMyOptionValues(): allow = ";
    //  for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
    //    *m_subScreenFile << " " << tmpAllow[i];
    //  }
    //  *m_subScreenFile << std::endl;
    //}

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_subScreenOutputAllowSet.insert((unsigned int) tmpAllow[i]);
      }
    }
  }

  if (m_allOptionsMap->count(m_option_verbosity.c_str())) {
    m_verbosity = (*m_allOptionsMap)[m_option_verbosity.c_str()].as<unsigned int>();
  }

  if (m_allOptionsMap->count(m_option_syncVerbosity.c_str())) {
    m_syncVerbosity = (*m_allOptionsMap)[m_option_syncVerbosity.c_str()].as<unsigned int>();
  }

  if (m_allOptionsMap->count(m_option_seed.c_str())) {
    m_seed = (*m_allOptionsMap)[m_option_seed.c_str()].as<int>();
  }

  //if (m_allOptionsMap->count(m_option_numDebugParams.c_str())) {
  //  m_seed = (*m_allOptionsMap)[m_option_numDebugParams.c_str()].as<unsigned int>();
  //}

  return;
}

void
uqFullEnvironmentClass::print(std::ostream& os) const
{
  os <<         m_option_numSubEnvironments      << " = " << m_numSubEnvironments
     << "\n" << m_option_subScreenOutputFileName << " = " << m_subScreenOutputFileName
     << "\n" << m_option_subScreenOutputAllowAll << " = ";
  os << "\n" << m_option_subScreenOutputAllow << " = ";
  for (std::set<unsigned int>::iterator setIt = m_subScreenOutputAllowSet.begin(); setIt != m_subScreenOutputAllowSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_verbosity               << " = " << m_verbosity
     << "\n" << m_option_syncVerbosity           << " = " << m_syncVerbosity
     << "\n" << m_option_seed                    << " = " << m_seed
   //<< "\n" << m_option_numDebugParams          << " = " << m_numDebugParams
     << std::endl;
  return;
}

std::ostream&
operator<<(std::ostream& os, const uqBaseEnvironmentClass& obj)
{
  obj.print(os);

  return os;
}
