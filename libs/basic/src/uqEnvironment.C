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
#include <sys/time.h>
#include <gsl/gsl_randist.h>

// Version "0.1"   on "Aug/11/2008"
// Version "0.11"  on "Aug/15/2008"
// Version "0.2"   on "Oct/01/2008"
// Version "0.21"  on "Oct/08/2008"
// Version "0.3.0" on "Feb/13/2009"
// Version "0.3.1" on "Feb/19/2009"
// Version "0.4.0" on "MMM/DD/2009"
#define QUESO_TOOLKIT_CURRENT_VERSION "0.4.0"
#define QUESO_TOOLKIT_RELEASE_DATE    "MMM/DD/2009"

uqEnvOptionsStruct::uqEnvOptionsStruct(
  unsigned int verbosity,
  int          seed)
  :
  m_numProcSubsets(UQ_ENV_NUM_PROC_SUBSETS_ODV),
  m_verbosity     (verbosity),
  m_seed          (seed),
  m_runName       (UQ_ENV_RUN_NAME_ODV),
  m_numDebugParams(0),
  m_debugParams   (0,0.)
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
  m_argc                 (0),
  m_argv                 (NULL),
  m_prefix               ((std::string)(prefix) + "env_"),
  m_fullComm             (NULL),
  m_fullRank             (0),
  m_fullCommSize         (1),
  m_argsWereProvided     (false),
  m_thereIsInputFile     (false),
  m_inputFileName        (""),
  m_allOptionsDesc       (NULL),
  m_envOptionsDesc       (NULL),
  m_allOptionsMap        (NULL),
  m_option_help          (m_prefix + "help"          ),
  m_option_numProcSubsets(m_prefix + "numProcSubsets"),
  m_option_verbosity     (m_prefix + "verbosity"     ),
  m_option_seed          (m_prefix + "seed"          ),
  m_option_runName       (m_prefix + "runName"       ),
  m_numProcSubsets       (UQ_ENV_NUM_PROC_SUBSETS_ODV),
  m_verbosity            (UQ_ENV_VERBOSITY_ODV),
  m_seed                 (UQ_ENV_SEED_ODV),
  m_runName              (UQ_ENV_RUN_NAME_ODV),
  m_numDebugParams       (UQ_ENV_NUM_DEBUG_PARAMS_ODV),
  m_debugParams          (m_numDebugParams,0.),
  m_subComm              (NULL),
  m_subRank              (0),
  m_subCommSize          (1),
  m_selfComm             (NULL),
  m_rng                  (NULL)
{
}

uqBaseEnvironmentClass::uqBaseEnvironmentClass(
  int&        argc,
  char**      &argv,
  MPI_Comm    inputComm,
  const char* prefix)
  :
  m_argc                 (argc),
  m_argv                 (argv),
  m_prefix               ((std::string)(prefix) + "env_"),
  m_fullComm             (NULL),
  m_fullRank             (0),
  m_fullCommSize         (1),
  m_argsWereProvided     (true),
  m_thereIsInputFile     (false),
  m_inputFileName        (""),
  m_allOptionsDesc       (NULL),
  m_envOptionsDesc       (NULL),
  m_allOptionsMap        (NULL),
  m_option_help          (m_prefix + "help"          ),
  m_option_numProcSubsets(m_prefix + "numProcSubsets"),
  m_option_verbosity     (m_prefix + "verbosity"     ),
  m_option_seed          (m_prefix + "seed"          ),
  m_option_runName       (m_prefix + "runName"       ),
  m_numProcSubsets       (UQ_ENV_NUM_PROC_SUBSETS_ODV),
  m_verbosity            (UQ_ENV_VERBOSITY_ODV),
  m_seed                 (UQ_ENV_SEED_ODV),
  m_runName              (UQ_ENV_RUN_NAME_ODV),
  m_numDebugParams       (UQ_ENV_NUM_DEBUG_PARAMS_ODV),
  m_debugParams          (m_numDebugParams,0.),
  m_subComm              (NULL),
  m_subRank              (0),
  m_subCommSize          (1),
  m_selfComm             (NULL),
  m_rng                  (NULL)
{
}

uqBaseEnvironmentClass::uqBaseEnvironmentClass(
  const uqEnvOptionsStruct& options,
  MPI_Comm                  inputComm,
  const char*               prefix)
  :
  m_argc                 (0),
  m_argv                 (NULL),
  m_prefix               ((std::string)(prefix) + "env_"),
  m_fullComm             (NULL),
  m_fullRank             (0),
  m_fullCommSize         (1),
  m_argsWereProvided     (false),
  m_thereIsInputFile     (false),
  m_inputFileName        (""),
  m_allOptionsDesc       (NULL),
  m_envOptionsDesc       (NULL),
  m_allOptionsMap        (NULL),
  m_option_help          (m_prefix + "help"          ),
  m_option_numProcSubsets(m_prefix + "numProcSubsets"),
  m_option_verbosity     (m_prefix + "verbosity"     ),
  m_option_seed          (m_prefix + "seed"          ),
  m_option_runName       (m_prefix + "runName"       ),
  m_numProcSubsets       (UQ_ENV_NUM_PROC_SUBSETS_ODV),
  m_verbosity            (options.m_verbosity),
  m_seed                 (options.m_seed),
  m_runName              (UQ_ENV_RUN_NAME_ODV),
  m_numDebugParams       (options.m_numDebugParams),
  m_debugParams          (options.m_debugParams),
  m_subComm              (NULL),
  m_subRank              (0),
  m_subCommSize          (1),
  m_selfComm             (NULL),
  m_rng                  (NULL)
{
}

uqBaseEnvironmentClass::uqBaseEnvironmentClass(const uqBaseEnvironmentClass& obj)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      this->rank(),
                      "uqBaseEnvironmentClass::constructor(), copy",
                      "code should not execute through here");
}

uqBaseEnvironmentClass::~uqBaseEnvironmentClass()
{
  //if (m_fullRank == 0) std::cout << "Entering uqBaseEnvironmentClass::destructor()"
  //                                << std::endl;

  if (m_allOptionsMap) {
    delete m_allOptionsMap;
    delete m_envOptionsDesc;
    delete m_allOptionsDesc;
  }

  if (m_rng) gsl_rng_free(m_rng);

  int iRC;
  struct timeval timevalNow;
  iRC = gettimeofday(&timevalNow, NULL);
  if (m_fullRank == 0) {
    std::cout << "Ending run at " << ctime(&timevalNow.tv_sec)
              << "Total run time = " << timevalNow.tv_sec - m_timevalBegin.tv_sec
              << " seconds"
              << std::endl;
  }

  //if (m_fullRank == 0) std::cout << "Leaving uqBaseEnvironmentClass::destructor()"
  //                                << std::endl;

  if (m_selfComm ) delete m_selfComm;
  if (m_subComm  ) delete m_subComm;
  if (m_fullComm ) delete m_fullComm;
}

uqBaseEnvironmentClass&
uqBaseEnvironmentClass::operator= (const uqBaseEnvironmentClass& rhs)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      this->rank(),
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

unsigned int
uqBaseEnvironmentClass::numProcSubsets() const
{
  return m_numProcSubsets;
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
  //std::cout << *m_allOptionsDesc
  //          << std::endl;
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

const std::string&
uqBaseEnvironmentClass::runName() const
{
  return m_runName;
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
  m_fullComm = new Epetra_MpiComm(inputComm);
  m_fullRank     = m_fullComm->MyPID();
  m_fullCommSize = m_fullComm->NumProc();
  int mpiRC = MPI_Comm_group(m_fullComm->Comm(), &m_fullGroup);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_fullRank,
                      "uqFullEnvironmentClass::commonConstructor()",
                      "failed MPI_Comm_group()");

  int iRC;
  iRC = gettimeofday(&m_timevalBegin, NULL);

  if ((this->verbosity() >= 5) && (this->rank() == 0)) {
    std::cout << "Entering uqFullEnvironmentClass::commonConstructor()"
              << std::endl;
  }

  if (m_fullRank == 0) {
    std::cout << "\n================================="
              << "\n QUESO toolkit, version " << QUESO_TOOLKIT_CURRENT_VERSION
              << ", released on "             << QUESO_TOOLKIT_RELEASE_DATE
              << "\n================================="
              << "\n"
              << std::endl;
  }

  if (m_fullRank == 0) {
    std::cout << "Beginning run at " << ctime(&m_timevalBegin.tv_sec)
              << std::endl;
  }

  m_allOptionsMap  = new po::variables_map();
  m_allOptionsDesc = new po::options_description("Allowed options");
  m_envOptionsDesc = new po::options_description("Environment options");

  if (m_argsWereProvided) readEventualInputFile();

  defineMyOptions          (*m_envOptionsDesc);
  scanInputFileForMyOptions(*m_envOptionsDesc);
  getMyOptionValues        (*m_envOptionsDesc);

  if (m_verbosity >= 1) {
    if (m_fullRank == 0) std::cout << "After getting option values, state of uqFullEnvironmentClass object is:"
                                    << "\n" << *this
                                    << std::endl;
  }

  // Deal with multiple subGroups
  unsigned int numRanksPerProcSubset = m_fullCommSize/m_numProcSubsets;

  m_subId = m_fullRank/numRanksPerProcSubset;
  char tmpSubId[16];
  sprintf(tmpSubId,"%d",m_subId);
  m_subIdString = tmpSubId;

  std::vector<int> fullRanksOfMyProcSubset(numRanksPerProcSubset,0);
  for (unsigned int i = 0; i < numRanksPerProcSubset; ++i) {
    fullRanksOfMyProcSubset[i] = m_subId * numRanksPerProcSubset + i;
  }
  mpiRC = MPI_Group_incl(m_fullGroup, (int) numRanksPerProcSubset, &fullRanksOfMyProcSubset[0], &m_subGroup);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_fullRank,
                      "uqFullEnvironmentClass::commonConstructor()",
                      "failed MPI_Group_incl()");
  mpiRC = MPI_Comm_create(m_fullComm->Comm(), m_subGroup, &m_subRawComm);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_fullRank,
                      "uqFullEnvironmentClass::commonConstructor()",
                      "failed MPI_Comm_group()");
  m_subComm = new Epetra_MpiComm(m_subRawComm);
  m_subRank     = m_subComm->MyPID();
  m_subCommSize = m_subComm->NumProc();

  m_selfComm = new Epetra_MpiComm(MPI_COMM_SELF);

#ifdef UQ_DEBUG_PARALLEL_RUNS_IN_DETAIL
  if (this->verbosity() >= 0) {
    for (int i = 0; i < m_fullCommSize; ++i) {
      if (i == m_fullRank) {
        std::cout << "In FullEnvironmentClass::commonConstructor()"
                  << ": fullRank " << m_fullRank
                  << " belongs to processor subset of id " << m_subId
                  << ", belongs to communicator with full ranks";
        for (unsigned int j = 0; j < numRanksPerProcSubset; ++j) {
          std::cout << " " << fullRanksOfMyProcSubset[j];
        }
        std::cout << ", and has subRank " << m_subRank
                  << std::endl;
      }
      m_fullComm->Barrier();
    }
    if (this->rank() == 0) std::cout << "Sleeping 3 seconds..."
                                     << std::endl;
    sleep(3);
  }
#endif

  // Deal with seed
  if (m_seed >= 0) {
    gsl_rng_default_seed = (unsigned long int) m_seed;
  }
  else {
    int iRC;
    struct timeval timevalNow;
    iRC = gettimeofday(&timevalNow, NULL);
    gsl_rng_default_seed = (unsigned long int) timevalNow.tv_usec;
  }

  m_rng = gsl_rng_alloc(gsl_rng_ranlxd2);
  UQ_FATAL_TEST_MACRO((m_rng == NULL),
                      m_fullRank,
                      "uqFullEnvironmentClass::commonConstructor()",
                      "null m_rng");

  if ((this->verbosity() >= 5) && (this->rank() == 0)) {
    std::cout << "In uqFullEnvironmentClass::commonConstructor():"
              << "\n  m_seed = "                                              << m_seed
              << "\n  internal seed = "                                       << gsl_rng_default_seed
              << "\n  first generated sample from uniform distribution = "    << gsl_rng_uniform(m_rng)
              << "\n  first generated sample from std normal distribution = " << gsl_ran_gaussian(m_rng,1.)
              << std::endl;
  }

  if ((this->verbosity() >= 5) && (this->rank() == 0)) {
    std::cout << "Leaving uqFullEnvironmentClass::commonConstructor()"
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
        //if (m_fullRank == 0) {
        //  int numLines = std::count(std::istreambuf_iterator<char>(*ifs),
        //                            std::istreambuf_iterator<char>(),'\n');
        //  std::cout << "Input file has " << numLines
        //             << " lines."
        //            << std::endl;
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
                                    << "\nAn application using the QUESO toolkit shall be executed by typing"
                                    << "\n  '<eventual mpi commands and options> <uqApplication> -i <uqInputFile>'"
                                    << "\nin the command line."
                                    << "\n"
                                    << std::endl;
    exit(1);
  }

  return;
}

void
uqFullEnvironmentClass::defineMyOptions(po::options_description& options) const
{
  options.add_options()
    (m_option_help.c_str(),                                                                                  "produce help message for uq environment")
    (m_option_numProcSubsets.c_str(), po::value<unsigned int>()->default_value(UQ_ENV_NUM_PROC_SUBSETS_ODV), "number of processor subsets"            )
    (m_option_verbosity.c_str(),      po::value<unsigned int>()->default_value(UQ_ENV_VERBOSITY_ODV),        "set verbosity"                          )
    (m_option_seed.c_str(),           po::value<int         >()->default_value(UQ_ENV_SEED_ODV),             "set seed"                               )
    (m_option_runName.c_str(),        po::value<std::string >()->default_value(UQ_ENV_RUN_NAME_ODV),         "set run name"                           )
  //(m_option_numDebugParams.c_str(), po::value<unsigned int>()->default_value(UQ_ENV_NUM_DEBUG_PARAMS_ODV), "set number of debug parameters"         )
  ;

  return;
}

void
uqFullEnvironmentClass::getMyOptionValues(po::options_description& optionsDesc)
{
  if (m_allOptionsMap->count(m_option_help.c_str())) {
    std::cout << optionsDesc
              << std::endl;
  }

  if (m_allOptionsMap->count(m_option_numProcSubsets.c_str())) {
    m_numProcSubsets = (*m_allOptionsMap)[m_option_numProcSubsets.c_str()].as<unsigned int>();
  }
  UQ_FATAL_TEST_MACRO((m_fullCommSize%m_numProcSubsets) != 0,
                      this->rank(),
                      "uqBaseEnvironmentClass::getMyOptionValues()",
                      "total number of processors in environment must be multiple of the specified number of processor subsets");

  if (m_allOptionsMap->count(m_option_verbosity.c_str())) {
    m_verbosity = (*m_allOptionsMap)[m_option_verbosity.c_str()].as<unsigned int>();
  }

  if (m_allOptionsMap->count(m_option_seed.c_str())) {
    m_seed = (*m_allOptionsMap)[m_option_seed.c_str()].as<int>();
  }

  if (m_allOptionsMap->count(m_option_runName.c_str())) {
    m_runName = (*m_allOptionsMap)[m_option_runName.c_str()].as<std::string >();
  }

  //if (m_allOptionsMap->count(m_option_numDebugParams.c_str())) {
  //  m_seed = (*m_allOptionsMap)[m_option_numDebugParams.c_str()].as<unsigned int>();
  //}

  return;
}

void
uqFullEnvironmentClass::print(std::ostream& os) const
{
  os <<         m_option_numProcSubsets << " = " << m_numProcSubsets
     << "\n" << m_option_verbosity      << " = " << m_verbosity
     << "\n" << m_option_seed           << " = " << m_seed
     << "\n" << m_option_runName        << " = " << m_runName
   //<< "\n" << m_option_numDebugParams << " = " << m_numDebugParams
     << std::endl;
  return;
}

std::ostream&
operator<<(std::ostream& os, const uqBaseEnvironmentClass& obj)
{
  obj.print(os);

  return os;
}
