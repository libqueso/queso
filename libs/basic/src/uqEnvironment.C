/* libs/basic/src/uqEnvironment.C
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

#include <uqEnvironment.h>
#include <uqDefines.h>
#include <sys/time.h>
#include <gsl/gsl_randist.h>

// Version "0.1"  on "Aug/11/2008"
// Version "0.11" on "Aug/15/2008"
// Version "0.2"  on "Oct/01/2008"
// Version "0.21" on "Oct/08/2008"
#define QUESO_TOOLKIT_CURRENT_VERSION "0.21"
#define QUESO_TOOLKIT_RELEASE_DATE    "Oct/08/2008"

uqEnvOptionsStruct::uqEnvOptionsStruct(
  unsigned int verbosity,
  int          seed)
  :
  m_verbosity     (verbosity),
  m_seed          (seed),
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
uqBaseEnvironmentClass::uqBaseEnvironmentClass()
  :
  m_argc            (0),
  m_argv            (NULL),
  m_comm            (NULL),
  m_rank            (0),
  m_commSize        (1),
  m_argsWereProvided(false),
  m_thereIsInputFile(false),
  m_inputFileName   (""),
  m_allOptionsDesc  (NULL),
  m_envOptionsDesc  (NULL),
  m_allOptionsMap   (NULL),
  m_verbosity       (UQ_ENV_VERBOSITY_ODV),
  m_seed            (UQ_ENV_SEED_ODV),
  m_numDebugParams  (UQ_ENV_NUM_DEBUG_PARAMS_ODV),
  m_debugParams     (m_numDebugParams,0.),
  m_rng             (NULL)
{
}

uqBaseEnvironmentClass::uqBaseEnvironmentClass(
  int&   argc,
  char** &argv)
  :
  m_argc            (argc),
  m_argv            (argv),
  m_comm            (NULL),
  m_rank            (0),
  m_commSize        (1),
  m_argsWereProvided(true),
  m_thereIsInputFile(false),
  m_inputFileName   (""),
  m_allOptionsDesc  (NULL),
  m_envOptionsDesc  (NULL),
  m_allOptionsMap   (NULL),
  m_verbosity       (UQ_ENV_VERBOSITY_ODV),
  m_seed            (UQ_ENV_SEED_ODV),
  m_numDebugParams  (UQ_ENV_NUM_DEBUG_PARAMS_ODV),
  m_debugParams     (m_numDebugParams,0.),
  m_rng             (NULL)
{
}

uqBaseEnvironmentClass::uqBaseEnvironmentClass(
  const uqEnvOptionsStruct& options)
  :
  m_argc            (0),
  m_argv            (NULL),
  m_comm            (NULL),
  m_rank            (0),
  m_commSize        (1),
  m_argsWereProvided(false),
  m_thereIsInputFile(false),
  m_inputFileName   (""),
  m_allOptionsDesc  (NULL),
  m_envOptionsDesc  (NULL),
  m_allOptionsMap   (NULL),
  m_verbosity       (options.m_verbosity),
  m_seed            (options.m_seed),
  m_numDebugParams  (options.m_numDebugParams),
  m_debugParams     (options.m_debugParams),
  m_rng             (NULL)
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
  //if (m_rank == 0) std::cout << "Entering uqBaseEnvironmentClass::destructor()"
  //                           << std::endl;

  if (m_allOptionsMap) {
    delete m_allOptionsMap;
    delete m_envOptionsDesc;
    delete m_allOptionsDesc;
  }

  if (m_rng) gsl_rng_free(m_rng);

  int iRC;
  struct timeval timevalNow;
  iRC = gettimeofday(&timevalNow, NULL);
  if (m_rank == 0) {
    std::cout << "Ending run at " << ctime(&timevalNow.tv_sec)
              << "Total run time = " << timevalNow.tv_sec - m_timevalBegin.tv_sec
              << " seconds"
              << std::endl;
  }

  //if (m_rank == 0) std::cout << "Leaving uqBaseEnvironmentClass::destructor()"
  //                           << std::endl;

  if (m_comm) delete m_comm;
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
  return m_rank;
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

void
uqBaseEnvironmentClass::barrier() const
{
  if (m_commSize > 1) {
    m_comm->Barrier();
  }
  return;
}

const Epetra_MpiComm&
uqBaseEnvironmentClass::comm() const
{
  return *m_comm;
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
  uqBaseEnvironmentClass()
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
uqFullEnvironmentClass::uqFullEnvironmentClass()
  :
  uqBaseEnvironmentClass()
{
  commonConstructor();
}

uqFullEnvironmentClass::uqFullEnvironmentClass(
  int&   argc,
  char** &argv)
  :
  uqBaseEnvironmentClass(argc,argv)
{
  commonConstructor();
}

uqFullEnvironmentClass::uqFullEnvironmentClass(
  const uqEnvOptionsStruct& options)
  :
  uqBaseEnvironmentClass(options)
{
  commonConstructor();
}

uqFullEnvironmentClass::~uqFullEnvironmentClass()
{
}

void
uqFullEnvironmentClass::commonConstructor()
{
  m_comm     = new Epetra_MpiComm(MPI_COMM_WORLD);
  m_rank     = m_comm->MyPID();
  m_commSize = m_comm->NumProc();

  int iRC;
  iRC = gettimeofday(&m_timevalBegin, NULL);

  if ((this->verbosity() >= 5) && (this->rank() == 0)) {
    std::cout << "Entering uqFullEnvironmentClass::commonConstructor()"
              << std::endl;
  }

  if (m_rank == 0) {
    std::cout << "\n================================="
              << "\n QUESO toolkit, version " << QUESO_TOOLKIT_CURRENT_VERSION
              << ", released on "             << QUESO_TOOLKIT_RELEASE_DATE
              << "\n================================="
              << "\n"
              << std::endl;
  }

  if (m_rank == 0) {
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
                      m_rank,
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
        //if (m_rank == 0) {
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
    if (m_rank == 0) std::cout << "Invalid command line parameters!"
                               << std::endl;
    displayHelpMessageAndExit = true;
  }

  if (displayHelpMessageAndExit) {
    if (m_rank == 0) std::cout << "\nThis is a help message of the UQ library."
                               << "\nAn UQ application using the MCMC tool ot the QUESO toolkit shall be executed by typing"
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
    ("uqEnv_help", "produce help message for uq environment")
    ("uqEnv_verbosity",      po::value<unsigned int>()->default_value(UQ_ENV_VERBOSITY_ODV),       "set verbosity"                 )
    ("uqEnv_seed",           po::value<int         >()->default_value(UQ_ENV_SEED_ODV),            "set seed"                      )
    //("uqEnv_numDebugParams", po::value<unsigned int>()->default_value(UQ_ENV_NUM_DEBUG_PARAMS_ODV),"set number of debug parameters")
  ;

  return;
}

void
uqFullEnvironmentClass::getMyOptionValues(po::options_description& optionsDesc)
{
  if (m_allOptionsMap->count("uqEnv_help")) {
    std::cout << optionsDesc
              << std::endl;
  }

  if (m_allOptionsMap->count("uqEnv_verbosity")) {
    m_verbosity = (*m_allOptionsMap)["uqEnv_verbosity"].as<unsigned int>();
  }

  if (m_allOptionsMap->count("uqEnv_seed")) {
    m_seed = (*m_allOptionsMap)["uqEnv_seed"].as<int>();
  }

  //if (m_allOptionsMap->count("uqEnv_numDebugParams")) {
  //  m_seed = (*m_allOptionsMap)["uqEnv_numDebugParams"].as<unsigned int>();
  //}

  if (m_verbosity >= 1) {
    if (m_rank == 0) std::cout << "After getting option values, state of uqFullEnvironmentClass object is:"
                               << "\n" << *this
                               << std::endl;
  }

  return;
}

void
uqFullEnvironmentClass::print(std::ostream& os) const
{
  os << "m_verbosity = "        << m_verbosity
     << "\nm_seed = "           << m_seed
    //<< "\nm_numDebugParams = " << m_numDebugParams
     << std::endl;
  return;
}

std::ostream&
operator<<(std::ostream& os, const uqBaseEnvironmentClass& obj)
{
  obj.print(os);

  return os;
}
