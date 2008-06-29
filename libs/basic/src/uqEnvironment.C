#include <uqEnvironment.h>
#include <uqDefines.h>

#define PECOS_TOOLKIT_CURRENT_VERSION "0.1"

uqEnvOptionsStruct::uqEnvOptionsStruct(
  unsigned int verbosity)
  :
  m_verbosity(verbosity)
{
}

uqEnvOptionsStruct::~uqEnvOptionsStruct()
{
}

uqEnvironmentClass::uqEnvironmentClass()
  :
  m_argc            (0),
  m_argv            (NULL),
#ifdef __UQ_USES_TRILINOS__
  m_comm            (NULL),
#endif
  m_rank            (0),
  m_commSize        (1),
  m_argsWereProvided(false),
  m_thereIsInputFile(false),
  m_inputFileName   (""),
  m_allOptionsDesc  (NULL),
  m_envOptionsDesc  (NULL),
  m_allOptionsMap   (NULL),
  m_verbosity       (UQ_ENV_VERBOSITY_DEFAULT_VALUE),
  m_rng             (NULL)
{
  commonConstructor();
}

uqEnvironmentClass::uqEnvironmentClass(
  int&   argc,
  char** &argv)
  :
  m_argc            (argc),
  m_argv            (argv),
#ifdef __UQ_USES_TRILINOS__
  m_comm            (NULL),
#endif
  m_rank            (0),
  m_commSize        (1),
  m_argsWereProvided(true),
  m_thereIsInputFile(false),
  m_inputFileName   (""),
  m_allOptionsDesc  (NULL),
  m_envOptionsDesc  (NULL),
  m_allOptionsMap   (NULL),
  m_verbosity       (UQ_ENV_VERBOSITY_DEFAULT_VALUE),
  m_rng             (NULL)
{
  //////////////////////////////////////////////////
  // Initialize MPI
  //////////////////////////////////////////////////
  commonConstructor();
}

uqEnvironmentClass::uqEnvironmentClass(
  const uqEnvOptionsStruct& options)
  :
  m_argc            (0),
  m_argv            (NULL),
#ifdef __UQ_USES_TRILINOS__
  m_comm            (NULL),
#endif
  m_rank            (0),
  m_commSize        (1),
  m_argsWereProvided(false),
  m_thereIsInputFile(false),
  m_inputFileName   (""),
  m_allOptionsDesc  (NULL),
  m_envOptionsDesc  (NULL),
  m_allOptionsMap   (NULL),
  m_verbosity       (options.m_verbosity),
  m_rng             (NULL)
{
  commonConstructor();
}

void
uqEnvironmentClass::commonConstructor()
{
#ifdef __UQ_USES_TRILINOS__
  m_comm     = new Epetra_MpiComm(MPI_COMM_WORLD);
  m_rank     = m_comm->MyPID();
  m_commSize = m_comm->Size();
#endif

  //if (m_rank == 0) std::cout << "Entering uqEnvironmentClass::commonConstructor()"
  //                           << std::endl;

  if (m_rank == 0) std::cout << "\n================================="
                             << "\n PECOS toolkit, version " << PECOS_TOOLKIT_CURRENT_VERSION
                             << "\n================================="
                             << "\n"
                             << std::endl;

  m_rng = gsl_rng_alloc(gsl_rng_ranlxd2);
  UQ_FATAL_TEST_MACRO((m_rng == NULL),
                      m_rank,
                      "uqEnvironmentClass::commonConstructor()",
                      "null m_rng");

  m_allOptionsMap  = new po::variables_map();
  m_allOptionsDesc = new po::options_description("Allowed options");
  m_envOptionsDesc = new po::options_description("Environment options");

  if (m_argsWereProvided) readEventualInputFile();

  defineMyOptions          (*m_envOptionsDesc);
  scanInputFileForMyOptions(*m_envOptionsDesc);
  getMyOptionValues        (*m_envOptionsDesc);

  //if (m_rank == 0) std::cout << "Leaving uqEnvironmentClass::commonConstructor()"
  //                           << std::endl;

  return;
}

void
uqEnvironmentClass::readEventualInputFile()
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
                               << "\nAn UQ application using the PECOS toolkit shall be executed by typing"
                               << "\n  '<eventual mpi commands and options> <uqApplication> -i <uqInputFile>'"
                               << "\nin the command line."
                               << "\n"
                               << std::endl;
    exit(1);
  }

  return;
}

uqEnvironmentClass::uqEnvironmentClass(const uqEnvironmentClass& obj)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      this->rank(),
                      "uqEnvironmentClass::constructor(), copy",
                      "code should not execute through here");

  m_argc           = obj.m_argc;
  m_argv           = obj.m_argv;
#ifdef __UQ_USES_TRILINOS__
  m_comm           = NULL;
#endif
  m_rank           = obj.m_rank;
  m_commSize       = obj.m_commSize;
  m_allOptionsDesc = NULL;
  m_envOptionsDesc = NULL;
  m_allOptionsMap  = obj.m_allOptionsMap;
  m_verbosity      = obj.m_verbosity;
  m_rng            = obj.m_rng;
}

uqEnvironmentClass::~uqEnvironmentClass()
{
  //if (m_rank == 0) std::cout << "Entering uqEnvironmentClass::destructor()"
  //                           << std::endl;

  if (m_allOptionsMap) {
    delete m_allOptionsMap;
    delete m_envOptionsDesc;
    delete m_allOptionsDesc;
  }

  if (m_rng) gsl_rng_free(m_rng);

  //if (m_rank == 0) std::cout << "Leaving uqEnvironmentClass::destructor()"
  //                           << std::endl;

  //////////////////////////////////////////////////
  // Finalize MPI
  //////////////////////////////////////////////////
#ifdef __UQ_USES_TRILINOS__
  if (m_comm)    delete m_comm;
#endif
}

uqEnvironmentClass&
uqEnvironmentClass::operator= (const uqEnvironmentClass& rhs)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      this->rank(),
                      "uqEnvironmentClass::operator=()",
                      "code should not execute through here");

  m_argc      = rhs.m_argc;
  m_argv      = rhs.m_argv;
#ifdef __UQ_USES_TRILINOS__
  m_comm      = rhs.m_comm;
#endif
  m_rank      = rhs.m_rank;
  m_commSize  = rhs.m_commSize;
  m_verbosity = rhs.m_verbosity;
  m_rng       = rhs.m_rng;

  return *this;
}

int
uqEnvironmentClass::rank() const
{
  return m_rank;
}

void
uqEnvironmentClass::barrier() const
{
#ifdef __UQ_USES_TRILINOS__
  if (m_commSize > 1) {
    m_comm->Barrier();
  }
#endif
  return;
}

#ifdef UQ_USES_COMMAND_LINE_OPTIONS
const po::options_description&
uqEnvironmentClass::allOptionsDesc() const
{
  return *m_allOptionsDesc;
}
#endif

po::variables_map&
uqEnvironmentClass::allOptionsMap() const
{
  return *m_allOptionsMap;
}

void
uqEnvironmentClass::defineMyOptions(po::options_description& options) const
{
  options.add_options()
    ("uqEnv_help", "produce help message for uq environment")
    ("uqEnv_verbosity", po::value<unsigned int>()->default_value(UQ_ENV_VERBOSITY_DEFAULT_VALUE), "set verbosity")
  ;

  return;
}

void
uqEnvironmentClass::scanInputFileForMyOptions(const po::options_description& optionsDesc) const
{
#ifdef UQ_USES_COMMAND_LINE_OPTIONS
  // If you want to use command line options, the following line does *not* work outside 'main.C',
  // e.g., in the constructor of uqEnvironmentClass:
  // Line: po::store(po::parse_command_line(argc, argv, *m_allOptionsDesc), *m_allOptionsMap);
  //
  // Instead, put the following three lines *immediately after* instantianting the UQ environment
  // variable "uqEnvironmentClass* env":
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
uqEnvironmentClass::getMyOptionValues(po::options_description& optionsDesc)
{
  if (m_allOptionsMap->count("uqEnv_help")) {
    std::cout << optionsDesc
              << std::endl;
  }

  if (m_allOptionsMap->count("uqEnv_verbosity")) {
    m_verbosity = (*m_allOptionsMap)["uqEnv_verbosity"].as<unsigned int>();
  }

  if (m_verbosity >= 1) {
    if (m_rank == 0) std::cout << "After getting option values, state of uqEnvironmentClass object is:"
                               << "\n" << *this
                               << std::endl;
  }

  return;
}

unsigned int
uqEnvironmentClass::verbosity() const
{
  return m_verbosity;
}

gsl_rng*
uqEnvironmentClass::rng() const
{
  return m_rng;
}

bool
uqEnvironmentClass::isThereInputFile() const
{
  return m_thereIsInputFile;
}

void
uqEnvironmentClass::print(std::ostream& os) const
{
  os << "m_verbosity = " << m_verbosity
     << std::endl;
  return;
}

std::ostream&
operator<<(std::ostream& os, const uqEnvironmentClass& obj)
{
  obj.print(os);

  return os;
}
