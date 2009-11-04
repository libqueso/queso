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
#include <uqEnvironmentOptions.h>
#include <uqMiscellaneous.h>
#include <sys/time.h>
#include <gsl/gsl_randist.h>
#include <hpct.h>

// Version "0.1"    on "Aug/11/2008"
// Version "0.11"   on "Aug/15/2008"
// Version "0.2"    on "Oct/01/2008"
// Version "0.21"   on "Oct/08/2008"
// Version "0.3.0"  on "Feb/13/2009"
// Version "0.3.1"  on "Feb/19/2009"
// Version "0.4.0"  on "Jul/22/2009"
// Version "0.4.1"  on "Sep/08/2009"
// Version "0.40.2" on "Sep/10/2009"
// Version "0.41.0" on "Oct/30/2009"
// Version "0.42.0" on "MMM/DD/20YY"
#define QUESO_LIBRARY_CURRENT_VERSION "0.42.0"
#define QUESO_LIBRARY_RELEASE_DATE    "MMM/DD/20YY"

//*****************************************************
// Base class
//*****************************************************
uqBaseEnvironmentClass::uqBaseEnvironmentClass(
  MPI_Comm    inputComm,
  const char* optionsInputFileName)
  :
  m_worldRank           (-1),
  m_fullRawComm         (inputComm),
  m_fullComm            (NULL),
  m_fullRank            (-1),
  m_fullCommSize        (1),
  m_optionsInputFileName(optionsInputFileName),
  m_allOptionsDesc      (NULL),
  m_allOptionsMap       (NULL),
  m_subComm             (NULL),
  m_subRank             (-1),
  m_subCommSize         (1),
  m_selfComm            (NULL),
  m_inter0Comm          (NULL),
  m_inter0Rank          (-1),
  m_inter0CommSize      (1),
  m_subDisplayFile      (NULL),
  m_rng                 (NULL),
  m_options             (NULL)
{
}

uqBaseEnvironmentClass::uqBaseEnvironmentClass(const uqBaseEnvironmentClass& obj)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      obj.fullRank(),
                      "uqBaseEnvironmentClass::constructor(), copy",
                      "code should not execute through here");
}

uqBaseEnvironmentClass::~uqBaseEnvironmentClass()
{
  //if (m_subDisplayFile) {
  //  *m_subDisplayFile << "Entering uqBaseEnvironmentClass::destructor()"
  //                          << std::endl;
  //}


  struct timeval timevalNow;
  /*int iRC = 0;*/
  /*iRC = */gettimeofday(&timevalNow, NULL);

  if( this->displayVerbosity() > 0 )
    {

      if (m_subDisplayFile) {
	*m_subDisplayFile << "Ending run at "    << ctime(&timevalNow.tv_sec)
			  << "Total run time = " << timevalNow.tv_sec - m_timevalBegin.tv_sec
			  << " seconds"
			  << std::endl;
      }

      if (m_fullRank == 0) {
	std::cout << "Ending run at "    << ctime(&timevalNow.tv_sec)
		  << "Total run time = " << timevalNow.tv_sec - m_timevalBegin.tv_sec
		  << " seconds"
		  << std::endl;
      }

    }

  if (m_options) delete m_options;

  if (m_allOptionsMap) {
    delete m_allOptionsMap;
    delete m_allOptionsDesc;
  }

  if (m_rng) gsl_rng_free(m_rng);

  //if (m_subDisplayFile) {
  //  *m_subDisplayFile << "Leaving uqBaseEnvironmentClass::destructor()"
  //                          << std::endl;
  //}

  if (m_subDisplayFile) delete m_subDisplayFile;
  if (m_inter0Comm    ) delete m_inter0Comm;
  if (m_selfComm      ) delete m_selfComm;
  if (m_subComm       ) delete m_subComm;
  if (m_fullComm      ) delete m_fullComm;
}

uqBaseEnvironmentClass&
uqBaseEnvironmentClass::operator= (const uqBaseEnvironmentClass& rhs)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      rhs.fullRank(),
                      "uqBaseEnvironmentClass::operator=()",
                      "code should not execute through here");
  return *this;
}

int
uqBaseEnvironmentClass::worldRank() const
{
  return m_worldRank;
}

int
uqBaseEnvironmentClass::fullRank() const
{
  return m_fullRank;
}

const Epetra_MpiComm&
uqBaseEnvironmentClass::fullComm() const
{
  UQ_FATAL_TEST_MACRO(m_fullComm == NULL,
                      m_fullRank,
                      "uqBaseEnvironmentClass::fullComm()",
                      "m_fullComm variable is NULL");
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
  UQ_FATAL_TEST_MACRO(m_subComm == NULL,
                      m_fullRank,
                      "uqBaseEnvironmentClass::subComm()",
                      "m_subComm variable is NULL");
  return *m_subComm;
}

const Epetra_MpiComm&
uqBaseEnvironmentClass::selfComm() const
{
  UQ_FATAL_TEST_MACRO(m_selfComm == NULL,
                      m_fullRank,
                      "uqBaseEnvironmentClass::selfComm()",
                      "m_selfComm variable is NULL");
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
  UQ_FATAL_TEST_MACRO(m_inter0Comm == NULL,
                      m_fullRank,
                      "uqBaseEnvironmentClass::inter0Comm()",
                      "m_inter0Comm variable is NULL");
  return *m_inter0Comm;
}

std::ofstream*
uqBaseEnvironmentClass::subDisplayFile() const
{
  return m_subDisplayFile;
}

unsigned int
uqBaseEnvironmentClass::numSubEnvironments() const
{
  UQ_FATAL_TEST_MACRO(m_options == NULL,
                      m_fullRank,
                      "uqBaseEnvironmentClass::numSubEnvironments()",
                      "m_options variable is NULL");
  return m_options->m_numSubEnvironments;
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
  //if (m_subDisplayFile) {
  //  *m_subDisplayFile << *m_allOptionsDesc
  //                          << std::endl;
  //}
  std::ifstream ifs(m_optionsInputFileName.c_str());
  po::store(po::parse_config_file(ifs, *m_allOptionsDesc, true), *m_allOptionsMap);
  po::notify(*m_allOptionsMap);
  ifs.close();

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
  UQ_FATAL_TEST_MACRO(m_allOptionsMap == NULL,
                      m_fullRank,
                      "uqBaseEnvironmentClass::allOptionsMap()",
                      "m_allOptionsMap variable is NULL");
  return *m_allOptionsMap;
}

unsigned int
uqBaseEnvironmentClass::displayVerbosity() const
{
  UQ_FATAL_TEST_MACRO(m_options == NULL,
                      m_fullRank,
                      "uqBaseEnvironmentClass::displayVerbosity()",
                      "m_options variable is NULL");
  return m_options->m_displayVerbosity;
}

unsigned int
uqBaseEnvironmentClass::syncVerbosity() const
{
  UQ_FATAL_TEST_MACRO(m_options == NULL,
                      m_fullRank,
                      "uqBaseEnvironmentClass::displayVerbosity()",
                      "m_options variable is NULL");
  return m_options->m_syncVerbosity;
}

const gsl_rng*
uqBaseEnvironmentClass::rng() const
{
  return m_rng;
}

void
uqBaseEnvironmentClass::syncPrintDebugMsg(const char* msg, unsigned int msgVerbosity, unsigned int numUSecs, const Epetra_MpiComm& commObj) const
{
  commObj.Barrier();
  if (this->syncVerbosity() >= msgVerbosity) {
    for (int i = 0; i < commObj.NumProc(); ++i) {
      if (i == commObj.MyPID()) {
        std::cout << msg
                  << ": fullRank "       << this->fullRank()
                  << ", subEnvironment " << this->subId()
                  << ", subRank "        << this->subRank()
                  << ", inter0Rank "     << this->inter0Rank()
                  << std::endl;
      }
      commObj.Barrier();
    }
    if (this->fullRank() == 0) std::cout << "Sleeping " << numUSecs << " microseconds..."
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
    if ((m_subDisplayFile) && (this->displayVerbosity() >= 10)) {
      *this->subDisplayFile() << "In openOutputFile()"
                              << ": no output file opened with base name '" << baseFileName
                              << "'"
                              << std::endl;
    }
  }
  else {
    if ((m_subDisplayFile) && (this->displayVerbosity() >= 10)) {
      *this->subDisplayFile() << "In openOutputFile()"
                              << ": opening output file with base name '" << baseFileName
                              << "'"
                              << std::endl;
    }

    if (this->subRank() == 0) {
#if 0
      std::cout << "In openOutputFile()"
                << ": opening output file with base name '" << baseFileName
                << "'"
                << ", writeOver = " << writeOver
                << std::endl;
#endif

      // Verify parent directory exists (for cases when a user
      // specifies a relative path for the desired output file).
#if 0
      std::cout << "checking " << baseFileName+"_sub"+this->subIdString()+"."+fileType << std::endl;
#endif
      int irtrn = hpct_check_file_path((baseFileName+"_sub"+this->subIdString()+"."+fileType).c_str());
      UQ_FATAL_TEST_MACRO(irtrn < 0,m_fullRank,"openOutputFile()","unable to verify output path");

      // Open file
      if (writeOver) {
        // Write over an eventual pre-existing file
        ofsvar = new std::ofstream((baseFileName+"_sub"+this->subIdString()+"."+fileType).c_str(), 
				   std::ofstream::out | std::ofstream::trunc);
      }
      else {
        // Write at the end of an eventual pre-existing file
#if 0
        ofsvar = new std::ofstream((baseFileName+"_sub"+this->subIdString()+"."+fileType).c_str(), 
				   std::ofstream::in);
	std::cout << "ofsvar(1) = " << ofsvar << std::endl;
        if (ofsvar) std::cout << "ofsvar(1)->is_open() = " << ofsvar->is_open() << std::endl;
        if (ofsvar) delete ofsvar;
#endif
        ofsvar = new std::ofstream((baseFileName+"_sub"+this->subIdString()+"."+fileType).c_str(), 
				   std::ofstream::out | std::ofstream::in | std::ofstream::app);
#if 0
	std::cout << "ofsvar(2) = " << ofsvar << std::endl;
        if (ofsvar) std::cout << "ofsvar(2)->is_open() = " << ofsvar->is_open() << std::endl;
#endif
        if ((ofsvar            == NULL ) ||
            (ofsvar->is_open() == false)) {
#if 0
	  std::cout << "Retrying 1..." << std::endl;
#endif
          delete ofsvar;
          ofsvar = new std::ofstream((baseFileName+"_sub"+this->subIdString()+"."+fileType).c_str(), 
				     std::ofstream::out | std::ofstream::trunc);
        }
      }
      if (ofsvar == NULL) {
        std::cerr << "In openOutputFile()"
                  << ": failed to open output file with base name '" << baseFileName
                  << "'"
                  << std::endl;
      }
      UQ_FATAL_TEST_MACRO((ofsvar && ofsvar->is_open()) == false,
                          this->fullRank(),
                          "openOutputFile()",
                          "failed to open output file");
    }
  }

  return;
}

void
uqBaseEnvironmentClass::openUnifiedOutputFile(
  const std::string&    baseFileName,
  const std::string&    fileType,
        bool            writeOver,
        std::ofstream*& ofsvar) const
{
  ofsvar = NULL;
  if (baseFileName == ".") {
    if ((m_subDisplayFile) && (this->displayVerbosity() >= 10)) {
      *this->subDisplayFile() << "In openUnifiedOutputFile()"
                              << ": no unified output file opened with base name '" << baseFileName
                              << "'"
                              << std::endl;
    }
  }
  else {
    if ((m_subDisplayFile) && (this->displayVerbosity() >= 10)) {
      *this->subDisplayFile() << "In openUnifiedOutputFile()"
                              << ": opening unified output file with base name '" << baseFileName
                              << "'"
                              << std::endl;
    }

    //if ((this->subRank   () == 0) &&
    //    (this->inter0Rank() == 0)) {
#if 0
      std::cout << "In openUnifiedOutputFile()"
                << ": opening output file with base name '" << baseFileName
                << "'"
                << ", writeOver = " << writeOver
                << std::endl;
#endif
      // Open file
      if (writeOver) {
        // Write over an eventual pre-existing file
        ofsvar = new std::ofstream((baseFileName+"."+fileType).c_str(), std::ofstream::out | std::ofstream::trunc);
      }
      else {
        // Write at the end of an eventual pre-existing file
        ofsvar = new std::ofstream((baseFileName+"."+fileType).c_str(), std::ofstream::out | std::ofstream::in | std::ofstream::app);
        if ((ofsvar            == NULL ) ||
            (ofsvar->is_open() == false)) {
#if 0
          std::cout << "Retrying 2..." << std::endl;
#endif
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
                          this->fullRank(),
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
uqFullEnvironmentClass::uqFullEnvironmentClass(
  MPI_Comm    inputComm,
  const char* optionsInputFileName,
  const char* prefix)
  :
  uqBaseEnvironmentClass(inputComm,optionsInputFileName)
{
  //////////////////////////////////////////////////
  // Initialize "full" communicator
  //////////////////////////////////////////////////
  int mpiRC = MPI_Comm_rank(MPI_COMM_WORLD,&m_worldRank);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      UQ_UNAVAILABLE_RANK,
                      "uqFullEnvironmentClass::commonConstructor()",
                      "failed to get world fullRank()");

  m_fullComm = new Epetra_MpiComm(m_fullRawComm);
  m_fullRank     = m_fullComm->MyPID();
  m_fullCommSize = m_fullComm->NumProc();
  mpiRC = MPI_Comm_group(m_fullComm->Comm(), &m_fullGroup);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_fullRank,
                      "uqFullEnvironmentClass::commonConstructor()",
                      "failed MPI_Comm_group()");

  //////////////////////////////////////////////////
  // Read options
  //////////////////////////////////////////////////
  m_allOptionsMap  = new po::variables_map();
  m_allOptionsDesc = new po::options_description("Allowed options");
  m_options = new uqEnvironmentOptionsClass(*this,prefix);

  readOptionsInputFile();

  m_options->scanOptionsValues();

  // Only display these messages if the user wants them
  // NOTE: This got moved below the Read Options section
  // because we need the options to be read to know what
  // the verbosity level is.
  if( this->displayVerbosity() > 0 )
    {

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
      
    }

  //////////////////////////////////////////////////
  // Deal with multiple subEnvironments: create the sub communicators, one for each subEnvironment
  //////////////////////////////////////////////////
  unsigned int numRanksPerSubEnvironment = m_fullCommSize/m_options->m_numSubEnvironments;

  m_subId = m_fullRank/numRanksPerSubEnvironment;
  char tmpSubId[16];
  sprintf(tmpSubId,"%u",m_subId);
  m_subIdString = tmpSubId;

  if (m_options->m_subDisplayAllowAll) {
    m_options->m_subDisplayAllowedSet.insert((unsigned int) m_subId);
  }

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
  std::vector<int> fullRanksOfInter0(m_options->m_numSubEnvironments,0);
  for (unsigned int i = 0; i < m_options->m_numSubEnvironments; ++i) {
    fullRanksOfInter0[i] = i * numRanksPerSubEnvironment;
  }
  mpiRC = MPI_Group_incl(m_fullGroup, (int) m_options->m_numSubEnvironments, &fullRanksOfInter0[0], &m_inter0Group);
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
  if ((m_subRank                                       == 0                                      ) &&
      (m_options->m_subDisplayFileName                 != UQ_ENV_FILENAME_FOR_NO_OUTPUT_FILE     ) &&
      (m_options->m_subDisplayAllowedSet.find(m_subId) != m_options->m_subDisplayAllowedSet.end())) {
    openFile = true;
  }

  if (openFile) {

    int irtrn = hpct_check_file_path((m_options->m_subDisplayFileName+"_sub"+m_subIdString+".txt").c_str());

    UQ_FATAL_TEST_MACRO(irtrn < 0,m_fullRank,"uqEnvironment::constructor()","unable to verify output path");
			

    // Always write over an eventual pre-existing file
    m_subDisplayFile = new std::ofstream((m_options->m_subDisplayFileName+"_sub"+m_subIdString+".txt").c_str(), std::ofstream::out | std::ofstream::trunc);
    UQ_FATAL_TEST_MACRO((m_subDisplayFile && m_subDisplayFile->is_open()) == false,
                        m_fullRank,
                        "uqEnvironment::constructor()",
                        "failed to open sub screen file");

    *m_subDisplayFile << "\n================================="
                      << "\n QUESO library, version " << QUESO_LIBRARY_CURRENT_VERSION
                      << ", released on "             << QUESO_LIBRARY_RELEASE_DATE
                      << "\n================================="
                      << "\n"
                      << std::endl;

    *m_subDisplayFile << "Beginning run at " << ctime(&m_timevalBegin.tv_sec)
                      << std::endl;
  }

  //////////////////////////////////////////////////
  // Debug message related to subEnvironments
  //////////////////////////////////////////////////
  if (this->displayVerbosity() >= 2) {
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
    //if (this->fullRank() == 0) std::cout << "Sleeping 3 seconds..."
    //                                 << std::endl;
    //sleep(3);
  }

  //////////////////////////////////////////////////
  // Deal with seed
  //////////////////////////////////////////////////
  if (m_options->m_seed >= 0) {
    gsl_rng_default_seed = (unsigned long int) m_options->m_seed;
  }
  else if (m_options->m_seed == -1) {
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

  if ((m_subDisplayFile)/* && (this->displayVerbosity() > 0)*/) {
    *m_subDisplayFile << "In uqFullEnvironmentClass::commonConstructor():"
                      << "\n  m_seed = "                                              << m_options->m_seed
                      << "\n  internal seed = "                                       << gsl_rng_default_seed
                    //<< "\n  first generated sample from uniform distribution = "    << gsl_rng_uniform(m_rng)
                    //<< "\n  first generated sample from std normal distribution = " << gsl_ran_gaussian(m_rng,1.)
                      << std::endl;
  }

  //////////////////////////////////////////////////
  // Leave commonConstructor()
  //////////////////////////////////////////////////
  if ((m_subDisplayFile) && (this->displayVerbosity() >= 5)) {
    *m_subDisplayFile << "Done with initializations at uqFullEnvironmentClass::commonConstructor()"
                      << std::endl;
  }

  return;
}

uqFullEnvironmentClass::~uqFullEnvironmentClass()
{
}

void
uqFullEnvironmentClass::readOptionsInputFile()
{
  std::ifstream* ifs = new std::ifstream(m_optionsInputFileName.c_str());
  if (ifs->is_open()) {
    //ifs->close();
    delete ifs;
  }
  else {
    if (m_fullRank == 0) std::cout << "An invalid input file has been passed to the 'environment' class constructor!"
                                   << std::endl;

    if (m_fullRank == 0) std::cout << "\nThis is a help message of the QUESO library."
                                   << "\nAn application using the QUESO library shall be executed by typing"
                                   << "\n  '<eventual mpi commands and options> <uqApplication> <uqInputFile>'"
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
uqFullEnvironmentClass::print(std::ostream& os) const
{
  os.flush(); // just to avoid icpc warnings
  return;
}

