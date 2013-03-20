//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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

#include <queso.h>
#include <uqEnvironment.h>
#include <uqEnvironmentOptions.h>
#include <uqRngGsl.h>
#include <uqRngBoost.h>
#include <uqMiscellaneous.h>
#include <sys/time.h>
#ifdef HAVE_GRVY
#include <grvy.h>
#endif

//-----------------------------
// Library versioning routines
//-----------------------------

namespace QUESO {

  void QUESO_version_print(std::ostream &os)
  {
    {
      os << "------------------------------------------------------------------------------------------" ;
      os << "--------------------" << std::endl;
      os << "QUESO Library: Version = " << QUESO_LIB_VERSION;
      os << " (" << QUESO_get_numeric_version() << ")" << std::endl << std::endl;
      
      os << QUESO_LIB_RELEASE << std::endl << std::endl;
      
      os << "Build Date   = " << QUESO_BUILD_DATE     << std::endl;
      os << "Build Host   = " << QUESO_BUILD_HOST     << std::endl;
      os << "Build User   = " << QUESO_BUILD_USER     << std::endl;
      os << "Build Arch   = " << QUESO_BUILD_ARCH     << std::endl;
      os << "Build Rev    = " << QUESO_BUILD_VERSION  << std::endl     << std::endl;
      
      os << "C++ Config   = " << QUESO_CXX << " " << QUESO_CXXFLAGS    << std::endl;
      os << std::endl;
      os << "Trilinos DIR = " << QUESO_TRILINOS_DIR << std::endl;
      os << "GSL Libs     = " << QUESO_GSL_DIR  << std::endl;
      os << "GRVY DIR     = " << QUESO_GRVY_DIR << std::endl;
      os << "GLPK DIR     = " << QUESO_GLPK_DIR << std::endl;
      os << "HDF5 DIR     = " << QUESO_HDF5_DIR << std::endl;
      os << "------------------------------------------------------------------------------------------" ;
      os << "--------------------" << std::endl;
    }
    
    return;
  }
  
  int QUESO_get_numeric_version()
  {
    // Note: return format follows the versioning convention xx.yy.z where
    //
    // xx = major version number 
    // yy = minor version number 
    // zz = micro version number
    //
    // For example:
    // v.    0.23  -> 002300 = 2300
    // v.   0.23.1 -> 002301 = 2301
    // v.  10.23.2 -> 102302
    
    int major_version = 0;
    int minor_version = 0;
    int micro_version = 0;
    
#ifdef QUESO_MAJOR_VERSION
    major_version = QUESO_MAJOR_VERSION;
#endif
    
#ifdef QUESO_MINOR_VERSION
    minor_version = QUESO_MINOR_VERSION;
#endif
    
#ifdef QUESO_MICRO_VERSION
    micro_version = QUESO_MICRO_VERSION;
#endif
    
    return(major_version*10000 + minor_version*100 + micro_version);
  }

}

uqFilePtrSetStruct::uqFilePtrSetStruct()
  :
  ofsVar(NULL),
  ifsVar(NULL)
{
}

uqFilePtrSetStruct::~uqFilePtrSetStruct()
{
}

//*****************************************************
// Base class
//*****************************************************
uqBaseEnvironmentClass::uqBaseEnvironmentClass(
  const char*                    passedOptionsInputFileName,
  const uqEnvOptionsValuesClass* alternativeOptionsValues)
  :
  m_fullEnvIsReady             (false),
  m_worldRank                  (-1),
  m_fullComm                   (NULL),
  m_fullRank                   (-1),
  m_fullCommSize               (1),
  m_optionsInputFileName       (""),
  m_optionsInputFileAccessState(true),
  m_allOptionsDesc             (NULL),
  m_allOptionsMap              (NULL),
  m_subComm                    (NULL),
  m_subRank                    (-1),
  m_subCommSize                (1),
  m_selfComm                   (NULL),
  m_inter0Comm                 (NULL),
  m_inter0Rank                 (-1),
  m_inter0CommSize             (1),
  m_subDisplayFile             (NULL),
  m_rngObject                  (NULL),
  m_exceptionalCircunstance    (false),
  m_alternativeOptionsValues   (),
  m_optionsObj                 (NULL)
{
  if (passedOptionsInputFileName) m_optionsInputFileName     = passedOptionsInputFileName;
  if (alternativeOptionsValues  ) m_alternativeOptionsValues = *alternativeOptionsValues;
}

uqBaseEnvironmentClass::uqBaseEnvironmentClass(const uqBaseEnvironmentClass& obj)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      obj.worldRank(),
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

  if (m_optionsObj) delete m_optionsObj;
  if (m_allOptionsMap) {
    delete m_allOptionsMap;
    delete m_allOptionsDesc;
  }

  if (m_rngObject) delete m_rngObject;

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
                      rhs.worldRank(),
                      "uqBaseEnvironmentClass::operator=()",
                      "code should not execute through here");
  return *this;
}

bool
uqBaseEnvironmentClass::fullEnvIsReady() const
{
  return m_fullEnvIsReady;
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

const uqMpiCommClass&
uqBaseEnvironmentClass::fullComm() const
{
  UQ_FATAL_TEST_MACRO(m_fullComm == NULL,
                      m_worldRank,
                      "uqBaseEnvironmentClass::fullComm()",
                      "m_fullComm variable is NULL");
  return *m_fullComm;
}

uqRawType_MPI_Group
uqBaseEnvironmentClass::subGroup() const
{
  return m_subGroup;
}

int
uqBaseEnvironmentClass::subRank() const
{
  return m_subRank;
}

const uqMpiCommClass&
uqBaseEnvironmentClass::subComm() const
{
  UQ_FATAL_TEST_MACRO(m_subComm == NULL,
                      m_worldRank,
                      "uqBaseEnvironmentClass::subComm()",
                      "m_subComm variable is NULL");
  return *m_subComm;
}

const uqMpiCommClass&
uqBaseEnvironmentClass::selfComm() const
{
  UQ_FATAL_TEST_MACRO(m_selfComm == NULL,
                      m_worldRank,
                      "uqBaseEnvironmentClass::selfComm()",
                      "m_selfComm variable is NULL");
  return *m_selfComm;
}

int
uqBaseEnvironmentClass::inter0Rank() const
{
  return m_inter0Rank;
}

const uqMpiCommClass&
uqBaseEnvironmentClass::inter0Comm() const
{
  UQ_FATAL_TEST_MACRO(m_inter0Comm == NULL,
                      m_worldRank,
                      "uqBaseEnvironmentClass::inter0Comm()",
                      "m_inter0Comm variable is NULL");
  return *m_inter0Comm;
}

std::ofstream*
uqBaseEnvironmentClass::subDisplayFile() const
{
  return m_subDisplayFile;
}

std::string
uqBaseEnvironmentClass::subDisplayFileName() const
{
  if (m_optionsObj == NULL) return ".";

  return m_optionsObj->m_ov.m_subDisplayFileName;
}

unsigned int
uqBaseEnvironmentClass::numSubEnvironments() const
{
  UQ_FATAL_TEST_MACRO(m_optionsObj == NULL,
                      m_worldRank,
                      "uqBaseEnvironmentClass::numSubEnvironments()",
                      "m_optionsObj variable is NULL");
  return m_optionsObj->m_ov.m_numSubEnvironments;
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

#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "Entering uqBaseEnvClass::scanInputFileForMyOptions()" << std::endl;
#endif

  UQ_FATAL_TEST_MACRO(m_allOptionsDesc == NULL,
                      m_worldRank,
                      "uqBaseEnvironmentClass::scanInputFileForMyOptions()",
                      "m_allOptionsDesc variable is NULL");
  m_allOptionsDesc->add(optionsDesc);
  //if (m_subDisplayFile) {
  //  *m_subDisplayFile << *m_allOptionsDesc
  //                    << std::endl;
  //}

  UQ_FATAL_TEST_MACRO(m_optionsInputFileName == "",
                      m_worldRank,
                      "uqBaseEnvironmentClass::scanInputFileForMyOptions()",
                      "m_optionsInputFileName is 'nothing'");
  //std::ifstream ifs(m_optionsInputFileName.c_str());
  std::ifstream* ifs = new std::ifstream(m_optionsInputFileName.c_str());
#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "in uqBaseEnvClass::scanInputFileForMyOptions(), before store(a)" << std::endl;
#endif

  UQ_FATAL_TEST_MACRO(m_allOptionsMap == NULL,
                      m_worldRank,
                      "uqBaseEnvironmentClass::scanInputFileForMyOptions()",
                      "m_allOptionsMap variable is NULL");
  po::store(po::parse_config_file(*ifs, *m_allOptionsDesc, true), *m_allOptionsMap);
#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "in uqBaseEnvClass::scanInputFileForMyOptions(), after store(a)" << std::endl;
#endif
  po::notify(*m_allOptionsMap);

  //ifs.close();
  delete ifs;
#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "Leaving uqBaseEnvClass::scanInputFileForMyOptions()" << std::endl;
#endif

  return;
}

std::string
uqBaseEnvironmentClass::optionsInputFileName() const
{
  if (m_optionsInputFileAccessState) {
    return m_optionsInputFileName;
  }
  else {
    return "";
  }
}

void
uqBaseEnvironmentClass::setOptionsInputFileAccessState(bool newState) const
{
  m_optionsInputFileAccessState = newState;

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
                      m_worldRank,
                      "uqBaseEnvironmentClass::allOptionsMap()",
                      "m_allOptionsMap variable is NULL");
  return *m_allOptionsMap;
}

unsigned int
uqBaseEnvironmentClass::displayVerbosity() const
{
  UQ_FATAL_TEST_MACRO(m_optionsObj == NULL,
                      m_worldRank,
                      "uqBaseEnvironmentClass::displayVerbosity()",
                      "m_optionsObj variable is NULL");
  return m_optionsObj->m_ov.m_displayVerbosity;
}

unsigned int
uqBaseEnvironmentClass::syncVerbosity() const
{
  UQ_FATAL_TEST_MACRO(m_optionsObj == NULL,
                      m_worldRank,
                      "uqBaseEnvironmentClass::displayVerbosity()",
                      "m_optionsObj variable is NULL");
  return m_optionsObj->m_ov.m_syncVerbosity;
}

unsigned int
uqBaseEnvironmentClass::checkingLevel() const
{
  UQ_FATAL_TEST_MACRO(m_optionsObj == NULL,
                      m_worldRank,
                      "uqBaseEnvironmentClass::checkingLevel()",
                      "m_optionsObj variable is NULL");
  return m_optionsObj->m_ov.m_checkingLevel;
}

const uqRngBaseClass*
uqBaseEnvironmentClass::rngObject() const
{
  return m_rngObject;
}

int
uqBaseEnvironmentClass::seed() const
{
  return m_rngObject->seed();
}

void
uqBaseEnvironmentClass::resetSeed(int newSeedOption)
{
  m_rngObject->resetSeed(newSeedOption);
  return;
}

std::string
uqBaseEnvironmentClass::platformName() const
{
  return m_optionsObj->m_ov.m_platformName;
}

std::string
uqBaseEnvironmentClass::identifyingString() const
{
  return m_optionsObj->m_ov.m_identifyingString;
}

void
uqBaseEnvironmentClass::resetIdentifyingString(const std::string& newString) const // Yes, const
{
  m_optionsObj->m_ov.m_identifyingString = newString;
  return;
}

struct timeval
uqBaseEnvironmentClass::timevalBegin() const
{
  return m_timevalBegin;
}

bool
uqBaseEnvironmentClass::openOutputFile(
  const std::string&            baseFileName,
  const std::string&            inputFileType,
  const std::set<unsigned int>& allowedSubEnvIds,
        bool                    writeOver,
        uqFilePtrSetStruct&     filePtrSet) const
{
  std::string fileType(inputFileType);
#ifdef QUESO_HAS_HDF5
  // Do nothing
#else
  if (fileType == UQ_FILE_EXTENSION_FOR_HDF_FORMAT) {
    if (m_subDisplayFile) {
      *this->subDisplayFile() << "WARNING in uqBaseEnvironmentClass::openOutputFile()"
                              << ": file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                              << "' has been requested, but this QUESO library has not been built with 'hdf5'"
                              << ". Code will therefore process the file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                              << "' instead..."
                              << std::endl;
    }
    if (this->subRank() == 0) {
      std::cerr << "WARNING in uqBaseEnvironmentClass::openOutputFile()"
                << ": file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                << "' has been requested, but this QUESO library has not been built with 'hdf5'"
                << ". Code will therefore process the file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                << "' instead..."
                << std::endl;
    }
    fileType = UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT;
  }
#endif

  bool returnValue = true;
  filePtrSet.ofsVar = NULL;
  if ((baseFileName                         == UQ_ENV_FILENAME_FOR_NO_OUTPUT_FILE) ||
      (allowedSubEnvIds.find(this->subId()) == allowedSubEnvIds.end()            )) {
    if ((m_subDisplayFile) && (this->displayVerbosity() > 10)) { // output debug
      *this->subDisplayFile() << "In uqBaseEnvironmentClass::openOutputFile()"
                              << ", subId = "     << this->subId()
                              << ": no output file opened with base name '" << baseFileName << "." << fileType
                              << "'"
                              << ", writeOver = " << writeOver
                              << std::endl;
    }
    returnValue = false;
  }
  else {
    //////////////////////////////////////////////////////////////////
    // Open file
    //////////////////////////////////////////////////////////////////
    if ((m_subDisplayFile) && (this->displayVerbosity() > 10)) { // output debug
      *this->subDisplayFile() << "In uqBaseEnvironmentClass::openOutputFile()"
                              << ", subId = "     << this->subId()
                              << ": opening output file with base name '" << baseFileName << "." << fileType
                              << "'"
                              << ", writeOver = " << writeOver
                              << std::endl;
    }

    if (this->subRank() == 0) {
#if 0
      std::cout << "In uqBaseEnvironmentClass::openOutputFile()"
                << ": opening output file with base name '" << baseFileName << "." << fileType
                << "'"
                << ", writeOver = " << writeOver
                << std::endl;
#endif

      ////////////////////////////////////////////////////////////////
      // Verify parent directory exists (for cases when a user
      // specifies a relative path for the desired output file).
      ////////////////////////////////////////////////////////////////
      // std::cout << "checking " << baseFileName+"_sub"+this->subIdString()+"."+fileType << std::endl;
      if ((m_subDisplayFile) && (this->displayVerbosity() > 10)) { // output debug
        *this->subDisplayFile() << "In uqBaseEnvironmentClass::openOutputFile()"
                                << ", subId = "     << this->subId()
                                << ", trying to open output file with base name '" << baseFileName << "." << fileType
                                << "'"
                                << ", writeOver = " << writeOver
                                << ": calling uqCheckFilePath()..." 
                                << std::endl;
      }
      int irtrn = uqCheckFilePath((baseFileName+"_sub"+this->subIdString()+"."+fileType).c_str());
      if ((m_subDisplayFile) && (this->displayVerbosity() > 10)) { // output debug
        *this->subDisplayFile() << "In uqBaseEnvironmentClass::openOutputFile()"
                                << ", subId = "     << this->subId()
                                << ", trying to open output file with base name '" << baseFileName << "." << fileType
                                << "'"
                                << ", writeOver = " << writeOver
                                << ": returned from uqCheckFilePath() with irtrn = " << irtrn
                                << std::endl;
      }
      UQ_FATAL_TEST_MACRO(irtrn < 0,
                          m_worldRank,
                          "uqBaseEnvironmentClass::openOutputFile()",
                          "unable to verify output path");

      if (writeOver) {
        //////////////////////////////////////////////////////////////
        // Write over an eventual pre-existing file
        //////////////////////////////////////////////////////////////
        if (fileType == UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT) {
          filePtrSet.ofsVar = new std::ofstream((baseFileName+"_sub"+this->subIdString()+"."+fileType).c_str(), 
                                                std::ofstream::out | std::ofstream::trunc);
        }
        else if (fileType == UQ_FILE_EXTENSION_FOR_HDF_FORMAT) {
          UQ_FATAL_TEST_MACRO(true,
                              m_worldRank,
                              "uqBaseEnvironmentClass::openOutputFile(), writeOver=true",
                              "hdf file type not supported yet");
        }
        else {
          UQ_FATAL_TEST_MACRO(true,
                              m_worldRank,
                              "uqBaseEnvironmentClass::openOutputFile(), writeOver=true",
                              "invalid file type");
        }
        if ((m_subDisplayFile) && (this->displayVerbosity() > 10)) { // output debug
          *this->subDisplayFile() << "In uqBaseEnvironmentClass::openOutputFile()"
                                  << ": just opened output file with base name '" << baseFileName << "." << fileType
                                  << "'"
                                  << ", writeOver = " << writeOver
                                  << ", options 'out|trunc'"
                                  << ", osfvar = " << filePtrSet.ofsVar
                                  << std::endl;
        }
      }
      else {
        //////////////////////////////////////////////////////////////
        // Write at the end of an eventual pre-existing file
        //////////////////////////////////////////////////////////////
        if (fileType == UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT) {
#if 0
          filePtrSet.ofsVar = new std::ofstream((baseFileName+"_sub"+this->subIdString()+"."+fileType).c_str(), 
                                                std::ofstream::in);
          std::cout << "filePtrSet.ofsVar(1) = " << filePtrSet.ofsVar << std::endl;
          if (filePtrSet.ofsVar) std::cout << "filePtrSet.ofsVar(1)->is_open() = " << filePtrSet.ofsVar->is_open() << std::endl;
          if (filePtrSet.ofsVar) delete filePtrSet.ofsVar;
#endif
          // 'uqm' and Ranger nodes behave differently on ofstream constructor... prudenci 2010/03/05
          filePtrSet.ofsVar = new std::ofstream((baseFileName+"_sub"+this->subIdString()+"."+fileType).c_str(), 
                                                std::ofstream::out /*| std::ofstream::in*/ | std::ofstream::app);
        }
        else if (fileType == UQ_FILE_EXTENSION_FOR_HDF_FORMAT) {
          UQ_FATAL_TEST_MACRO(true,
                              m_worldRank,
                              "uqBaseEnvironmentClass::openOutputFile(), writeOver=false",
                              "hdf file type not supported yet");
        }
        else {
          UQ_FATAL_TEST_MACRO(true,
                              m_worldRank,
                              "uqBaseEnvironmentClass::openOutputFile(), writeOver=false",
                              "invalid file type");
        }
        if ((m_subDisplayFile) && (this->displayVerbosity() > 10)) { // output debug
          *this->subDisplayFile() << "In uqBaseEnvironmentClass::openOutputFile()"
                                  << ": just opened output file with base name '" << baseFileName << "." << fileType
                                  << "'"
                                  << ", writeOver = " << writeOver
                                  << ", options 'out|in|app'"
                                  << ", osfvar = " << filePtrSet.ofsVar
                                  << std::endl;
        }

        if (fileType == UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT) {
          //std::cout << "filePtrSet.ofsVar(2) = " << filePtrSet.ofsVar << std::endl;
          //if (filePtrSet.ofsVar) std::cout << "filePtrSet.ofsVar(2)->is_open() = " << filePtrSet.ofsVar->is_open() << std::endl;
          if ((filePtrSet.ofsVar            == NULL ) ||
              (filePtrSet.ofsVar->is_open() == false)) {
            //std::cout << "Retrying 1..." << std::endl;
            delete filePtrSet.ofsVar;
            filePtrSet.ofsVar = new std::ofstream((baseFileName+"_sub"+this->subIdString()+"."+fileType).c_str(), 
                                                  std::ofstream::out | std::ofstream::trunc);
            if ((m_subDisplayFile) && (this->displayVerbosity() > 10)) { // output debug
              *this->subDisplayFile() << "In uqBaseEnvironmentClass::openOutputFile()"
                                      << ": just opened output file with base name '" << baseFileName << "." << fileType
                                      << "'"
                                      << ", writeOver = " << writeOver
                                      << ", options 'out|trunc'"
                                      << ", osfvar = " << filePtrSet.ofsVar
                                      << std::endl;
            }
          }
        } // only for matlab formats
      }
      if (filePtrSet.ofsVar != NULL) {
        if ((m_subDisplayFile) && (this->displayVerbosity() > 10)) { // output debug
          *this->subDisplayFile() << "In uqBaseEnvironmentClass::openOutputFile()"
                                  << ", subId = "     << this->subId()
                                  << ": succeeded on opening output file with base name '" << baseFileName << "." << fileType
                                  << "'"
                                  << ", writeOver = " << writeOver
                                  << std::endl;
        }
      }
      else {
        std::cerr << "In uqBaseEnvironmentClass::openOutputFile()"
                  << ": failed to open output file with base name '" << baseFileName << "." << fileType
                  << "'"
                  << std::endl;
      }
      UQ_FATAL_TEST_MACRO((filePtrSet.ofsVar && filePtrSet.ofsVar->is_open()) == false,
                          this->worldRank(),
                          "openOutputFile()",
                          "failed to open output file");
    }
    else {
      returnValue = false;
    }
    //this->subComm().Barrier(); // prudenci-2011-01-17
  }

  return returnValue;
}

bool
uqBaseEnvironmentClass::openUnifiedOutputFile(
  const std::string&        baseFileName,
  const std::string&        inputFileType,
        bool                writeOver,
        uqFilePtrSetStruct& filePtrSet) const
{
  std::string fileType(inputFileType);
#ifdef QUESO_HAS_HDF5
  // Do nothing
#else
  if (fileType == UQ_FILE_EXTENSION_FOR_HDF_FORMAT) {
    if (m_subDisplayFile) {
      *this->subDisplayFile() << "WARNING in uqBaseEnvironmentClass::openUnifiedOutputFile()"
                              << ": file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                              << "' has been requested, but this QUESO library has not been built with 'hdf5'"
                              << ". Code will therefore process the file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                              << "' instead..."
                              << std::endl;
    }
    if (this->subRank() == 0) {
      std::cerr << "WARNING in uqBaseEnvironmentClass::openUnifiedOutputFile()"
                << ": file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                << "' has been requested, but this QUESO library has not been built with 'hdf5'"
                << ". Code will therefore process the file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                << "' instead..."
                << std::endl;
    }
    fileType = UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT;
  }
#endif

  bool returnValue = true;
  filePtrSet.ofsVar = NULL;
  if (baseFileName == ".") {
    if ((m_subDisplayFile) && (this->displayVerbosity() > 10)) { // output debug
      *this->subDisplayFile() << "In uqBaseEnvironmentClass::openUnifiedOutputFile()"
                              << ": no unified output file opened with base name '" << baseFileName << "." << fileType
                              << "'"
                              << ", writeOver = " << writeOver
                              << std::endl;
    }
    returnValue = false;
  }
  else {
    //////////////////////////////////////////////////////////////////
    // Open file
    //////////////////////////////////////////////////////////////////
    if ((m_subDisplayFile) && (this->displayVerbosity() > 10)) { // output debug
      *this->subDisplayFile() << "In uqBaseEnvironmentClass::openUnifiedOutputFile()"
                              << ": opening unified output file with base name '" << baseFileName << "." << fileType
                              << "'"
                              << ", writeOver = " << writeOver
                              << std::endl;
    }


    //if ((this->subRank   () == 0) &&
    //    (this->inter0Rank() == 0)) {
#if 0
      std::cout << "In uqBaseEnvironmentClass::openUnifiedOutputFile()"
                << ": opening output file with base name '" << baseFileName << "." << fileType
                << "'"
                << ", writeOver = " << writeOver
                << std::endl;
#endif
      ////////////////////////////////////////////////////////////////
      // Verify parent directory exists (for cases when a user
      // specifies a relative path for the desired output file). prudenci 2010/06/26
      ////////////////////////////////////////////////////////////////
      // std::cout << "checking " << baseFileName+"."+fileType << std::endl;
      int irtrn = uqCheckFilePath((baseFileName+"."+fileType).c_str());
      UQ_FATAL_TEST_MACRO(irtrn < 0,
                          m_worldRank,
                          "uqBaseEnvironmentClass::openUnifiedOutputFile()",
                          "unable to verify output path");

      if (writeOver) {
        ////////////////////////////////////////////////////////////////
        // Write over an eventual pre-existing file
        ////////////////////////////////////////////////////////////////
        if (fileType == UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT) {
          filePtrSet.ofsVar = new std::ofstream((baseFileName+"."+fileType).c_str(),
                                                std::ofstream::out | std::ofstream::trunc);
        }
#ifdef QUESO_HAS_HDF5
        else if (fileType == UQ_FILE_EXTENSION_FOR_HDF_FORMAT) {
          filePtrSet.h5Var = H5Fcreate((baseFileName+"."+fileType).c_str(),
                                       H5F_ACC_TRUNC,
                                       H5P_DEFAULT,
                                       H5P_DEFAULT);
        }
#endif
        else {
          UQ_FATAL_TEST_MACRO(true,
                              m_worldRank,
                              "uqBaseEnvironmentClass::openUnifiedOutputFile(), writeOver=true",
                              "invalid file type");
        }
        if ((m_subDisplayFile) && (this->displayVerbosity() > 10)) { // output debug
          *this->subDisplayFile() << "In uqBaseEnvironmentClass::openUnifiedOutputFile()"
                                  << ": just opened output file with base name '" << baseFileName << "." << fileType
                                  << "'"
                                  << ", writeOver = " << writeOver
                                  << ", options 'out|trunc'"
                                  << ", osfvar = " << filePtrSet.ofsVar
                                  << std::endl;
        }
      }
      else {
        ////////////////////////////////////////////////////////////////
        // Write at the end of an eventual pre-existing file
        // 'uqm' and Ranger nodes behave differently on ofstream constructor... prudenci 2010/03/05
        ////////////////////////////////////////////////////////////////
        if (fileType == UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT) {
          filePtrSet.ofsVar = new std::ofstream((baseFileName+"."+fileType).c_str(),
                                                std::ofstream::out /*| std::ofstream::in*/ | std::ofstream::app);
        }
#ifdef QUESO_HAS_HDF5
        else if (fileType == UQ_FILE_EXTENSION_FOR_HDF_FORMAT) {
          filePtrSet.h5Var = H5Fcreate((baseFileName+"."+fileType).c_str(), // TEMPORARY - FIX ME
                                       H5F_ACC_TRUNC,
                                       H5P_DEFAULT,
                                       H5P_DEFAULT);
          //UQ_FATAL_TEST_MACRO(true,
          //                    m_worldRank,
          //                    "uqBaseEnvironmentClass::openUnifiedOutputFile(), writeOver=false",
          //                    "hdf file type not supported yet");
        }
#endif
        else {
          UQ_FATAL_TEST_MACRO(true,
                              m_worldRank,
                              "uqBaseEnvironmentClass::openUnifiedOutputFile(), writeOver=false",
                              "invalid file type");
        }
        if ((m_subDisplayFile) && (this->displayVerbosity() > 10)) { // output debug
          *this->subDisplayFile() << "In uqBaseEnvironmentClass::openUnifiedOutputFile()"
                                  << ": just opened output file with base name '" << baseFileName << "." << fileType
                                  << "'"
                                  << ", writeOver = " << writeOver
                                  << ", options 'out|in|app'"
                                  << ", osfvar = " << filePtrSet.ofsVar
                                  << std::endl;
        }
        if ((filePtrSet.ofsVar            == NULL ) ||
            (filePtrSet.ofsVar->is_open() == false)) {
#if 0
          std::cout << "Retrying 2..." << std::endl;
#endif
          delete filePtrSet.ofsVar;
          filePtrSet.ofsVar = new std::ofstream((baseFileName+"."+fileType).c_str(),
                                                std::ofstream::out | std::ofstream::trunc);
          if ((m_subDisplayFile) && (this->displayVerbosity() > 10)) { // output debug
            *this->subDisplayFile() << "In uqBaseEnvironmentClass::openUnifiedOutputFile()"
                                    << ": just opened output file with base name '" << baseFileName << "." << fileType
                                    << "'"
                                    << ", writeOver = " << writeOver
                                    << ", options 'out|trunc'"
                                    << ", osfvar = " << filePtrSet.ofsVar
                                    << std::endl;
          }
        }
      }
      if (filePtrSet.ofsVar == NULL) {
        std::cerr << "In uqBaseEnvironmentClass::openUnifiedOutputFile()"
                  << ": failed to open unified output file with base name '" << baseFileName << "." << fileType
                  << "'"
                  << std::endl;
      }
      UQ_FATAL_TEST_MACRO((filePtrSet.ofsVar && filePtrSet.ofsVar->is_open()) == false,
                          this->worldRank(),
                          "openUnifiedOutputFile()",
                          "failed to open output file");
    //}
  }

  return returnValue;
}

bool
uqBaseEnvironmentClass::openInputFile(
  const std::string&            baseFileName,
  const std::string&            inputFileType,
  const std::set<unsigned int>& allowedSubEnvIds,
        uqFilePtrSetStruct&     filePtrSet) const
{
  std::string fileType(inputFileType);
#ifdef QUESO_HAS_HDF5
  // Do nothing
#else
  if (fileType == UQ_FILE_EXTENSION_FOR_HDF_FORMAT) {
    if (m_subDisplayFile) {
      *this->subDisplayFile() << "WARNING in uqBaseEnvironmentClass::openInputFile()"
                              << ": file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                              << "' has been requested, but this QUESO library has not been built with 'hdf5'"
                              << ". Code will therefore process the file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                              << "' instead..."
                              << std::endl;
    }
    if (this->subRank() == 0) {
      std::cerr << "WARNING in uqBaseEnvironmentClass::openInputFile()"
                << ": file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                << "' has been requested, but this QUESO library has not been built with 'hdf5'"
                << ". Code will therefore process the file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                << "' instead..."
                << std::endl;
    }
    fileType = UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT;
  }
#endif

  bool returnValue = true;
  filePtrSet.ifsVar = NULL;
  if ((baseFileName                         == UQ_ENV_FILENAME_FOR_NO_INPUT_FILE) ||
      (allowedSubEnvIds.find(this->subId()) == allowedSubEnvIds.end()           )) {
    if ((m_subDisplayFile) && (this->displayVerbosity() >= 10)) {
      *this->subDisplayFile() << "In uqBaseEnvironmentClass::openInputFile()"
                              << ": no input file opened with base name '" << baseFileName << "." << fileType
                              << "'"
                              << std::endl;
    }
    returnValue = false;
  }
  else {
    //////////////////////////////////////////////////////////////////
    // Open file
    //////////////////////////////////////////////////////////////////
    if ((m_subDisplayFile) && (this->displayVerbosity() >= 10)) {
      *this->subDisplayFile() << "In uqBaseEnvironmentClass::openInputFile()"
                              << ": opening input file with base name '" << baseFileName << "." << fileType
                              << "'"
                              << std::endl;
    }
    if (this->subRank() == 0) {
      ////////////////////////////////////////////////////////////////
      // Verify parent directory exists (for cases when a user
      // specifies a relative path for the desired output file). prudenci 2010/06/26
      ////////////////////////////////////////////////////////////////
      // std::cout << "checking " << baseFileName+"."+fileType << std::endl;
      int irtrn = uqCheckFilePath((baseFileName+"."+fileType).c_str());
      UQ_FATAL_TEST_MACRO(irtrn < 0,
                          m_worldRank,
                          "uqBaseEnvironmentClass::openInputFile()",
                          "unable to verify input path");

      if (fileType == UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT) {
        filePtrSet.ifsVar = new std::ifstream((baseFileName+"."+fileType).c_str(),
                                              std::ofstream::in);
        if ((filePtrSet.ifsVar == NULL) || (filePtrSet.ifsVar->is_open() == false)) {
          std::cerr << "In uqBaseEnvironmentClass::openInputFile()"
                    << ": failed to open input file with base name '" << baseFileName << "." << fileType
                    << "'"
                    << std::endl;
        }
        UQ_FATAL_TEST_MACRO((filePtrSet.ifsVar == NULL) || (filePtrSet.ifsVar->is_open() == false),
                            this->worldRank(),
                            "uqBaseEnvironmentClass::openInputFile()",
                            "file with fileName could not be found");
      }
#ifdef QUESO_HAS_HDF5
      else if (fileType == UQ_FILE_EXTENSION_FOR_HDF_FORMAT) {
        filePtrSet.h5Var = H5Fopen((baseFileName+"."+fileType).c_str(),
                                   H5F_ACC_RDONLY,
                                   H5P_DEFAULT);
      }
#endif
      else {
        UQ_FATAL_TEST_MACRO(true,
                            m_worldRank,
                            "uqBaseEnvironmentClass::openInputFile()",
                            "invalid file type");
      }
    }
    else {
      returnValue = false;
    }
    //this->subComm().Barrier(); // prudenci-2011-01-17
  }

  return returnValue;
}

bool
uqBaseEnvironmentClass::openUnifiedInputFile(
  const std::string&        baseFileName,
  const std::string&        inputFileType,
        uqFilePtrSetStruct& filePtrSet) const
{
  std::string fileType(inputFileType);
#ifdef QUESO_HAS_HDF5
  // Do nothing
#else
  if (fileType == UQ_FILE_EXTENSION_FOR_HDF_FORMAT) {
    if (m_subDisplayFile) {
      *this->subDisplayFile() << "WARNING in uqBaseEnvironmentClass::openUnifiedInputFile()"
                              << ": file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                              << "' has been requested, but this QUESO library has not been built with 'hdf5'"
                              << ". Code will therefore process the file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                              << "' instead..."
                              << std::endl;
    }
    if (this->subRank() == 0) {
      std::cerr << "WARNING in uqBaseEnvironmentClass::openUnifiedInputFile()"
                << ": file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                << "' has been requested, but this QUESO library has not been built with 'hdf5'"
                << ". Code will therefore process the file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                << "' instead..."
                << std::endl;
    }
    fileType = UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT;
  }
#endif

  bool returnValue = true;
  filePtrSet.ifsVar = NULL;
  if (baseFileName == UQ_ENV_FILENAME_FOR_NO_INPUT_FILE) {
    if ((m_subDisplayFile) && (this->displayVerbosity() >= 10)) {
      *this->subDisplayFile() << "In uqBaseEnvironmentClass::openUnifiedInputFile()"
                              << ": no input file opened with base name '" << baseFileName << "." << fileType
                              << "'"
                              << std::endl;
    }
    returnValue = false;
  }
  else {
    //////////////////////////////////////////////////////////////////
    // Open file
    //////////////////////////////////////////////////////////////////
    if ((m_subDisplayFile) && (this->displayVerbosity() >= 10)) {
      *this->subDisplayFile() << "In uqBaseEnvironmentClass::openUnifiedInputFile()"
                              << ": opening input file with base name '" << baseFileName << "." << fileType
                              << "'"
                              << std::endl;
    }
    if (this->subRank() == 0) { // Needed ???????? prudenci 2010-11-11
      ////////////////////////////////////////////////////////////////
      // Verify parent directory exists (for cases when a user
      // specifies a relative path for the desired output file). prudenci 2010/06/26
      ////////////////////////////////////////////////////////////////
      // std::cout << "checking " << baseFileName+"."+fileType << std::endl;
      int irtrn = uqCheckFilePath((baseFileName+"."+fileType).c_str());
      UQ_FATAL_TEST_MACRO(irtrn < 0,
                          m_worldRank,
                          "uqBaseEnvironmentClass::openUnifiedInputFile()",
                          "unable to verify input path");

      if (fileType == UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT) {
        filePtrSet.ifsVar = new std::ifstream((baseFileName+"."+fileType).c_str(),
                                              std::ofstream::in);
        if ((filePtrSet.ifsVar == NULL) || (filePtrSet.ifsVar->is_open() == false)) {
          std::cerr << "In uqBaseEnvironmentClass::openUnifiedInputFile()"
                    << ": failed to open input file with base name '" << baseFileName << "." << fileType
                    << "'"
                    << std::endl;
        }
        UQ_FATAL_TEST_MACRO((filePtrSet.ifsVar == NULL) || (filePtrSet.ifsVar->is_open() == false),
                            this->worldRank(),
                            "uqBaseEnvironmentClass::openUnifiedInputFile()",
                            "file with fileName could not be found");
      }
#ifdef QUESO_HAS_HDF5
      else if (fileType == UQ_FILE_EXTENSION_FOR_HDF_FORMAT) {
        filePtrSet.h5Var = H5Fopen((baseFileName+"."+fileType).c_str(),
                                   H5F_ACC_RDONLY,
                                   H5P_DEFAULT);
      }
#endif
      else {
        UQ_FATAL_TEST_MACRO(true,
                            m_worldRank,
                            "uqBaseEnvironmentClass::openUnifiedInputFile()",
                            "invalid file type");
      }
    }
    //else {
    //  returnValue = false;
    //}
    //this->subComm().Barrier();
  }

  return returnValue;
}

void
uqBaseEnvironmentClass::closeFile(
  uqFilePtrSetStruct& filePtrSet,
  const std::string&  inputFileType) const
{
  std::string fileType(inputFileType);
#ifdef QUESO_HAS_HDF5
  // Do nothing
#else
  if (fileType == UQ_FILE_EXTENSION_FOR_HDF_FORMAT) {
    if (m_subDisplayFile) {
      *this->subDisplayFile() << "WARNING in uqBaseEnvironmentClass::closeFile()"
                              << ": file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                              << "' has been requested, but this QUESO library has not been built with 'hdf5'"
                              << ". Code will therefore process the file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                              << "' instead..."
                              << std::endl;
    }
    if (this->subRank() == 0) {
      std::cerr << "WARNING in uqBaseEnvironmentClass::closeFile()"
                << ": file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                << "' has been requested, but this QUESO library has not been built with 'hdf5'"
                << ". Code will therefore process the file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                << "' instead..."
                << std::endl;
    }
    fileType = UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT;
  }
#endif

  if (fileType == UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT) {
    //filePtrSet.ofsVar->close(); // close() crashes on Mac; need to use delete(); why? prudenci 2010/June
    delete filePtrSet.ofsVar;
    filePtrSet.ofsVar = NULL;

    //filePtrSet.ifsVar->close();
    delete filePtrSet.ifsVar;
    filePtrSet.ifsVar = NULL;
  }
#ifdef QUESO_HAS_HDF5
  else if (fileType == UQ_FILE_EXTENSION_FOR_HDF_FORMAT) {
    H5Fclose(filePtrSet.h5Var);
  }
#endif
  else {
    UQ_FATAL_TEST_MACRO(true,
                        m_worldRank,
                        "uqBaseEnvironmentClass::closeFile()",
                        "invalid file type");
  }

  return;
}

void
uqBaseEnvironmentClass::setExceptionalCircunstance(bool value) const
{
  m_exceptionalCircunstance = value;
  return;
}

bool
uqBaseEnvironmentClass::exceptionalCircunstance() const
{
  return m_exceptionalCircunstance;
}


//*****************************************************
// Empty Environment
//*****************************************************
uqEmptyEnvironmentClass::uqEmptyEnvironmentClass()
  :
  uqBaseEnvironmentClass("",NULL)
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
  uqRawType_MPI_Comm             inputComm,
  const char*                    passedOptionsInputFileName,
  const char*                    prefix,
  const uqEnvOptionsValuesClass* alternativeOptionsValues)
  :
  uqBaseEnvironmentClass(passedOptionsInputFileName,alternativeOptionsValues)
{
#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "Entering uqFullEnvClass" << std::endl;
#endif

  //////////////////////////////////////////////////
  // Initialize "full" communicator
  //////////////////////////////////////////////////
#ifdef QUESO_HAS_MPI
  int mpiRC = MPI_Comm_rank(MPI_COMM_WORLD,&m_worldRank);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      UQ_UNAVAILABLE_RANK,
                      "uqFullEnvironmentClass::commonConstructor()",
                      "failed to get world fullRank()");
#else
  m_worldRank = 0;
#endif

  m_fullComm = new uqMpiCommClass(*this,inputComm);
  m_fullRank     = m_fullComm->MyPID();
  m_fullCommSize = m_fullComm->NumProc();
#ifdef QUESO_HAS_MPI
  mpiRC = MPI_Comm_group(m_fullComm->Comm(), &m_fullGroup);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_worldRank,
                      "uqFullEnvironmentClass::commonConstructor()",
                      "failed MPI_Comm_group()");
#else
  m_fullGroup = 0;
#endif

#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "In uqFullEnvClass, finished dealing with MPI initially" << std::endl;
#endif

  //////////////////////////////////////////////////
  // Read options
  //////////////////////////////////////////////////
  if (m_optionsInputFileName == "") {
    m_optionsObj = new uqEnvironmentOptionsClass(*this,prefix,m_alternativeOptionsValues);
  }
  else {
    m_allOptionsMap  = new po::variables_map();
    m_allOptionsDesc = new po::options_description("Allowed options");
    m_optionsObj = new uqEnvironmentOptionsClass(*this,prefix);

    readOptionsInputFile();

    m_optionsObj->scanOptionsValues();
  }

#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "In uqFullEnvClass, finished scanning options" << std::endl;
#endif

  // Only display these messages if the user wants them
  // NOTE: This got moved below the Read Options section
  // because we need the options to be read to know what
  // the verbosity level is.
  if (this->displayVerbosity() > 0) {
    //////////////////////////////////////////////////
    // Display main initial messages
    // 'std::cout' is for: main trace messages + synchronized trace messages + error messages prior to 'exit()' or 'abort()'
    //////////////////////////////////////////////////
    /*int iRC = 0;*/
    /*iRC = */gettimeofday(&m_timevalBegin, NULL);
      
    if (m_fullRank == 0) {
      QUESO::QUESO_version_print(std::cout);
    }
      
    if (m_fullRank == 0) {
      std::cout << "Beginning run at " << ctime(&m_timevalBegin.tv_sec)
                << std::endl;
    }
  }

  //////////////////////////////////////////////////
  // Deal with multiple subEnvironments: create the sub communicators, one for each subEnvironment
  //////////////////////////////////////////////////
  unsigned int numRanksPerSubEnvironment = m_fullCommSize/m_optionsObj->m_ov.m_numSubEnvironments;

  m_subId = m_fullRank/numRanksPerSubEnvironment;
  char tmpSubId[16];
  sprintf(tmpSubId,"%u",m_subId);
  m_subIdString = tmpSubId;

  if (m_optionsObj->m_ov.m_subDisplayAllowAll) {
    m_optionsObj->m_ov.m_subDisplayAllowedSet.insert((unsigned int) m_subId);
  }

  std::vector<int> fullRanksOfMySubEnvironment(numRanksPerSubEnvironment,0);
  for (unsigned int i = 0; i < numRanksPerSubEnvironment; ++i) {
    fullRanksOfMySubEnvironment[i] = m_subId * numRanksPerSubEnvironment + i;
  }
#ifdef QUESO_HAS_MPI
  mpiRC = MPI_Group_incl(m_fullGroup, (int) numRanksPerSubEnvironment, &fullRanksOfMySubEnvironment[0], &m_subGroup);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_worldRank,
                      "uqFullEnvironmentClass::commonConstructor()",
                      "failed MPI_Group_incl() for a subEnvironment");
#else
  m_subGroup = 0;
#endif
  uqRawType_MPI_Comm subRawComm;
#ifdef QUESO_HAS_MPI
  mpiRC = MPI_Comm_create(m_fullComm->Comm(), m_subGroup, &subRawComm);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_worldRank,
                      "uqFullEnvironmentClass::commonConstructor()",
                      "failed MPI_Comm_group() for a subEnvironment");
#else
  subRawComm = 0;
#endif
  m_subComm = new uqMpiCommClass(*this,subRawComm);
  m_subRank     = m_subComm->MyPID();
  m_subCommSize = m_subComm->NumProc();

  //////////////////////////////////////////////////
  // Deal with multiple subEnvironments: create the self communicator
  //////////////////////////////////////////////////
  m_selfComm = new uqMpiCommClass(*this,uqRawValue_MPI_COMM_SELF);

  //////////////////////////////////////////////////
  // Deal with multiple subEnvironments: create the inter0 communicator
  //////////////////////////////////////////////////
  std::vector<int> fullRanksOfInter0(m_optionsObj->m_ov.m_numSubEnvironments,0);
  for (unsigned int i = 0; i < m_optionsObj->m_ov.m_numSubEnvironments; ++i) {
    fullRanksOfInter0[i] = i * numRanksPerSubEnvironment;
  }
#ifdef QUESO_HAS_MPI
  mpiRC = MPI_Group_incl(m_fullGroup, (int) m_optionsObj->m_ov.m_numSubEnvironments, &fullRanksOfInter0[0], &m_inter0Group);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_worldRank,
                      "uqFullEnvironmentClass::commonConstructor()",
                      "failed MPI_Group_incl() for inter0");
#else
  m_inter0Group = 0;
#endif
  uqRawType_MPI_Comm inter0RawComm;
#ifdef QUESO_HAS_MPI
  mpiRC = MPI_Comm_create(m_fullComm->Comm(), m_inter0Group, &inter0RawComm);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      m_worldRank,
                      "uqFullEnvironmentClass::commonConstructor()",
                      "failed MPI_Comm_group() for inter0");
#else
  inter0RawComm = 0;
#endif
  if (m_fullRank%numRanksPerSubEnvironment == 0) {
    m_inter0Comm = new uqMpiCommClass(*this,inter0RawComm);
    m_inter0Rank     = m_inter0Comm->MyPID();
    m_inter0CommSize = m_inter0Comm->NumProc();
  }

  if (m_optionsObj->m_ov.m_subDisplayAllowAll) {
    // This situation has been already taken care of above
  }
  else if (m_optionsObj->m_ov.m_subDisplayAllowInter0) {
    if (m_inter0Rank >= 0) {
      m_optionsObj->m_ov.m_subDisplayAllowedSet.insert((unsigned int) m_subId);
    }
  }


  //////////////////////////////////////////////////
  // Open "screen" file
  //////////////////////////////////////////////////
  bool openFile = false;
  if ((m_subRank                                               == 0                                              ) &&
      (m_optionsObj->m_ov.m_subDisplayFileName                 != UQ_ENV_FILENAME_FOR_NO_OUTPUT_FILE             ) &&
      (m_optionsObj->m_ov.m_subDisplayAllowedSet.find(m_subId) != m_optionsObj->m_ov.m_subDisplayAllowedSet.end())) {
    openFile = true;
  }

  if (openFile && m_worldRank == 0) {
    //////////////////////////////////////////////////////////////////
    // Verify parent directory exists (for cases when a user
    // specifies a relative path for the desired output file).
    //////////////////////////////////////////////////////////////////
    int irtrn = uqCheckFilePath((m_optionsObj->m_ov.m_subDisplayFileName+"_sub"+m_subIdString+".txt").c_str());
    UQ_FATAL_TEST_MACRO(irtrn < 0,
                        m_worldRank,
                        "uqEnvironment::constructor()",
                        "unable to verify output path");
  }

  ////////////////////////////////////////////////////////////////////
  // Ensure that rank 0 has created path, if necessary, before other tasks use it
  ////////////////////////////////////////////////////////////////////
  m_fullComm->Barrier(); 

  if (openFile) {
    //////////////////////////////////////////////////////////////////
    // Always write over an eventual pre-existing file
    //////////////////////////////////////////////////////////////////
    m_subDisplayFile = new std::ofstream((m_optionsObj->m_ov.m_subDisplayFileName+"_sub"+m_subIdString+".txt").c_str(),
                                         std::ofstream::out | std::ofstream::trunc);
    UQ_FATAL_TEST_MACRO((m_subDisplayFile && m_subDisplayFile->is_open()) == false,
                        m_worldRank,
                        "uqEnvironment::constructor()",
                        "failed to open sub screen file");

    QUESO::QUESO_version_print(*m_subDisplayFile);

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
    //                                     << std::endl;
    //sleep(3);
  }

  //////////////////////////////////////////////////
  // Deal with seed
  //////////////////////////////////////////////////
  if (m_optionsObj->m_ov.m_rngType == "gsl") {
    m_rngObject = new uqRngGslClass(m_optionsObj->m_ov.m_seed,m_worldRank);
  }
  else if (m_optionsObj->m_ov.m_rngType == "boost") {
    m_rngObject = new uqRngBoostClass(m_optionsObj->m_ov.m_seed,m_worldRank);
  }
  else {
    std::cerr << "In uqEnvironment::constructor()"
              << ": rngType = " << m_optionsObj->m_ov.m_rngType
              << std::endl;
    UQ_FATAL_TEST_MACRO(true,
                        m_worldRank,
                        "uqEnvironment::constructor()",
                        "the requested 'rngType' is not supported yet");
  }

  //////////////////////////////////////////////////
  // Leave commonConstructor()
  //////////////////////////////////////////////////
  m_fullComm->Barrier();
  m_fullEnvIsReady = true;

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
                                   << ": name of file is '" << m_optionsInputFileName.c_str() << "'"
                                   << std::endl;

    if (m_fullRank == 0) std::cout << "\nThis is a help message of the QUESO library."
                                   << "\nAn application using the QUESO library shall be executed by typing"
                                   << "\n  '<eventual mpi commands and options> <uqApplication> <uqInputFile>'"
                                   << "\nin the command line."
                                   << "\n"
                                   << std::endl;
#ifdef QUESO_HAS_MPI
    /*int mpiRC = 0;*/
    /*mpiRC = */MPI_Abort(m_fullComm->Comm(),-999);
#endif
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

