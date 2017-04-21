//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2017 The PECOS Development Team
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

#include <queso/Environment.h>

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

#include <queso/getpot.h>

#include <queso/config_queso.h>
#include <queso/EnvironmentOptions.h>
#include <queso/RngGsl.h>
#include <queso/RngBoost.h>
#include <queso/RngCXX11.h>
#include <queso/BasicPdfsGsl.h>
#include <queso/BasicPdfsBoost.h>
#include <queso/BasicPdfsCXX11.h>
#include <queso/Miscellaneous.h>
#include <sys/time.h>
#ifdef HAVE_GRVY
#include <grvy.h>
#endif

// queso error handling
#include <queso/asserts.h>

//-----------------------------
// Library versioning routines
//-----------------------------

namespace QUESO {

  void QUESO_version_print(std::ostream &os)
  {
    {
      os << "------------------------------------------------------------------------------------------" ;
      os << "--------------------" << std::endl;
      os << "QUESO Library: Version = " << QUESO_VERSION;
      os << " (" << QUESO_get_numeric_version() << ")" << std::endl << std::endl;

      os << QUESO_BUILD_DEVSTATUS << std::endl << std::endl;

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

FilePtrSetStruct::FilePtrSetStruct()
  :
  ofsVar(NULL),
  ifsVar(NULL)
#ifdef QUESO_HAS_HDF5
  ,  // lol
  h5Var(-1)
#endif
{
  queso_deprecated();
}

FilePtrSetStruct::~FilePtrSetStruct()
{
  queso_deprecated();
}

  //
  // queso terminate handler will be invoked for unhandled exceptions
  // needs to be invoked in FullEnvironment constructor
  //
  std::terminate_handler old_terminate_handler;

//*****************************************************
// Base class
//*****************************************************
// Default constructor --------------------------------
BaseEnvironment::BaseEnvironment(
  const char*                    passedOptionsInputFileName,
  EnvOptionsValues* alternativeOptionsValues)
  :
  m_fullEnvIsReady             (false),
  m_worldRank                  (-1),
  m_fullComm                   (),
  m_fullRank                   (-1),
  m_fullCommSize               (1),
  m_optionsInputFileName       (""),
  m_optionsInputFileAccessState(true),
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  m_allOptionsDesc             (),
  m_allOptionsMap              (),
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
  m_input                      (new GetPot),
  m_subComm                    (),
  m_subRank                    (-1),
  m_subCommSize                (1),
  m_selfComm                   (),
  m_inter0Comm                 (),
  m_inter0Rank                 (-1),
  m_inter0CommSize             (1),
  m_subDisplayFile             (),
  m_rngObject                  (),
  m_basicPdfs                  (),
  m_exceptionalCircumstance    (false),
  m_optionsObj                 ()
{
  if (passedOptionsInputFileName) m_optionsInputFileName     = passedOptionsInputFileName;

  // If the user passed in an options object pointer, we really shouldn't let
  // ScopedPtr delete their object, so we make a copy.  That way, the dtor
  // will kill this local copy and leave the user's object in tact.
  if (alternativeOptionsValues != NULL) {
    m_optionsObj.reset(new EnvOptionsValues(*alternativeOptionsValues));
  }
}

BaseEnvironment::BaseEnvironment(
  const std::string&             passedOptionsInputFileName,
  EnvOptionsValues* alternativeOptionsValues)
  :
  m_fullEnvIsReady             (false),
  m_worldRank                  (-1),
  m_fullComm                   (),
  m_fullRank                   (-1),
  m_fullCommSize               (1),
  m_optionsInputFileName       (passedOptionsInputFileName),
  m_optionsInputFileAccessState(true),
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  m_allOptionsDesc             (),
  m_allOptionsMap              (),
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
  m_input                      (new GetPot),
  m_subComm                    (),
  m_subRank                    (-1),
  m_subCommSize                (1),
  m_selfComm                   (),
  m_inter0Comm                 (),
  m_inter0Rank                 (-1),
  m_inter0CommSize             (1),
  m_subDisplayFile             (),
  m_rngObject                  (),
  m_basicPdfs                  (),
  m_exceptionalCircumstance    (false),
  m_optionsObj                 ()
{
  // If the user passed in an options object pointer, we really shouldn't let
  // ScopedPtr delete their object, so we make a copy.  That way, the dtor
  // will kill this local copy and leave the user's object in tact.
  if (alternativeOptionsValues != NULL) {
    m_optionsObj.reset(new EnvOptionsValues(*alternativeOptionsValues));
  }
}

// Destructor -------------------------------------------
BaseEnvironment::~BaseEnvironment()
{
  //if (m_subDisplayFile) {
  //  *m_subDisplayFile << "Entering BaseEnvironment::destructor()"
  //                          << std::endl;
  //}

  struct timeval timevalNow;
  /*int iRC = 0;*/
  /*iRC = */gettimeofday(&timevalNow, NULL);

  if (this->displayVerbosity() > 0) {
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

  //if (m_subDisplayFile) {
  //  *m_subDisplayFile << "Leaving BaseEnvironment::destructor()"
  //                          << std::endl;
  //}
}
// Environment, Communicator and Options Input File methods
bool
BaseEnvironment::fullEnvIsReady() const
{
  return m_fullEnvIsReady;
}
//-------------------------------------------------------
int
BaseEnvironment::worldRank() const
{
  return m_worldRank;
}
//-------------------------------------------------------
int
BaseEnvironment::fullRank() const
{
  return m_fullRank;
}
//-------------------------------------------------------
const MpiComm&
BaseEnvironment::fullComm() const
{
  queso_require_msg(m_fullComm, "m_fullComm variable is NULL");
  return *m_fullComm;
}
//-------------------------------------------------------
RawType_MPI_Group
BaseEnvironment::subGroup() const
{
  return m_subGroup;
}
//-------------------------------------------------------
int
BaseEnvironment::subRank() const
{
  return m_subRank;
}
//-------------------------------------------------------
const MpiComm&
BaseEnvironment::subComm() const
{
  queso_require_msg(m_subComm, "m_subComm variable is NULL");
  return *m_subComm;
}
//-------------------------------------------------------
const MpiComm&
BaseEnvironment::selfComm() const
{
  queso_require_msg(m_selfComm, "m_selfComm variable is NULL");
  return *m_selfComm;
}
//-------------------------------------------------------
int
BaseEnvironment::inter0Rank() const
{
  return m_inter0Rank;
}
//-------------------------------------------------------
const MpiComm&
BaseEnvironment::inter0Comm() const
{
  queso_require_msg(m_inter0Comm, "m_inter0Comm variable is NULL");
  return *m_inter0Comm;
}
//-------------------------------------------------------
std::ofstream*
BaseEnvironment::subDisplayFile() const
{
  // Potentially dangerous?  The user might delete it...
  return m_subDisplayFile.get();
}
//-------------------------------------------------------
std::string
BaseEnvironment::subDisplayFileName() const
{
  if (m_optionsObj == NULL) return ".";

  return m_optionsObj->m_subDisplayFileName;
}
//-------------------------------------------------------
unsigned int
BaseEnvironment::numSubEnvironments() const
{
  queso_require_msg(m_optionsObj, "m_optionsObj variable is NULL");
  return m_optionsObj->m_numSubEnvironments;
}
//-------------------------------------------------------
unsigned int
BaseEnvironment::subId() const
{
  return m_subId;
}
//-------------------------------------------------------
const std::string&
BaseEnvironment::subIdString() const
{
  return m_subIdString;
}
//-------------------------------------------------------
std::string
BaseEnvironment::optionsInputFileName() const
{
  if (m_optionsInputFileAccessState) {
    return m_optionsInputFileName;
  }
  else {
    return "";
  }
}
//-------------------------------------------------------
void
BaseEnvironment::setOptionsInputFileAccessState(bool newState) const
{
  m_optionsInputFileAccessState = newState;

  return;
}
//-------------------------------------------------------
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
#ifdef UQ_USES_COMMAND_LINE_OPTIONS
const boost::program_options::options_description&
BaseEnvironment::allOptionsDesc() const
{
  queso_deprecated();

  return *m_allOptionsDesc;
}
#endif
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
//-------------------------------------------------------
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
boost::program_options::variables_map&
BaseEnvironment::allOptionsMap() const
{
  queso_deprecated();

  queso_require_msg(m_allOptionsMap, "m_allOptionsMap variable is NULL");
  return *m_allOptionsMap;
}
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
//-------------------------------------------------------
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
void
BaseEnvironment::scanInputFileForMyOptions(const boost::program_options::options_description& optionsDesc) const
{
  queso_deprecated();

#ifdef UQ_USES_COMMAND_LINE_OPTIONS
  // If you want to use command line options, the following line does *not* work outside 'main.C',
  // e.g., in the constructor of FullEnvironment:
  // Line: boost::program_options::store(boost::program_options::parse_command_line(argc, argv, *m_allOptionsDesc), *m_allOptionsMap);
  //
  // Instead, put the following three lines *immediately after* instantianting the UQ environment
  // variable "FullEnvironment* env":
  // Line 1: boost::program_options::store(boost::program_options::parse_command_line(argc, argv, env->allOptionsDesc()), env->allOptionsMap());
  // Line 2: boost::program_options::notify(env->allOptionsMap());
  // Line 3: env->getMyOptionValues();
#endif

#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "Entering BaseEnv::scanInputFileForMyOptions()" << std::endl;
#endif

  queso_require_msg(m_allOptionsDesc, "m_allOptionsDesc variable is NULL");
  m_allOptionsDesc->add(optionsDesc);
  //if (m_subDisplayFile) {
  //  *m_subDisplayFile << *m_allOptionsDesc
  //                    << std::endl;
  //}

  queso_require_not_equal_to_msg(m_optionsInputFileName, std::string(""),
                                 std::string("m_optionsInputFileName is 'nothing'"));
  //std::ifstream ifs(m_optionsInputFileName.c_str());
  std::ifstream* ifs = new std::ifstream(m_optionsInputFileName.c_str());
#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "in BaseEnv::scanInputFileForMyOptions(), before store(a)" << std::endl;
#endif

  queso_require_msg(m_allOptionsMap, "m_allOptionsMap variable is NULL");
  boost::program_options::store(boost::program_options::parse_config_file(*ifs, *m_allOptionsDesc, true), *m_allOptionsMap);
#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "in BaseEnv::scanInputFileForMyOptions(), after store(a)" << std::endl;
#endif
  boost::program_options::notify(*m_allOptionsMap);

  //ifs.close();
  delete ifs;
#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "Leaving BaseEnv::scanInputFileForMyOptions()" << std::endl;
#endif

  return;
}
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
//-----------------------------------------------------
unsigned int
BaseEnvironment::displayVerbosity() const
{
  queso_require_msg(m_optionsObj, "m_optionsObj variable is NULL");
  return m_optionsObj->m_displayVerbosity;
}
//-------------------------------------------------------
unsigned int
BaseEnvironment::syncVerbosity() const
{
  queso_require_msg(m_optionsObj, "m_optionsObj variable is NULL");
  return m_optionsObj->m_syncVerbosity;
}
//-------------------------------------------------------
unsigned int
BaseEnvironment::checkingLevel() const
{
  queso_require_msg(m_optionsObj, "m_optionsObj variable is NULL");
  return m_optionsObj->m_checkingLevel;
}
//-------------------------------------------------------
const RngBase*
BaseEnvironment::rngObject() const
{
  return m_rngObject.get();
}
//-------------------------------------------------------
int
BaseEnvironment::seed() const
{
  return m_rngObject->seed();
}
//-------------------------------------------------------
void
BaseEnvironment::resetSeed(int newSeedOption)
{
  m_rngObject->resetSeed(newSeedOption);
  return;
}
//-------------------------------------------------------
const BasicPdfsBase*
BaseEnvironment::basicPdfs() const
{
  return m_basicPdfs.get();
}
//-------------------------------------------------------
std::string
BaseEnvironment::platformName() const
{
  return m_optionsObj->m_platformName;
}
//-------------------------------------------------------
std::string
BaseEnvironment::identifyingString() const
{
  return m_optionsObj->m_identifyingString;
}
//-------------------------------------------------------
void
BaseEnvironment::resetIdentifyingString(const std::string& newString)
{
  m_optionsObj->m_identifyingString = newString;
  return;
}
//-------------------------------------------------------
struct timeval
BaseEnvironment::timevalBegin() const
{
  return m_timevalBegin;
}
//-------------------------------------------------------
bool
BaseEnvironment::openOutputFile(
  const std::string&            baseFileName,
  const std::string&            inputFileType,
  const std::set<unsigned int>& allowedSubEnvIds,
        bool                    writeOver,
        FilePtrSetStruct&     filePtrSet) const
{
  std::string fileType(inputFileType);
#ifdef QUESO_HAS_HDF5
  // Do nothing
#else
  if (fileType == UQ_FILE_EXTENSION_FOR_HDF_FORMAT) {
    if (m_subDisplayFile) {
      *this->subDisplayFile() << "WARNING in BaseEnvironment::openOutputFile()"
                              << ": file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                              << "' has been requested, but this QUESO library has not been built with 'hdf5'"
                              << ". Code will therefore process the file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                              << "' instead..."
                              << std::endl;
    }
    if (this->subRank() == 0) {
      std::cerr << "WARNING in BaseEnvironment::openOutputFile()"
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
      *this->subDisplayFile() << "In BaseEnvironment::openOutputFile()"
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
      *this->subDisplayFile() << "In BaseEnvironment::openOutputFile()"
                              << ", subId = "     << this->subId()
                              << ": opening output file with base name '" << baseFileName << "." << fileType
                              << "'"
                              << ", writeOver = " << writeOver
                              << std::endl;
    }

    if (this->subRank() == 0) {
#if 0
      std::cout << "In BaseEnvironment::openOutputFile()"
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
        *this->subDisplayFile() << "In BaseEnvironment::openOutputFile()"
                                << ", subId = "     << this->subId()
                                << ", trying to open output file with base name '" << baseFileName << "." << fileType
                                << "'"
                                << ", writeOver = " << writeOver
                                << ": calling CheckFilePath()..."
                                << std::endl;
      }
      int irtrn = CheckFilePath((baseFileName+"_sub"+this->subIdString()+"."+fileType).c_str());
      if ((m_subDisplayFile) && (this->displayVerbosity() > 10)) { // output debug
        *this->subDisplayFile() << "In BaseEnvironment::openOutputFile()"
                                << ", subId = "     << this->subId()
                                << ", trying to open output file with base name '" << baseFileName << "." << fileType
                                << "'"
                                << ", writeOver = " << writeOver
                                << ": returned from CheckFilePath() with irtrn = " << irtrn
                                << std::endl;
      }
      queso_require_greater_equal_msg(irtrn, 0, "unable to verify output path");

      if (writeOver) {
        //////////////////////////////////////////////////////////////
        // Write over an eventual pre-existing file
        //////////////////////////////////////////////////////////////
        if ((fileType == UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT) ||
            (fileType == UQ_FILE_EXTENSION_FOR_TXT_FORMAT)) {
          filePtrSet.ofsVar = new std::ofstream((baseFileName+"_sub"+this->subIdString()+"."+fileType).c_str(),
                                                std::ofstream::out | std::ofstream::trunc);
        }
        else if (fileType == UQ_FILE_EXTENSION_FOR_HDF_FORMAT) {
          queso_error_msg("hdf file type not supported yet");
        }
        else {
          queso_error_msg("invalid file type");
        }
        if ((m_subDisplayFile) && (this->displayVerbosity() > 10)) { // output debug
          *this->subDisplayFile() << "In BaseEnvironment::openOutputFile()"
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
        if ((fileType == UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT) ||
            (fileType == UQ_FILE_EXTENSION_FOR_TXT_FORMAT)) {
#if 0
          filePtrSet.ofsVar = new std::ofstream((baseFileName+"_sub"+this->subIdString()+"."+fileType).c_str(),
                                                std::ofstream::in);
          std::cout << "filePtrSet.ofsVar(1) = " << filePtrSet.ofsVar << std::endl;
          if (filePtrSet.ofsVar) std::cout << "filePtrSet.ofsVar(1)->is_open() = " << filePtrSet.ofsVar->is_open() << std::endl;
          if (filePtrSet.ofsVar) delete filePtrSet.ofsVar;
#endif
          // 'm' and Ranger nodes behave differently on ofstream constructor... prudenci 2010/03/05
          filePtrSet.ofsVar = new std::ofstream((baseFileName+"_sub"+this->subIdString()+"."+fileType).c_str(),
                                                std::ofstream::out /*| std::ofstream::in*/ | std::ofstream::app);
        }
#ifdef QUESO_HAS_HDF5
        else if (fileType == UQ_FILE_EXTENSION_FOR_HDF_FORMAT) {
          std::string fullFileName =
            baseFileName+"_sub"+this->subIdString()+"."+fileType;

          // Use H5F_ACC_EXCL because not overwriting, so fail on existing file
          filePtrSet.h5Var = H5Fcreate(fullFileName.c_str(),
                                       H5F_ACC_EXCL,
                                       H5P_DEFAULT,
                                       H5P_DEFAULT);

          queso_require_greater_equal_msg(
              filePtrSet.h5Var, 0,
              "error opening file `" << fullFileName << "`");
        }
#endif
        else {
          queso_error_msg("invalid file type");
        }
        if ((m_subDisplayFile) && (this->displayVerbosity() > 10)) { // output debug
          *this->subDisplayFile() << "In BaseEnvironment::openOutputFile()"
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
              *this->subDisplayFile() << "In BaseEnvironment::openOutputFile()"
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

      // Check the file actually opened
      if ((fileType == UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT) ||
          (fileType == UQ_FILE_EXTENSION_FOR_TXT_FORMAT)) {
        queso_require_msg(
            (filePtrSet.ofsVar && filePtrSet.ofsVar->is_open()),
            "failed to open output file");
      }
    }
    else {
      returnValue = false;
    }
    //this->subComm().Barrier(); // prudenci-2011-01-17
  }

  return returnValue;
}
//-------------------------------------------------------
bool
BaseEnvironment::openUnifiedOutputFile(
  const std::string&        baseFileName,
  const std::string&        inputFileType,
        bool                writeOver,
        FilePtrSetStruct& filePtrSet) const
{
  std::string fileType(inputFileType);
#ifdef QUESO_HAS_HDF5
  // Do nothing
#else
  if (fileType == UQ_FILE_EXTENSION_FOR_HDF_FORMAT) {
    if (m_subDisplayFile) {
      *this->subDisplayFile() << "WARNING in BaseEnvironment::openUnifiedOutputFile()"
                              << ": file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                              << "' has been requested, but this QUESO library has not been built with 'hdf5'"
                              << ". Code will therefore process the file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                              << "' instead..."
                              << std::endl;
    }
    if (this->subRank() == 0) {
      std::cerr << "WARNING in BaseEnvironment::openUnifiedOutputFile()"
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
      *this->subDisplayFile() << "In BaseEnvironment::openUnifiedOutputFile()"
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
      *this->subDisplayFile() << "In BaseEnvironment::openUnifiedOutputFile()"
                              << ": opening unified output file with base name '" << baseFileName << "." << fileType
                              << "'"
                              << ", writeOver = " << writeOver
                              << std::endl;
    }


    //if ((this->subRank   () == 0) &&
    //    (this->inter0Rank() == 0)) {
#if 0
      std::cout << "In BaseEnvironment::openUnifiedOutputFile()"
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
      int irtrn = CheckFilePath((baseFileName+"."+fileType).c_str());
      queso_require_greater_equal_msg(irtrn, 0, "unable to verify output path");

      if (writeOver) {
        ////////////////////////////////////////////////////////////////
        // Write over an eventual pre-existing file
        ////////////////////////////////////////////////////////////////
        if ((fileType == UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT) ||
            (fileType == UQ_FILE_EXTENSION_FOR_TXT_FORMAT)) {
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
          queso_error_msg("invalid file type");
        }
        if ((m_subDisplayFile) && (this->displayVerbosity() > 10)) { // output debug
          *this->subDisplayFile() << "In BaseEnvironment::openUnifiedOutputFile()"
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
        // 'm' and Ranger nodes behave differently on ofstream constructor... prudenci 2010/03/05
        ////////////////////////////////////////////////////////////////
        if ((fileType == UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT) ||
            (fileType == UQ_FILE_EXTENSION_FOR_TXT_FORMAT)) {
          filePtrSet.ofsVar = new std::ofstream((baseFileName+"."+fileType).c_str(),
                                                std::ofstream::out /*| std::ofstream::in*/ | std::ofstream::app);
        }
#ifdef QUESO_HAS_HDF5
        else if (fileType == UQ_FILE_EXTENSION_FOR_HDF_FORMAT) {
          filePtrSet.h5Var = H5Fcreate((baseFileName+"."+fileType).c_str(), // TEMPORARY - FIX ME
                                       H5F_ACC_TRUNC,
                                       H5P_DEFAULT,
                                       H5P_DEFAULT);

          //                    m_worldRank,
          //                    "BaseEnvironment::openUnifiedOutputFile(), writeOver=false",
          //                    "hdf file type not supported yet");
        }
#endif
        else {
          queso_error_msg("invalid file type");
        }
        if ((m_subDisplayFile) && (this->displayVerbosity() > 10)) { // output debug
          *this->subDisplayFile() << "In BaseEnvironment::openUnifiedOutputFile()"
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
            *this->subDisplayFile() << "In BaseEnvironment::openUnifiedOutputFile()"
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
        std::cerr << "In BaseEnvironment::openUnifiedOutputFile()"
                  << ": failed to open unified output file with base name '" << baseFileName << "." << fileType
                  << "'"
                  << std::endl;
      }
      queso_require_msg((filePtrSet.ofsVar && filePtrSet.ofsVar->is_open()), "failed to open output file");
    //}
  }

  return returnValue;
}
//-------------------------------------------------------
bool
BaseEnvironment::openInputFile(
  const std::string&            baseFileName,
  const std::string&            inputFileType,
  const std::set<unsigned int>& allowedSubEnvIds,
        FilePtrSetStruct&     filePtrSet) const
{
  std::string fileType(inputFileType);
#ifdef QUESO_HAS_HDF5
  // Do nothing
#else
  if (fileType == UQ_FILE_EXTENSION_FOR_HDF_FORMAT) {
    if (m_subDisplayFile) {
      *this->subDisplayFile() << "WARNING in BaseEnvironment::openInputFile()"
                              << ": file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                              << "' has been requested, but this QUESO library has not been built with 'hdf5'"
                              << ". Code will therefore process the file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                              << "' instead..."
                              << std::endl;
    }
    if (this->subRank() == 0) {
      std::cerr << "WARNING in BaseEnvironment::openInputFile()"
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
      *this->subDisplayFile() << "In BaseEnvironment::openInputFile()"
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
      *this->subDisplayFile() << "In BaseEnvironment::openInputFile()"
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
      int irtrn = CheckFilePath((baseFileName+"."+fileType).c_str());
      queso_require_greater_equal_msg(irtrn, 0, "unable to verify input path");

      if ((fileType == UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT) ||
          (fileType == UQ_FILE_EXTENSION_FOR_TXT_FORMAT)) {
        filePtrSet.ifsVar = new std::ifstream((baseFileName+"."+fileType).c_str(),
                                              std::ofstream::in);
        if ((filePtrSet.ifsVar == NULL) || (filePtrSet.ifsVar->is_open() == false)) {
          std::cerr << "In BaseEnvironment::openInputFile()"
                    << ": failed to open input file with base name '" << baseFileName << "." << fileType
                    << "'"
                    << std::endl;
        }
        queso_require_msg(!((filePtrSet.ifsVar == NULL) || (filePtrSet.ifsVar->is_open() == false)), "file with fileName could not be found");
      }
#ifdef QUESO_HAS_HDF5
      else if (fileType == UQ_FILE_EXTENSION_FOR_HDF_FORMAT) {
        filePtrSet.h5Var = H5Fopen((baseFileName+"."+fileType).c_str(),
                                   H5F_ACC_RDONLY,
                                   H5P_DEFAULT);
      }
#endif
      else {
        queso_error_msg("invalid file type");
      }
    }
    else {
      returnValue = false;
    }
    //this->subComm().Barrier(); // prudenci-2011-01-17
  }

  return returnValue;
}
//-------------------------------------------------------
bool
BaseEnvironment::openUnifiedInputFile(
  const std::string&        baseFileName,
  const std::string&        inputFileType,
        FilePtrSetStruct& filePtrSet) const
{
  std::string fileType(inputFileType);
#ifdef QUESO_HAS_HDF5
  // Do nothing
#else
  if (fileType == UQ_FILE_EXTENSION_FOR_HDF_FORMAT) {
    if (m_subDisplayFile) {
      *this->subDisplayFile() << "WARNING in BaseEnvironment::openUnifiedInputFile()"
                              << ": file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                              << "' has been requested, but this QUESO library has not been built with 'hdf5'"
                              << ". Code will therefore process the file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                              << "' instead..."
                              << std::endl;
    }
    if (this->subRank() == 0) {
      std::cerr << "WARNING in BaseEnvironment::openUnifiedInputFile()"
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
      *this->subDisplayFile() << "In BaseEnvironment::openUnifiedInputFile()"
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
      *this->subDisplayFile() << "In BaseEnvironment::openUnifiedInputFile()"
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
      int irtrn = CheckFilePath((baseFileName+"."+fileType).c_str());
      queso_require_greater_equal_msg(irtrn, 0, "unable to verify input path");

      if ((fileType == UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT) ||
          (fileType == UQ_FILE_EXTENSION_FOR_TXT_FORMAT)) {
        filePtrSet.ifsVar = new std::ifstream((baseFileName+"."+fileType).c_str(),
                                              std::ofstream::in);
        if ((filePtrSet.ifsVar == NULL) || (filePtrSet.ifsVar->is_open() == false)) {
          std::cerr << "In BaseEnvironment::openUnifiedInputFile()"
                    << ": failed to open input file with base name '" << baseFileName << "." << fileType
                    << "'"
                    << std::endl;
        }
        queso_require_msg(!((filePtrSet.ifsVar == NULL) || (filePtrSet.ifsVar->is_open() == false)), "file with fileName could not be found");
      }
#ifdef QUESO_HAS_HDF5
      else if (fileType == UQ_FILE_EXTENSION_FOR_HDF_FORMAT) {
        filePtrSet.h5Var = H5Fopen((baseFileName+"."+fileType).c_str(),
                                   H5F_ACC_RDONLY,
                                   H5P_DEFAULT);
      }
#endif
      else {
        queso_error_msg("invalid file type");
      }
    }
    //else {
    //  returnValue = false;
    //}
    //this->subComm().Barrier();
  }

  return returnValue;
}
//-------------------------------------------------------
void
BaseEnvironment::closeFile(
  FilePtrSetStruct& filePtrSet,
  const std::string&  inputFileType) const
{
  std::string fileType(inputFileType);
#ifdef QUESO_HAS_HDF5
  // Do nothing
#else
  if (fileType == UQ_FILE_EXTENSION_FOR_HDF_FORMAT) {
    if (m_subDisplayFile) {
      *this->subDisplayFile() << "WARNING in BaseEnvironment::closeFile()"
                              << ": file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                              << "' has been requested, but this QUESO library has not been built with 'hdf5'"
                              << ". Code will therefore process the file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                              << "' instead..."
                              << std::endl;
    }
    if (this->subRank() == 0) {
      std::cerr << "WARNING in BaseEnvironment::closeFile()"
                << ": file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                << "' has been requested, but this QUESO library has not been built with 'hdf5'"
                << ". Code will therefore process the file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                << "' instead..."
                << std::endl;
    }
    fileType = UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT;
  }
#endif

  if ((fileType == UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT) ||
      (fileType == UQ_FILE_EXTENSION_FOR_TXT_FORMAT)) {
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
    queso_error_msg("invalid file type");
  }

  return;
}
//-------------------------------------------------------
void
BaseEnvironment::setExceptionalCircumstance(bool value) const
{
  m_exceptionalCircumstance = value;
  return;
}
//-------------------------------------------------------
bool
BaseEnvironment::exceptionalCircumstance() const
{
  return m_exceptionalCircumstance;
}

const GetPot &
BaseEnvironment::input() const
{
  return *m_input;
}


//*****************************************************
// Empty Environment
//*****************************************************
EmptyEnvironment::EmptyEnvironment()
  :
  BaseEnvironment("",NULL)
{
}
//-------------------------------------------------------
EmptyEnvironment::~EmptyEnvironment()
{
}
//-------------------------------------------------------
void
EmptyEnvironment::print(std::ostream& os) const
{
  os.flush(); // just to avoid icpc warnings
  return;
}

//*****************************************************
// Full Environment
//*****************************************************
#ifdef QUESO_HAS_MPI
FullEnvironment::FullEnvironment(
  RawType_MPI_Comm             inputComm,
  const char*                    passedOptionsInputFileName,
  const char*                    prefix,
  EnvOptionsValues* alternativeOptionsValues)
  :
  BaseEnvironment(passedOptionsInputFileName,alternativeOptionsValues)
{
  this->construct(inputComm, prefix);
}

FullEnvironment::FullEnvironment(
  RawType_MPI_Comm             inputComm,
  const std::string&           passedOptionsInputFileName,
  const std::string&           prefix,
  EnvOptionsValues* alternativeOptionsValues)
  :
  BaseEnvironment(passedOptionsInputFileName,alternativeOptionsValues)
{
  this->construct(inputComm, prefix.c_str());
}

void
FullEnvironment::construct (RawType_MPI_Comm inputComm,
                            const char *prefix)
{
#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "Entering FullEnv" << std::endl;
#endif

  ///////////////////////////////////////////////////////////////////////
  // Initialize "full" communicator -- Not necessarily MPI_COMM_WORLD
  ///////////////////////////////////////////////////////////////////////
  int mpiRC = MPI_Comm_rank(inputComm,&m_worldRank);
  queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, "failed to get world fullRank()");

  m_fullComm.reset(new MpiComm(*this,inputComm));

  m_fullRank     = m_fullComm->MyPID();
  m_fullCommSize = m_fullComm->NumProc();

  mpiRC = MPI_Comm_group(m_fullComm->Comm(), &m_fullGroup);
  queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, "failed MPI_Comm_group()");

  // saving old uncaught exception handler, invoking queso_terminate
  old_terminate_handler = std::set_terminate(queso_terminate_handler);

#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "In FullEnv, finished dealing with MPI initially" << std::endl;
#endif

  //////////////////////////////////////////////////
  // Read options
  //////////////////////////////////////////////////
  // If NULL, we create one
  if (m_optionsObj == NULL) {
    // If there's an input file, we grab the options from there.  Otherwise the
    // defaults are used
    if (m_optionsInputFileName != "") {
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
      m_allOptionsMap.reset(new boost::program_options::variables_map());
      m_allOptionsDesc.reset(new boost::program_options::options_description("Allowed options"));
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

      readOptionsInputFile();

      m_input->parse_input_file(m_optionsInputFileName);
    }

    m_optionsObj.reset(new EnvOptionsValues(this, prefix));
  }

  // If help option was supplied, print info
  if (m_optionsObj->m_help != "") {
    // We write to std::cout because subDisplayFile() isn't ready yet?
    std::cout << (*m_optionsObj) << std::endl;
  }

  queso_require_equal_to_msg(
      fullComm().NumProc() % m_optionsObj->m_numSubEnvironments, 0,
      "total number of processors in environment must be multiple of the specified number of subEnvironments");

#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "In FullEnv, finished scanning options" << std::endl;
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
      QUESO_version_print(std::cout);
    }

    if (m_fullRank == 0) {
      std::cout << "Beginning run at " << ctime(&m_timevalBegin.tv_sec)
                << std::endl;
    }
  }

  //////////////////////////////////////////////////
  // Deal with multiple subEnvironments: create the sub communicators, one for each subEnvironment
  //////////////////////////////////////////////////
  unsigned int numRanksPerSubEnvironment = m_fullCommSize/m_optionsObj->m_numSubEnvironments;

  m_subId = m_fullRank/numRanksPerSubEnvironment;
  char tmpSubId[16];
  sprintf(tmpSubId,"%u",m_subId);
  m_subIdString = tmpSubId;

  if (m_optionsObj->m_subDisplayAllowAll) {
    m_optionsObj->m_subDisplayAllowedSet.insert((unsigned int) m_subId);
  }

  std::vector<int> fullRanksOfMySubEnvironment(numRanksPerSubEnvironment,0);
  for (unsigned int i = 0; i < numRanksPerSubEnvironment; ++i) {
    fullRanksOfMySubEnvironment[i] = m_subId * numRanksPerSubEnvironment + i;
  }

  mpiRC = MPI_Group_incl(m_fullGroup, (int) numRanksPerSubEnvironment, &fullRanksOfMySubEnvironment[0], &m_subGroup);
  queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, "failed MPI_Group_incl() for a subEnvironment");

  RawType_MPI_Comm subRawComm;
  mpiRC = MPI_Comm_create(m_fullComm->Comm(), m_subGroup, &subRawComm);
  queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, "failed MPI_Comm_group() for a subEnvironment");
  m_subComm.reset(new MpiComm(*this,subRawComm));
  m_subRank     = m_subComm->MyPID();
  m_subCommSize = m_subComm->NumProc();

  //////////////////////////////////////////////////
  // Deal with multiple subEnvironments: create the self communicator
  //////////////////////////////////////////////////
  m_selfComm.reset(new MpiComm(*this,RawValue_MPI_COMM_SELF));

  //////////////////////////////////////////////////
  // Deal with multiple subEnvironments: create the inter0 communicator
  //////////////////////////////////////////////////
  std::vector<int> fullRanksOfInter0(m_optionsObj->m_numSubEnvironments,0);
  for (unsigned int i = 0; i < m_optionsObj->m_numSubEnvironments; ++i) {
    fullRanksOfInter0[i] = i * numRanksPerSubEnvironment;
  }
  mpiRC = MPI_Group_incl(m_fullGroup, (int) m_optionsObj->m_numSubEnvironments, &fullRanksOfInter0[0], &m_inter0Group);
  queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, "failed MPI_Group_incl() for inter0");
  RawType_MPI_Comm inter0RawComm;
  mpiRC = MPI_Comm_create(m_fullComm->Comm(), m_inter0Group, &inter0RawComm);
  queso_require_equal_to_msg(mpiRC, MPI_SUCCESS, "failed MPI_Comm_group() for inter0");
  if (m_fullRank%numRanksPerSubEnvironment == 0) {
    m_inter0Comm.reset(new MpiComm(*this,inter0RawComm));
    m_inter0Rank     = m_inter0Comm->MyPID();
    m_inter0CommSize = m_inter0Comm->NumProc();
  }

  if (m_optionsObj->m_subDisplayAllowAll) {
    // This situation has been already taken care of above
  }
  else if (m_optionsObj->m_subDisplayAllowInter0) {
    if (m_inter0Rank >= 0) {
      m_optionsObj->m_subDisplayAllowedSet.insert((unsigned int) m_subId);
    }
  }


  //////////////////////////////////////////////////
  // Open "screen" file
  //////////////////////////////////////////////////
  bool openFile = false;
  if ((m_subRank                                               == 0                                              ) &&
      (m_optionsObj->m_subDisplayFileName                 != UQ_ENV_FILENAME_FOR_NO_OUTPUT_FILE             ) &&
      (m_optionsObj->m_subDisplayAllowedSet.find(m_subId) != m_optionsObj->m_subDisplayAllowedSet.end())) {
    openFile = true;
  }

  if (openFile && m_worldRank == 0) {
    //////////////////////////////////////////////////////////////////
    // Verify parent directory exists (for cases when a user
    // specifies a relative path for the desired output file).
    //////////////////////////////////////////////////////////////////
    int irtrn = CheckFilePath((m_optionsObj->m_subDisplayFileName+"_sub"+m_subIdString+".txt").c_str());
    queso_require_greater_equal_msg(irtrn, 0, "unable to verify output path");
  }

  ////////////////////////////////////////////////////////////////////
  // Ensure that rank 0 has created path, if necessary, before other tasks use it
  ////////////////////////////////////////////////////////////////////
  m_fullComm->Barrier();

  if (openFile) {
    //////////////////////////////////////////////////////////////////
    // Always write over an eventual pre-existing file
    //////////////////////////////////////////////////////////////////
    m_subDisplayFile.reset(new std::ofstream((m_optionsObj->m_subDisplayFileName+"_sub"+m_subIdString+".txt").c_str(),
                                         std::ofstream::out | std::ofstream::trunc));
    queso_require_msg((m_subDisplayFile && m_subDisplayFile->is_open()), "failed to open sub screen file");

    QUESO_version_print(*m_subDisplayFile);

    *m_subDisplayFile << "Beginning run at " << ctime(&m_timevalBegin.tv_sec)
                      << std::endl;
  }

  //////////////////////////////////////////////////
  // Debug message related to subEnvironments
  //////////////////////////////////////////////////
  if (this->displayVerbosity() >= 2) {
    for (int i = 0; i < m_fullCommSize; ++i) {
      if (i == m_fullRank) {
        //std::cout << "In FullEnvironment::commonConstructor()"
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
  if (m_optionsObj->m_rngType == "gsl") {
    m_rngObject.reset(new RngGsl(m_optionsObj->m_seed,m_worldRank));
    m_basicPdfs.reset(new BasicPdfsGsl(m_worldRank));
  }
  else if (m_optionsObj->m_rngType == "boost") {
    m_rngObject.reset(new RngBoost(m_optionsObj->m_seed,m_worldRank));
    m_basicPdfs.reset(new BasicPdfsBoost(m_worldRank));
  }
  else if (m_optionsObj->m_rngType == "cxx11") {
#ifdef QUESO_HAVE_CXX11
    m_rngObject.reset(new RngCXX11(m_optionsObj->m_seed, m_worldRank));
    m_basicPdfs.reset(new BasicPdfsCXX11(m_worldRank));
#else
    queso_error_msg("C++11 RNGs requested, but QUESO wasn't compiled with C++11 support");
#endif
  }
  else {
    std::cerr << "In Environment::constructor()"
              << ": rngType = " << m_optionsObj->m_rngType
              << std::endl;
    queso_error_msg("the requested 'rngType' is not supported yet");
  }

  //////////////////////////////////////////////////
  // Leave commonConstructor()
  //////////////////////////////////////////////////
  m_fullComm->Barrier();
  m_fullEnvIsReady = true;

  if ((m_subDisplayFile) && (this->displayVerbosity() >= 5)) {
    *m_subDisplayFile << "Done with initializations at FullEnvironment::commonConstructor()"
                      << std::endl;
  }

  return;
}
#endif  // QUESO_HAS_MPI

FullEnvironment::FullEnvironment(
  const char*                    passedOptionsInputFileName,
  const char*                    prefix,
  EnvOptionsValues* alternativeOptionsValues)
  :
  BaseEnvironment(passedOptionsInputFileName,alternativeOptionsValues)
{
  this->construct(prefix);
}

FullEnvironment::FullEnvironment(
  const std::string&             passedOptionsInputFileName,
  const std::string&             prefix,
  EnvOptionsValues* alternativeOptionsValues)
  :
  BaseEnvironment(passedOptionsInputFileName,alternativeOptionsValues)
{
  this->construct(prefix.c_str());
}

void
FullEnvironment::construct (const char *prefix)
{
#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "Entering FullEnv" << std::endl;
#endif

  m_worldRank = 0;

  m_fullComm.reset(new MpiComm(*this));
  m_fullRank     = 0;
  m_fullCommSize = 1;

#ifndef QUESO_HAS_MPI
  m_fullGroup = 0;
#endif

  // saving old uncaught exception handler, invoking queso_terminate
  old_terminate_handler = std::set_terminate(queso_terminate_handler);

#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "In FullEnv, finished dealing with MPI initially" << std::endl;
#endif

  //////////////////////////////////////////////////
  // Read options
  //////////////////////////////////////////////////
  // If NULL, we create one
  if (m_optionsObj == NULL) {
    // If there's an input file, we grab the options from there.  Otherwise the
    // defaults are used
    if (m_optionsInputFileName != "") {
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
      m_allOptionsMap.reset(new boost::program_options::variables_map());
      m_allOptionsDesc.reset(new boost::program_options::options_description("Allowed options"));
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

      readOptionsInputFile();

      m_input->parse_input_file(m_optionsInputFileName);
    }

    m_optionsObj.reset(new EnvOptionsValues(this, prefix));
  }

  // If help option was supplied, print info
  if (m_optionsObj->m_help != "") {
    // We write to std::cout because subDisplayFile() isn't ready yet?
    std::cout << (*m_optionsObj) << std::endl;
  }

  queso_require_equal_to_msg(
      fullComm().NumProc() % m_optionsObj->m_numSubEnvironments, 0,
      "total number of processors in environment must be multiple of the specified number of subEnvironments");

#ifdef QUESO_MEMORY_DEBUGGING
  std::cout << "In FullEnv, finished scanning options" << std::endl;
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

    QUESO_version_print(std::cout);

    std::cout << "Beginning run at " << ctime(&m_timevalBegin.tv_sec)
              << std::endl;
  }

  m_subId = 0;
  char tmpSubId[16];
  sprintf(tmpSubId,"%u",m_subId);
  m_subIdString = tmpSubId;

  if (m_optionsObj->m_subDisplayAllowAll) {
    m_optionsObj->m_subDisplayAllowedSet.insert((unsigned int) m_subId);
  }

  int fullRanksOfMySubEnvironment = 1;

#ifndef QUESO_HAS_MPI
  m_subGroup = 0;
#endif

  m_subComm.reset(new MpiComm(*this));
  m_subRank     = 0;
  m_subCommSize = 1;

  //////////////////////////////////////////////////
  // Deal with multiple subEnvironments: create the self communicator
  //////////////////////////////////////////////////
  m_selfComm.reset(new MpiComm(*this));

  //////////////////////////////////////////////////
  // Deal with multiple subEnvironments: create the inter0 communicator
  //////////////////////////////////////////////////
  int fullRanksOfInter0 = 0;
#ifndef QUESO_HAS_MPI
  m_inter0Group = 0;
#endif
  m_inter0Comm.reset(new MpiComm(*this));
  m_inter0Rank     = 0;
  m_inter0CommSize = 1;

  if (m_optionsObj->m_subDisplayAllowAll) {
    // This situation has been already taken care of above
  }
  else if (m_optionsObj->m_subDisplayAllowInter0) {
    if (m_inter0Rank >= 0) {
      m_optionsObj->m_subDisplayAllowedSet.insert((unsigned int) m_subId);
    }
  }


  //////////////////////////////////////////////////
  // Open "screen" file
  //////////////////////////////////////////////////
  bool openFile = false;
  if ((m_subRank                                               == 0                                              ) &&
      (m_optionsObj->m_subDisplayFileName                 != UQ_ENV_FILENAME_FOR_NO_OUTPUT_FILE             ) &&
      (m_optionsObj->m_subDisplayAllowedSet.find(m_subId) != m_optionsObj->m_subDisplayAllowedSet.end())) {
    openFile = true;
  }

  if (openFile && m_worldRank == 0) {
    //////////////////////////////////////////////////////////////////
    // Verify parent directory exists (for cases when a user
    // specifies a relative path for the desired output file).
    //////////////////////////////////////////////////////////////////
    int irtrn = CheckFilePath((m_optionsObj->m_subDisplayFileName+"_sub"+m_subIdString+".txt").c_str());
    queso_require_greater_equal_msg(irtrn, 0, "unable to verify output path");
  }

  ////////////////////////////////////////////////////////////////////
  // Ensure that rank 0 has created path, if necessary, before other tasks use it
  ////////////////////////////////////////////////////////////////////
  m_fullComm->Barrier();

  if (openFile) {
    //////////////////////////////////////////////////////////////////
    // Always write over an eventual pre-existing file
    //////////////////////////////////////////////////////////////////
    m_subDisplayFile.reset(new std::ofstream((m_optionsObj->m_subDisplayFileName+"_sub"+m_subIdString+".txt").c_str(),
                                         std::ofstream::out | std::ofstream::trunc));
    queso_require_msg((m_subDisplayFile && m_subDisplayFile->is_open()), "failed to open sub screen file");

    QUESO_version_print(*m_subDisplayFile);

    *m_subDisplayFile << "Beginning run at " << ctime(&m_timevalBegin.tv_sec)
                      << std::endl;
  }

  //////////////////////////////////////////////////
  // Debug message related to subEnvironments
  //////////////////////////////////////////////////
  if (this->displayVerbosity() >= 2) {
    for (int i = 0; i < m_fullCommSize; ++i) {
      if (i == m_fullRank) {
        //std::cout << "In FullEnvironment::commonConstructor()"
        std::cout << "MPI node of worldRank "             << m_worldRank
                  << " has fullRank "                     << m_fullRank
                  << ", belongs to subEnvironment of id " << m_subId
                  << ", and has subRank "                 << m_subRank
                  << std::endl;

        std::cout << "MPI node of worldRank " << m_worldRank
                  << " belongs to sub communicator with full ranks";
          std::cout << " " << fullRanksOfMySubEnvironment;
  std::cout << "\n";

        if (m_inter0Comm) {
          std::cout << "MPI node of worldRank " << m_worldRank
                    << " also belongs to inter0 communicator with full ranks";
            std::cout << " " << fullRanksOfInter0;
          std::cout << ", and has inter0Rank " << m_inter0Rank;
        }
  std::cout << "\n";

  std::cout << std::endl;
      }
      m_fullComm->Barrier();
    }
  }

  //////////////////////////////////////////////////
  // Deal with seed
  //////////////////////////////////////////////////
  if (m_optionsObj->m_rngType == "gsl") {
    m_rngObject.reset(new RngGsl(m_optionsObj->m_seed,m_worldRank));
    m_basicPdfs.reset(new BasicPdfsGsl(m_worldRank));
  }
  else if (m_optionsObj->m_rngType == "boost") {
    m_rngObject.reset(new RngBoost(m_optionsObj->m_seed,m_worldRank));
    m_basicPdfs.reset(new BasicPdfsBoost(m_worldRank));
  }
  else if (m_optionsObj->m_rngType == "cxx11") {
#ifdef QUESO_HAVE_CXX11
    m_rngObject.reset(new RngCXX11(m_optionsObj->m_seed, m_worldRank));
    m_basicPdfs.reset(new BasicPdfsCXX11(m_worldRank));
#else
    queso_error_msg("C++11 RNGs requested, but QUESO wasn't compiled with C++11 support");
#endif
  }
  else {
    std::cerr << "In Environment::constructor()"
              << ": rngType = " << m_optionsObj->m_rngType
              << std::endl;
    queso_error_msg("the requested 'rngType' is not supported yet");
  }

  //////////////////////////////////////////////////
  // Leave commonConstructor()
  //////////////////////////////////////////////////
  m_fullComm->Barrier();
  m_fullEnvIsReady = true;

  if ((m_subDisplayFile) && (this->displayVerbosity() >= 5)) {
    *m_subDisplayFile << "Done with initializations at FullEnvironment::commonConstructor()"
                      << std::endl;
  }

  return;
}

//-------------------------------------------------------
FullEnvironment::~FullEnvironment()
{
}
//-------------------------------------------------------
void
FullEnvironment::print(std::ostream& os) const
{
  os.flush(); // just to avoid icpc warnings
  return;
}

void queso_terminate_handler()
{
#ifdef QUESO_HAS_MPI
  int mpi_initialized;
  MPI_Initialized (&mpi_initialized);

  if (mpi_initialized)
    {
      //MPI_Abort(m_fullComm->Comm(), 1);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  else
    {
      // The system terminate_handler may do useful things like printing
      // uncaught exception information, or the user may have created
      // their own terminate handler that we want to call.
      old_terminate_handler();
    }
#else
  old_terminate_handler();
#endif
  exit(1);
}


//-------------------------------------------------------
void
FullEnvironment::readOptionsInputFile()
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
                                   << "\n  '<eventual mpi commands and options> <Application> <InputFile>'"
                                   << "\nin the command line."
                                   << "\n"
                                   << std::endl;
    queso_error();
  }

  return;
}

}  // End namespace QUESO
