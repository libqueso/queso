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

#ifndef UQ_ENVIRONMENT_H
#define UQ_ENVIRONMENT_H

#include <queso/Defines.h>
#undef UQ_USES_COMMAND_LINE_OPTIONS

#include <queso/MpiComm.h>
#include <queso/ScopedPtr.h>

#ifdef QUESO_HAS_HDF5
#include <hdf5.h>
#endif
#include <iostream>
#include <fstream>


#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
// Forward declarations
namespace boost {
  namespace program_options {
    class options_description;
    class variables_map;
    }
}
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

namespace QUESO {

// Forward declarations
class GetPot;
class EnvironmentOptions;
class EnvOptionsValues;
class BasicPdfsBase;
class RngBase;


  /*! queso_terminate_handler
   *  \brief Function for unhandled exceptions in Queso
   *
   *  This function deals with unhandled exceptions encountered in Queso.
   *  It provides a call to MPI_abort using the global communicator.
   */
  void queso_terminate_handler();

/*! \struct FilePtrSetStruct
 *  \brief Struct for handling data input and output from files.
 *
 *  This struct deals with data input and output from files.
 *  It encapsulates the input/output stream class std:: fstream.
 */

//!
struct FilePtrSetStruct {

  //! Struct constructor
  FilePtrSetStruct();

  //! Destructor
  ~FilePtrSetStruct();

  //! Provides a stream interface to write data to files.
  std::ofstream* ofsVar;

  //! Provides a stream interface to read data from files.
  std::ifstream* ifsVar;
#ifdef QUESO_HAS_HDF5
  hid_t  h5Var;
#endif
};

//------------------------------------------------------------------------
// Library versioning routines: we include them in a QUESO namespace
// here so that multiple classes can use them as required (and so we
// can have a standalone versioning binary which does not require a full
// QUESO MPI environment.
//------------------------------------------------------------------------

  void QUESO_version_print       (std::ostream &os);
  int  QUESO_get_numeric_version ();

//*****************************************************
// Base class
//*****************************************************
/*! \file Environment.h
 *  \brief Class to set up a QUESO environment.
 *  \class BaseEnvironment
 *  \brief This (virtual) class sets up the environment underlying the use of the QUESO library by an executable.
 */

/*! This class sets up the environment underlying the use of the QUESO library by an executable. It:
<list type=number>
<item> assigns rank numbers, other than the world rank, to nodes participating in a parallel job,
<item> provides communicators for generating a sequence of vectors in a distributed way,
<item> provides functionality to read options from the 'options input file' (whose name is passed
in the constructor of this environment class),
<item> opens output files for messages that would otherwise be written to the screen (one output
file per allowed rank is opened and allowed ranks can be specified through the 'options input file').
</list>
-------------------------------------------------------------*/
/*! This class is virtual. It is inherited by 'EmptyEnvironment' and 'FullEnvironment'.
    The QUESO environment class is instantiated at the application level, right after 'MPI_Init(&argc,&argv)'.
    The QUESO environment is required by reference by many constructors in the QUESO library,
    and is available by reference from many classes as well.
-------------------------------------------------------------*/
/*! Throughout QUESO, there are five classes whose constructors check options in the 'options input file':
<list type=number>
<item> BaseEnvironment
<item> StatisticalInverseProblem
<item> StatisticalForwardProblem
<item> MetropolisHastingsSG ('SG' stands for 'sequence generator')
<item> MonteCarloSG
</list>
*/
/*! These classes rely on 'options classes' to read their options from the input file.
    The options classes are, respectively:
<list type=number>
<item> EnvironmentOptions
<item> StatisticalInverseProblemOptions
<item> StatisticalForwardProblemOptions
<item> MetropolisHastingsSGOptions
<item> MonteCarloSGOptions
</list>
    The last two classes also rely on SequenceStatisticalOptions for reading the
    options specifying which statistics have to be computed on the sequences of vectors
    involved.
-------------------------------------------------------------*/

/*! The QUESO environment class manages five types of communicators. Let:
<list type=number>
<item> 'W >= 1' be the size of whole world communicator involved in a parallel run;
<item> 'N >= 1' be the size of the communicator passed to the QUESO environment constructor;
<item> 'S >= 1' be the number of statistical problems a QUESO environment will be handling
at the same time, in parallel.
</list>
    Usually 'W'='N', but such equality is not necessary.
    The number 'S' is equal to the QUESO environment option 'm_numSubEnvironments', and is equal to
    1 by default. The number 'N' must be a multiple of 'S', otherwise the QUESO class prints a fatal
    error message and MPI aborts. The five types of communicators that QUESO manages are referred to as:
<list type=number>
<item> world = MPI_WORLD_COMM, of size W;
<item> full = communicator passed to the constructor of BaseEnvironment, of size N and usually equal to the world communicator;
<item> sub = communicator of size N/S that contains the number of MPI nodes necessary to solve a statistical inverse problem or a statistical forward problem.
<item> self = MPI_SELF_COMM, of size 1;
<item> inter0 = communicator of size S formed by all MPI nodes that have 'sub' rank 0 in their respective 'sub' communicators.
</list>
    So, any given node has potentially five different ranks. Of course, if the user is solving just one statistical problem with just one MPI node, then all ranks are equal to zero.

-------------------------------------------------------------*/

/*! In the QUESO library terminology, one might refer to a QUESO "full" environment composed of
 * 'S' QUESO "sub" environments. Each sub environment is assigned a "sub" id varying from 0 (zero)
 * to S-1. Each sub environment is able to generate a statistical inverse problem and/or a statistical
 * forward problem. That is, each sub environment is able to handle a "sub" Markov chain (a sequence)
 * of vectors and/or a "sub" Monte Carlo sequence of output vectors. The "sub" sequences can be seen
 * as forming a "unified" sequence in a distributed way. Indeed, the virtual class 'VectorSequence'
 * provides "sub" and "unified" statistical operations.
 *
 *  A QUESO "sub" environment eventually prints messages to its own output file. In order for that to
 * happen, the requirements are:
<list type=number>
<item> option 'm_subDisplayFileName', a string, must be different than the default value ".";
<item> option 'm_subDisplayAllowedSet', a set of sub ids, must contain the id of the sub environment
wanting to write a message to the output file;
<item> the previous requirement is automatically satisfied if the option 'm_subDisplayAllowAll',
a boolean, is set to 1 (the default value is 0);
<item> the processor wanting to write a message to the output file must have sub rank 0 (zero).
</list>
   If all requirements are satisfied, then QUESO will generate a file with name '\<m_subDisplayFileName\>_sub\<sub id\>.txt'.
    For instance, if 'm_subDisplayFileName' is 'pROblem_775_' then a node of sub rank 0 in sub environment 17
    will write a message to the file 'pROblem_775_sub17.txt'.
*/


class BaseEnvironment {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  BaseEnvironment(const char* passedOptionsInputFileName, EnvOptionsValues* alternativeOptionsValues);

  BaseEnvironment(const std::string & passedOptionsInputFileName,
                  EnvOptionsValues* alternativeOptionsValues);

  //! Destructor
  /*! It deallocates memory and does other cleanup for the class object and its class members when
   * the object is destroyed. It displays the total run time of the combo QUESO + application using
     the function gettimeofday() from a struct timeval (as specified in <sys/time.h>). */
  virtual ~BaseEnvironment();
  //@}

  //! @name Environment, Communicator and Options Input File methods
  //@{
  //! Returns whether the full environment class is ready (constructor has successfully been called).
  bool    fullEnvIsReady() const;

  //! Returns the same thing as fullRank()
  /*!
   * This is the same thing as fullRank(), since QUESO's 'world' communicator
   * is not MPI_COMM_WORLD, but the communicator that user passed to it when
   * creating the environment.
   */
  int     worldRank     () const;

  //! Returns the rank of the MPI process in QUESO's full communicator
  /*!
   * Returns the rank of the MPI process in the communicator returned by
   * fullComm().
   *
   * See fullComm() for what the full communicator is.
   */
  int     fullRank      () const;

  //! Access function for the communicator that was passed to QUESO's environment
  /*!
   * The 'full' communicator is the MPI communicator that the user passed when
   * creating the QUESO FullEnvironment.  This is usually MPI_COMM_WORLD, but
   * the user is permitted to pass any MPI communicator smaller than
   * MPI_COMM_WORLD.
   */
  const MpiComm&   fullComm      () const;

  //! Access function for sub-group.
  RawType_MPI_Group     subGroup      () const;

  //! Returns the rank of the MPI process in the sub-communicator subComm()
  /*!
   * Example, if the calling MPI process has fullRank() equal to 3, the size of
   * fullComm() is 6, and the user asked for two sub-environments, then this
   * method will return 0. Here's why.
   *
   * fullComm() has MPI processes with these ranks:
   * 0 1 2 3 4 5
   *
   * QUESO divides the first three (ranks 0, 1, 2) of these into a
   * sub-communicator for sub-environment 0.  Inside the sub-communicator their
   * ranks are 0, 1, 2, respectively.
   *
   * QUESO divides the second three (ranks 3, 4, 5) of these into a
   * sub-communicator for sub-environment 1.  Inside the sub-communicator their
   * ranks are 0, 1, 2, respectively.
   *
   * It should be clear, now, that if fullRank() is 3 then subRank() is 0.
   */
  int     subRank       () const;

  //! Access function for each sub-environment's communicator.
  /*!
   * Let's say QUESO was passed a fullComm() communicator of size N.  The ranks
   * of each process in this communicator are:
   *
   * 0 1 2 ... N-2 N-1
   *
   * If the user asks for M sub-environments (chains) then, assuming M divides
   * N, QUESO partitions the processes in the fullComm() communicator
   * into M sub-communicators like so:
   *
   * Sub-environment 0 contains processes with fullRank()
   * 0 1 ... M-1
   *
   * Sub-environment 1 contains processes with fullRank()
   * M M+1 ... 2M-1
   *
   * et cetera
   *
   * Sub-environment M-1 contains processes with fullRank()
   * N-M N-M+1 ... N-1
   *
   * subComm() returns the sub-communicator corresponding to the
   * sub-environment the calling MPI process belongs to.  For example, if I am
   * an MPI process calling this function and I live in sub-environment \c k,
   * then this method returns the sub-communicator for sub-environment k.
   */
  const MpiComm&   subComm       () const;

  //! Access function for MpiComm self-communicator.
  /*!
   * This communicator is exactly MPI_COMM_SELF.
   */
  const MpiComm&   selfComm      () const;

  //! Returns the process inter0 rank.
  int     inter0Rank    () const;

  //! Access function for MpiComm communicator for processes with subRank() 0
  /*
   * This communicator contains all the processes that have subRank() equal to
   * 0.
   *
   * Their corresponding fullRank() values will be 0, M, 2M, ..., N-M,
   * where M is the number of sub-environments the user asked for and N is the
   * size of fullComm().
   */
  const MpiComm&   inter0Comm    () const;

  //! Access function for m_subDisplayFile (displays file on stream).
  std::ofstream*  subDisplayFile() const;

  //! Access function for m_subDisplayFileName (displays filename on stream).
  std::string     subDisplayFileName    () const;

  //! Access function to the number of sub-environments.
  unsigned int    numSubEnvironments    () const;

  //! Access function to the number of each sub-environment Id: m_subId.
  unsigned int    subId () const;

  //! Access to the attribute m_subIdString; which stores the string for the sub-environment, and it will be used, for instance,    to create the output files for each sub-environment.
  const std::string&      subIdString   () const;

  //TODO Not implemented?
  void    checkTheParallelEnvironment   () const;

  //! Access to the attribute m_optionsInputFileName, which stores the  name of the input file passed by the user to QUESO.
  std::string     optionsInputFileName  () const;

  void    setOptionsInputFileAccessState(bool newState) const; // Yes, 'const'


#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
#ifdef UQ_USES_COMMAND_LINE_OPTIONS
  const boost::program_options::options_description& allOptionsDesc () const;
#endif
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  //! Access function to private attribute m_allOptionsMap. It is an instance of boost::program_options::variables_map(), which
  //! allows concrete variables to map which store variables in real map.
  boost::program_options::variables_map&      allOptionsMap () const;
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS


#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  //! This method scans the input file provided by the user to QUESO.
  /*! It checks if no input file is passed and updates the private attribute m_allOptionsDesc, which
   * keeps all the options.*/
  void    scanInputFileForMyOptions(const boost::program_options::options_description& optionsDesc) const;
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

  //! Access function to private attribute m_displayVerbosity. It manages how much information will be
  //! release during the use of the QUESO library.
  unsigned int    displayVerbosity () const;

  //! Access function to private attribute m_syncVerbosity.
  unsigned int    syncVerbosity    () const;

  //! Access function to private attribute m_checkingLevel.
  unsigned int    checkingLevel    () const;

  //! Access to the RNG object.
  const RngBase* rngObject  () const;

  //! Reset RNG seed.
  void                  resetSeed  (int newSeedOption);

  //! Access to the RNG seed.
  int                   seed       () const;

  //! Access to Basic PDFs.
  const BasicPdfsBase* basicPdfs() const;

  //! Access to the platform name.
  std::string platformName      () const;

   //! Access function to private attribute m_identifyingString: identifying string.
  std::string identifyingString () const;

  //! Reset private attribute m_identifyingString with the value \c newString.
  void    resetIdentifyingString(const std::string& newString);

  //! //TODO Not implemented? it is called in examples/validationCycle/tests_old/results_5_25/uqTgaEx4.h.
  bool    isThereInputFile      () const;

  //! Used to save the time when the combo `QUESO+user's application' started to run.
  struct timeval  timevalBegin  () const;
  //@}

  //! @name I/O methods
  //@{

  //! Opens an output file for each sub-environment that was chosen to send data to the file.
  bool    openOutputFile(const std::string& fileName, const std::string& fileType,
                         const std::set<unsigned int>& allowedSubEnvIds, bool writeOver,
                         FilePtrSetStruct& filePtrSet) const;

  //! Opens a unified output file, that will contain data from all sub-environments.
  bool    openUnifiedOutputFile (const std::string& fileName, const std::string& fileType,
                                 bool writeOver, FilePtrSetStruct& filePtrSet) const;

  //! Opens an input file.
  bool    openInputFile (const std::string& fileName, const std::string& fileType,
                         const std::set<unsigned int>& allowedSubEnvIds,
                         FilePtrSetStruct& filePtrSet) const;

  //! Opens the unified input file.
  bool    openUnifiedInputFile  (const std::string& fileName, const std::string& fileType,
                                 FilePtrSetStruct& filePtrSet) const;

  //! Closes the file.
  void    closeFile     (FilePtrSetStruct& filePtrSet, const std::string& fileType) const;

  //! Set an exceptional circumstance.
  void    setExceptionalCircumstance    (bool value) const;

    //! Decides whether there is an exceptional circumstance.
  bool    exceptionalCircumstance       () const;

  //! The GetPot input file parser
  const GetPot & input() const;


  virtual void    print (std::ostream& os) const = 0;

  //@}
protected:
  bool m_fullEnvIsReady;
  int m_worldRank;

  ScopedPtr<MpiComm>::Type m_fullComm;
  int m_fullRank;
  int m_fullCommSize;
  RawType_MPI_Group m_fullGroup;

  std::string m_optionsInputFileName;
  mutable bool m_optionsInputFileAccessState; // Yes, 'mutable'
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  ScopedPtr<boost::program_options::options_description>::Type m_allOptionsDesc;
  ScopedPtr<boost::program_options::variables_map>::Type m_allOptionsMap;
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
  ScopedPtr<GetPot>::Type m_input;

  unsigned int m_subId;
  std::string m_subIdString;
  RawType_MPI_Group m_subGroup;
  ScopedPtr<MpiComm>::Type m_subComm;
  int m_subRank;
  int m_subCommSize;

  ScopedPtr<MpiComm>::Type m_selfComm;

  RawType_MPI_Group m_inter0Group;
  ScopedPtr<MpiComm>::Type m_inter0Comm;
  int m_inter0Rank;
  int m_inter0CommSize;

  mutable ScopedPtr<std::ofstream>::Type m_subDisplayFile;
  ScopedPtr<RngBase>::Type m_rngObject;
  ScopedPtr<BasicPdfsBase>::Type m_basicPdfs;
  struct timeval m_timevalBegin;
  mutable bool m_exceptionalCircumstance;

  ScopedPtr<EnvOptionsValues>::Type m_optionsObj;
};

//*****************************************************
// Empty Environment
//*****************************************************
/*!  \class EmptyEnvironment
 *  \brief This class sets up the environment underlying the use of the QUESO library by an executable.
 */
class EmptyEnvironment : public BaseEnvironment {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor. Does nothing.
  /*! It initialized BaseEnvironment with no input file and a NULL pointer for the alternativeOptionsValues.*/
  EmptyEnvironment();

  //! Destructor
 ~EmptyEnvironment();
  //@}

void print(std::ostream& os) const;
};

//*****************************************************
// Full Environment
//*****************************************************
/*!  \class FullEnvironment
 *  \brief This class sets up the full environment underlying the use of the QUESO library by an executable.
 *
 * This is the class that is actually used during a QUESO+application run.
 */

class FullEnvironment : public BaseEnvironment {
public:
    //! @name Constructor/Destructor methods
  //@{
  //! Parallel constructor.
  /*!
   * Initializes the full communicator, reads the options, deals with multiple
   * sub-environments, e.g. dealing with sub/self/inter0-communicators, handles
   * path for output files.
   */
#ifdef QUESO_HAS_MPI
  FullEnvironment(RawType_MPI_Comm inputComm, const char* passedOptionsInputFileName, const char* prefix, EnvOptionsValues* alternativeOptionsValues);

  FullEnvironment(RawType_MPI_Comm inputComm,
                  const std::string& passedOptionsInputFileName,
                  const std::string& prefix,
                  EnvOptionsValues* alternativeOptionsValues);
#endif

  //! Serial constructor.
  /*!
   * No communicator is passed. Output path handling is exactly as in the
   * parallel ctor.
   */
  FullEnvironment(const char* passedOptionsInputFileName, const char* prefix, EnvOptionsValues* alternativeOptionsValues);

  FullEnvironment(const std::string& passedOptionsInputFileName,
                  const std::string& prefix,
                  EnvOptionsValues* alternativeOptionsValues);

  //! Destructor
 ~FullEnvironment();
  //@}

  //! @name I/O methods
  //@{
  //! Sends the environment options to the stream.
  void        print       (std::ostream& os) const;
  //@}

private:
#ifdef QUESO_HAS_MPI
  //! Named constructor backend for multiple constructor overloads
  void        construct(RawType_MPI_Comm inputComm,
                  const char *prefix);
#endif

  //! Named constructor backend for multiple constructor overloads
  void        construct(const char *prefix);

  //! Checks the options input file and reads the options.
  void        readOptionsInputFile();
  //void        queso_terminate_handler();

};

std::ostream& operator<<(std::ostream& os, const BaseEnvironment& obj);

}  // End namespace QUESO

#endif // UQ_ENVIRONMENT_H
