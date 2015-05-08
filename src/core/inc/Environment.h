//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
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
#ifdef QUESO_HAS_HDF5
#include <hdf5.h>
#endif
#include <iostream>
#include <fstream>

#include <queso/RngBase.h>
#include <queso/BasicPdfsBase.h>

// Forward declarations
namespace boost {
  namespace program_options {
    class options_description;
    class variables_map;
    }
}

namespace QUESO {

// Forward declarations
class EnvironmentOptions;
class EnvOptionsValues;


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

  //! Returns the process world rank.
  int     worldRank     () const;

  //! Returns the process full rank.
  int     fullRank      () const;

  //! Access function for MpiComm full communicator.
  const MpiComm&   fullComm      () const;

  //! Access function for sub-group.
  RawType_MPI_Group     subGroup      () const;

  //! Access function for sub-rank.
  int     subRank       () const;

  //! Access function for MpiComm sub communicator.
  const MpiComm&   subComm       () const;

  //! Access function for MpiComm self-communicator.
  const MpiComm&   selfComm      () const;

  //! Returns the process inter0 rank.
  int     inter0Rank    () const;

  //! Access function for MpiComm inter0-communicator.
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


#ifdef UQ_USES_COMMAND_LINE_OPTIONS
  const boost::program_options::options_description& allOptionsDesc () const;
#endif

  //! Access function to private attribute m_allOptionsMap. It is an instance of boost::program_options::variables_map(), which
  //! allows concrete variables to map which store variables in real map.
  boost::program_options::variables_map&      allOptionsMap () const;


  //! This method scans the input file provided by the user to QUESO.
  /*! It checks if no input file is passed and updates the private attribute m_allOptionsDesc, which
   * keeps all the options.*/
  void    scanInputFileForMyOptions(const boost::program_options::options_description& optionsDesc) const;

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


  virtual void    print (std::ostream& os) const = 0;

  //@}
protected:
  bool       		     m_fullEnvIsReady;
  int 	     		     m_worldRank;

  MpiComm*    	     m_fullComm;
  int                        m_fullRank;
  int                        m_fullCommSize;
  RawType_MPI_Group        m_fullGroup;

  std::string		     m_optionsInputFileName;
  mutable bool       	     m_optionsInputFileAccessState; // Yes, 'mutable'
  boost::program_options::options_description*   m_allOptionsDesc;
  boost::program_options::variables_map* 	     m_allOptionsMap;

  unsigned int               m_subId;
  std::string 		     m_subIdString;
  RawType_MPI_Group        m_subGroup;
  MpiComm*            m_subComm;
  int			     m_subRank;
  int			     m_subCommSize;

  MpiComm*            m_selfComm;

  RawType_MPI_Group        m_inter0Group;
  MpiComm*            m_inter0Comm;
  int	                     m_inter0Rank;
  int                        m_inter0CommSize;

  mutable std::ofstream*     m_subDisplayFile;
  RngBase*    	     m_rngObject;
  BasicPdfsBase*      m_basicPdfs;
  struct timeval             m_timevalBegin;
  mutable bool       	     m_exceptionalCircumstance;

  EnvOptionsValues * m_optionsObj;
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
  //! Default constructor.
  /*! It initializes the full communicator, reads the options, deals with multiple sub-environments,
   * e.g. dealing with sub/self/inter0-communicators, handles path for output files. */
  FullEnvironment(RawType_MPI_Comm inputComm, const char* passedOptionsInputFileName, const char* prefix, EnvOptionsValues* alternativeOptionsValues);

  //! Destructor
 ~FullEnvironment();
  //@}

  //! @name I/O methods
  //@{
  //! Sends the environment options to the stream.
  void	print       (std::ostream& os) const;
  //@}

private:
  //! Checks the options input file and reads the options.
  void	readOptionsInputFile();
  //void	queso_terminate_handler();

};

std::ostream& operator<<(std::ostream& os, const BaseEnvironment& obj);

}  // End namespace QUESO

#endif // UQ_ENVIRONMENT_H
