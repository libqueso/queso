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

#ifndef __UQ_ENVIRONMENT_H__
#define __UQ_ENVIRONMENT_H__

#include <uqDefines.h>
class uqEnvironmentOptionsClass;

#undef UQ_USES_COMMAND_LINE_OPTIONS

#include <uqMpiComm.h>
#ifdef QUESO_HAS_HDF5
#include <hdf5.h>
#endif
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <iostream>
#include <fstream>

#include <uqRngBase.h>

struct uqFilePtrSetStruct {
  uqFilePtrSetStruct();
 ~uqFilePtrSetStruct();

  std::ofstream* ofsVar;
  std::ifstream* ifsVar;
#ifdef QUESO_HAS_HDF5
  hid_t  h5Var;
#endif
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
/*! \file uqEnvironment.h
    \brief Class to set up a QUESO environment.
*/


/*!  \class uqBaseEnvironmentClass
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
/*! This class is virtual. It is inherited by 'uqEmptyEnvironmentClass' and 'uqFullEnvironmentClass'.
    The QUESO environment class is instantiated at the application level, right after 'MPI_Init(&argc,&argv)'. 
    The QUESO environment is required by reference by many constructors in the QUESO library, 
    and is available by reference from many classes as well.
-------------------------------------------------------------*/
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
<item> full = communicator passed to the constructor of uqBaseEnvironmentClass, of size N and usually equal to the world communicator;
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
 * as forming a "unified" sequence in a distributed way. Indeed, the virtual class 'uqVectorSequenceClass' 
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


class uqBaseEnvironmentClass {
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  uqBaseEnvironmentClass(const char* passedOptionsInputFileName, const uqEnvOptionsValuesClass* alternativeOptionsValues);
  
  //! Copy constructor. It should not be used be the user.
  uqBaseEnvironmentClass(const uqBaseEnvironmentClass& obj);
  
  //! Destructor
  /*! It deallocates memory and does other cleanup for the class object and its class members when 
   * the object is destroyed. It displays the total run time of the combo QUESO + application using
     the function gettimeofday() from a struct timeval (as specified in <sys/time.h>). */
  virtual ~uqBaseEnvironmentClass();
  //@}
  
  //! @name Set methods
  //@{
  //! Assignment operator. It should not be used be the user.
  uqBaseEnvironmentClass& operator= (const uqBaseEnvironmentClass& rhs);
  //@}

  //! @name Environment, Communicator and Options Input File methods
  //@{
  //! Returns whether the full environment class is ready (constructor has successfully been called).
  bool    fullEnvIsReady() const;
  
  //! Returns the process world rank.
  int     worldRank     () const;

  //! Returns the process full rank.
  int     fullRank      () const;
  
  //! Access function for uqMpiComm full communicator. 
  const uqMpiCommClass&   fullComm      () const; 

  //! Access function for sub-group. 
  uqRawType_MPI_Group     subGroup      () const;
  
  //! Access function for sub-rank. 
  int     subRank       () const;
  
  //! Access function for uqMpiComm sub communicator. 
  const uqMpiCommClass&   subComm       () const; 

  //! Access function for uqMpiComm self-communicator. 
  const uqMpiCommClass&   selfComm      () const; 

  //! Returns the process inter0 rank.
  int     inter0Rank    () const;
  
  //! Access function for uqMpiComm inter0-communicator. 
  const uqMpiCommClass&   inter0Comm    () const;

  //! Access function for m_subDisplayFile (displays file on stream).
  std::ofstream*  subDisplayFile() const;
  
  //! Access function for m_subDisplayFileName (displays filename on stream).
  std::string     subDisplayFileName    () const;

  //! Access function to the number of sub-environments.
  unsigned int    numSubEnvironments    () const;
  
  //! Access function to the number of each sub-enviroment Id: m_subId.
  unsigned int    subId () const;
  
  //! Access to the attribute m_subIdString; which stores the string for the sub-enviroment, and it will be used, for instance,    to create the output files for each sub-enviroment.
  const std::string&      subIdString   () const;
  
  //TODO Not implemented?
  void    checkTheParallelEnvironment   () const;

  //! Access to the attribute m_optionsInputFileName, which stores the  name of the input file passed by the user to QUESO.
  std::string     optionsInputFileName  () const;
  
  void    setOptionsInputFileAccessState(bool newState) const; // Yes, 'const'
	    

#ifdef UQ_USES_COMMAND_LINE_OPTIONS
  const po::options_description& allOptionsDesc () const;
#endif
  
  //! Access function to private attribute m_allOptionsMap. It is an instance of po::variables_map(), which
  //! allows concrete variables to map which store variables in real map.
  po::variables_map&      allOptionsMap () const;
  
  
  //! This method scans the input file provided by the user to QUESO.
  /*! It checks if no input file is passed and updates the private attribute m_allOptionsDesc, which
   * keeps all the options.*/
  void    scanInputFileForMyOptions     (const po::options_description& optionsDesc) const;
  
  //! Access function to private attribute m_displayVerbosity. It manages how much information will be
  //! release during the use of the QUESO library.
  unsigned int    displayVerbosity      () const;
  
  //! Access function to private attribute m_syncVerbosity.
  unsigned int    syncVerbosity () const;
  
  //! Access function to private attribute m_checkingLevel.
  unsigned int    checkingLevel () const;
  
  //! Access to the RNG object.
  const uqRngBaseClass*   rngObject     () const;
  
  //! Reset RNG seed.
  void    resetSeed     (int newSeedOption);
  
  //! Access to the RNG seed.
  int     seed  () const;
  
  //! Access to the platform name.
  std::string     platformName  () const;
  
   //! Access function to private attribute m_identifyingString: identifying string.
  std::string     identifyingString     () const;
  
  //! Reset private attribute m_identifyingString wit the value \c newString.
  void    resetIdentifyingString(const std::string& newString) const; // Yes, const
  
  //! //TODO Not implemented? Whether or not there is an option input file.
  bool    isThereInputFile      () const;
  
  //! Used to save the time when the combo `QUESO+user's application' started to run.
  struct timeval  timevalBegin  () const;
  //@}
  
  //! @name I/O methods
  //@{
    
  //! Opens an output file for each sub-environment that was chosen to send data to the file.  
  bool    openOutputFile(const std::string& fileName, const std::string& fileType, 
			 const std::set<unsigned int>& allowedSubEnvIds, bool writeOver,
			 uqFilePtrSetStruct& filePtrSet) const;

  //! Opens a unified output file, that will contain data from all sub-environments.
  bool    openUnifiedOutputFile (const std::string& fileName, const std::string& fileType,
				 bool writeOver, uqFilePtrSetStruct& filePtrSet) const;
				 
  //! Opens an input file.
  bool    openInputFile (const std::string& fileName, const std::string& fileType,
			 const std::set<unsigned int>& allowedSubEnvIds,
			 uqFilePtrSetStruct& filePtrSet) const;
			 
  //! Opens the unified input file.			 
  bool    openUnifiedInputFile  (const std::string& fileName, const std::string& fileType,
				 uqFilePtrSetStruct& filePtrSet) const;
				 
  //! Closes the file.				 
  void    closeFile     (uqFilePtrSetStruct& filePtrSet, const std::string& fileType) const; 
  
  //! Set an exceptional circunstance.
  void    setExceptionalCircunstance    (bool value) const;
  
    //! Decides whether there is an exceptional circunstance.
  bool    exceptionalCircunstance       () const;
  

  virtual void    print (std::ostream& os) const = 0;

  //@}
protected:
  bool       		m_fullEnvIsReady;
  int 	     		m_worldRank;

  uqMpiCommClass*    	m_fullComm;
  int 			m_fullRank;
  int 			m_fullCommSize;
  uqRawType_MPI_Group 	m_fullGroup;

  std::string		m_optionsInputFileName;
  mutable bool       	m_optionsInputFileAccessState; // Yes, 'mutable'
  po::options_description*   m_allOptionsDesc;
  po::variables_map* 	m_allOptionsMap;

  unsigned int          m_subId;
  std::string 		m_subIdString;
  uqRawType_MPI_Group 	m_subGroup;
  uqMpiCommClass*       m_subComm;
  int			m_subRank;
  int			m_subCommSize;

  uqMpiCommClass*    	m_selfComm;

  uqRawType_MPI_Group	m_inter0Group;
  uqMpiCommClass*    	m_inter0Comm;
  int			m_inter0Rank;
  int			m_inter0CommSize;

  mutable std::ofstream*     m_subDisplayFile;
  uqRngBaseClass*    	m_rngObject;
  struct timeval     	m_timevalBegin;
  mutable bool       	m_exceptionalCircunstance;

  uqEnvOptionsValuesClass    m_alternativeOptionsValues;
  uqEnvironmentOptionsClass* m_optionsObj;
};

//*****************************************************
// Empty Environment
//*****************************************************
/*!  \class uqEmptyEnvironmentClass
 *  \brief This class sets up the environment underlying the use of the QUESO library by an executable.
 */
class uqEmptyEnvironmentClass : public uqBaseEnvironmentClass {
public:
      //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor. Does nothing.
  /*! It initialized uqBaseEnvironmentClass with no input file and a NULL pointer for the alternativeOptionsValues.*/ 
  uqEmptyEnvironmentClass();
  
  //! Destructor
 ~uqEmptyEnvironmentClass();
  //@}
 
void print(std::ostream& os) const;
};

//*****************************************************
// Full Environment
//*****************************************************
/*!  \class uqFullEnvironmentClass
 *  \brief This class sets up the full environment underlying the use of the QUESO library by an executable.
 * 
 * This is the class that is actually used during a QUESO+application run.
 */

class uqFullEnvironmentClass : public uqBaseEnvironmentClass {
public:
    //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  /*! It initializes the full communicator, reads the options, deals with multiple subenvironments, 
   * e.g. dealing with sub/self/inter0-communicators, handles path for output files. */
  uqFullEnvironmentClass(uqRawType_MPI_Comm inputComm, const char* passedOptionsInputFileName, const char* prefix, const uqEnvOptionsValuesClass* alternativeOptionsValues);
 
  //! Destructor
 ~uqFullEnvironmentClass();
  //@}
 
  //! @name I/O methods
  //@{ 
  //! Sends the environment options to the stream.
  void	print       (std::ostream& os) const;
  //@}

private:
  //! Checks the options input file and reads the options.
void	readOptionsInputFile();
};

std::ostream& operator<<(std::ostream& os, const uqBaseEnvironmentClass& obj);

#endif // __UQ_ENVIRONMENT_H__
