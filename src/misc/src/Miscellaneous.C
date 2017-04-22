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

#include <cstring>
#include <queso/Defines.h>
#include <queso/Miscellaneous.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <libgen.h>
#include <sys/stat.h>
#include <cmath>

namespace QUESO {

void
MiscReadDoublesFromString(
  const std::string&         inputString,
        std::vector<double>& outputDoubles)
{
  //std::cout << "In MiscReadDoublesFromString()"
  //          << ": inputString = " << inputString
  //          << std::endl;
  outputDoubles.clear();
  bool aDoubleIsBeingRead = false;
  std::string::size_type positionOfFirstChar = 0;
  std::string::size_type numberOfChars = 0;
  for (std::string::size_type i = 0; i < inputString.size(); ++i) {
    queso_require_not_equal_to_msg(inputString[i], '\0', "character '\0' should not be found!");
    if (inputString[i] == ' ') {
      if (aDoubleIsBeingRead == true) {
        // We just finished reading the current string/double. Convert string to double now.
        char tmpVar[numberOfChars+1];
        for (std::string::size_type j = 0; j < numberOfChars; ++j) {
          tmpVar[j] = inputString[positionOfFirstChar+j];
        }
        tmpVar[numberOfChars] = '\0';
        outputDoubles.push_back(strtod(tmpVar,NULL));

        // Continue loop
        aDoubleIsBeingRead = false;
        positionOfFirstChar = 0;
        numberOfChars = 0;
      }
    }
    else {
      if (aDoubleIsBeingRead == false) {
        aDoubleIsBeingRead = true;
        positionOfFirstChar = i;
      }
      numberOfChars++;
    }
  } // for
  if (aDoubleIsBeingRead == true) {
    // We just finished reading the current string/double. Convert string to double now.
    char tmpVar[numberOfChars+1];
    for (std::string::size_type j = 0; j < numberOfChars; ++j) {
      tmpVar[j] = inputString[positionOfFirstChar+j];
    }
    tmpVar[numberOfChars] = '\0';
    outputDoubles.push_back(strtod(tmpVar,NULL));
  }
  std::vector<double>(outputDoubles).swap(outputDoubles);

  return;
}

void
MiscReadWordsFromString(
  const std::string&        inputString,
  std::vector<std::string>& outputWords)
{
  //std::cout << "In MiscReadWordsFromString()"
  //          << ": inputString = " << inputString
  //          << std::endl;
  outputWords.clear();
  bool aWordIsBeingRead = false;
  std::string::size_type positionOfFirstChar = 0;
  std::string::size_type numberOfChars = 0;
  for (std::string::size_type i = 0; i < inputString.size(); ++i) {
    queso_require_not_equal_to_msg(inputString[i], '\0', "character '\0' should not be found!");
    if (inputString[i] == ' ') {
      if (aWordIsBeingRead == true) {
        // We just finished reading the current string/word.
        char tmpVar[numberOfChars+1];
        for (std::string::size_type j = 0; j < numberOfChars; ++j) {
          tmpVar[j] = inputString[positionOfFirstChar+j];
        }
        tmpVar[numberOfChars] = '\0';
        outputWords.push_back(tmpVar);

        // Continue loop
        aWordIsBeingRead = false;
        positionOfFirstChar = 0;
        numberOfChars = 0;
      }
    }
    else {
      if (aWordIsBeingRead == false) {
        aWordIsBeingRead = true;
        positionOfFirstChar = i;
      }
      numberOfChars++;
    }
  } // for
  if (aWordIsBeingRead == true) {
    // We just finished reading the current string/word.
    char tmpVar[numberOfChars+1];
    for (std::string::size_type j = 0; j < numberOfChars; ++j) {
      tmpVar[j] = inputString[positionOfFirstChar+j];
    }
    tmpVar[numberOfChars] = '\0';
    outputWords.push_back(tmpVar);
  }
  std::vector<std::string>(outputWords).swap(outputWords);

  return;
}

//void
//MiscExtractDoubleFromString(
//  std::string& inputString,
//  double&      outputDouble)
//{
//  return;
//}

//void
//MiscExtractWordFromString(
//  std::string& inputString,
//  std::string& outputWord)
//{
//  return;
//}

int
MiscReadStringAndDoubleFromFile(
  std::ifstream& ifs,
  std::string&   termString,
  double*        termValue)
{
  int iRC = UQ_OK_RC;

  ifs >> termString;
  if ((ifs.rdstate() & std::ifstream::failbit)) {
    iRC = UQ_FAILED_READING_FILE_RC;
  }
  else if (termValue) {
    if (termString == std::string("inf")) {
      *termValue = INFINITY;
    }
    else if (termString == std::string("-inf")) {
      *termValue = -INFINITY;
    }
    else if (termString == std::string("nan")) {
      *termValue = nan("");
    }
    else {
      *termValue = strtod(termString.c_str(),NULL);
    }
  }
  //if (!iRC) std::cout << "Read termString = " << termString << std::endl;

  return iRC;
}

int
MiscReadCharsAndDoubleFromFile(
  std::ifstream& ifs,
  std::string&   termString,
  double*        termValue,
  bool&          endOfLineAchieved)
{
  int iRC = UQ_OK_RC;
  endOfLineAchieved = false;

  char c = ' ';
  while (c == ' ') {
    ifs.get(c);
    if ((ifs.rdstate() & std::ifstream::failbit)) {
      iRC = UQ_FAILED_READING_FILE_RC;
      break;
    }
  };

  char term[512];
  unsigned int pos = 0;

  if (!iRC) {
    while ((pos < 512) && (c != '\n') && (c != '\0') && (c != ' ')) {
      term[pos++] = c;
      if ((ifs.rdstate() & std::ifstream::failbit)) {
        iRC = UQ_FAILED_READING_FILE_RC;
        break;
      }
      ifs.get(c);
    };
  }

  if (!iRC) {
    if (c == '\n') endOfLineAchieved = true;
    term[pos] = '\0';
    termString = term;
    //std::cout << "Read chars = " << termString << std::endl;
    if (termValue) {
      if (termString == std::string("inf")) {
        *termValue = INFINITY;
      }
      else if (termString == std::string("-inf")) {
        *termValue = -INFINITY;
      }
      else if (termString == std::string("nan")) {
        *termValue = nan("");
      }
      else {
        *termValue = strtod(termString.c_str(),NULL);
      }
    }
  }

  return iRC;
}

double
MiscGammar(
  double                a,
  double                b,
  const RngBase* rngObject)
{
  double result = 0.;
  if (a < 1.) {
    result = MiscGammar(1.+a,b,rngObject)*std::pow( rngObject->uniformSample(),1./a );
  }
  else {
    double d = a-1./3.;
    double c = 1./std::sqrt(9.*d);
    double x = 0.;
    double w = 0.;
    while (1) {
      while (1) {
        x = rngObject->gaussianSample(1.);
        w = 1.+c*x;
        if (w > 0.) break;
      }
      w = std::pow(w,3.);
      double u = rngObject->uniformSample();
      double compValue = 1.-0.0331*std::pow(x,4.);
      if (u < compValue) break;
      compValue = 0.5*std::pow(x,2.)+d*(1.-w+log(w));
      if (log(u) < compValue) break;
    }
    result = b*d*w;
  }

  return result;
}

double
MiscGetEllapsedSeconds(struct timeval *timeval0)
{
  double result = 0.;

  struct timeval timevalNow;
  /*int iRC;*/
  /*iRC = */gettimeofday(&timevalNow, NULL);

  result  = (double) (timevalNow.tv_sec  - timeval0->tv_sec );
  result *= 1.e+6;
  result += (double) (timevalNow.tv_usec - timeval0->tv_usec);
  result *= 1.e-6;

  return result;
}

double MiscHammingWindow(unsigned int N, unsigned int j)
{
  double angle = 2.*M_PI*((double) j)/((double) N);
  double result = 0.53836 - 0.46164*cos(angle);

  return result;
}

double MiscGaussianDensity(double x, double mu, double sigma)
{
  double sigma2 = sigma*sigma;
  double diff   = x-mu;

  return (1./std::sqrt(2*M_PI*sigma2))*std::exp(-.5*diff*diff/sigma2);
}

unsigned int MiscUintDebugMessage(
  unsigned int value,
  const char*  message)
{
  if (message) {
    std::cout << "Passing in MiscUintDebugMessage(), value = " << value << ", message = " << message << std::endl;
  }
  return value;
}

int MiscIntDebugMessage(
  int         value,
  const char* message)
{
  if (message) {
    std::cout << "Passing in MiscIntDebugMessage(), value = " << value << ", message = " << message << std::endl;
  }
  return value;
}

double MiscDoubleDebugMessage(
  double     value,
  const char* message)
{
  if (message) {
    std::cout << "Passing in MiscDoubleDebugMessage(), value = " << value << ", message = " << message << std::endl;
  }
  return value;
}

///int CheckFilePath(const char *path)
///{
///
///  // verify parent directories in path exist (and create if not).
///
///#ifdef HAVE_GRVY
///  return(grvy_check_file_path(path));
///#else
///
///
///  return 0;
///#endif
///}

// ------------------------------------------
// Following routines borrowed from libGRVY
// ------------------------------------------

int CheckFilePath(const char *pathname)
  {

    // verify parent directories in path exist (and create if not).

#ifdef HAVE_GRVY
    return(grvy_check_file_path(pathname));
#else

    const int MAX_DEPTH = 50;

    char *pathlocal;
    char *parents;
    char *dirstring;
    char *token;
    int depth = 0;

    // Save a copy of pathname and look for the parent directories.

    pathlocal = strdup(pathname);
    dirstring = strdup(pathname);
    parents   = dirname(pathlocal);

    if(strcmp(parents,".") == 0)
      {
	free(pathlocal);
	free(dirstring);
	return 0;
      }

    // Deal with the possibility of an absolute path being provided

    bool abs_path = false;

    std::string leading_char("");
    std::string path_to_check;

    if(strncmp(parents,"/",1) == 0)
      {
	leading_char = "/";
	abs_path     = true;
      }

    // Verify existence of top-level directory

    if( (token = strtok(parents,"/")) != NULL )
      {
	path_to_check += leading_char + token;

	if ( GRVY_CheckDir(path_to_check.c_str()) )
	  {
	    free(pathlocal);
	    free(dirstring);
	    return -1;
	  }

	// Now, search for any remaining parent directories.

	if(abs_path)
	  sprintf(dirstring,"/%s",token);
	else
	  sprintf(dirstring,"%s",token);

	while ( (token = strtok(0,"/")) && (depth < MAX_DEPTH) )
	  {
	    dirstring = strcat(dirstring,"/");

	    if(GRVY_CheckDir(strcat(dirstring,token)))
	      {
		free(pathlocal);
		free(dirstring);
		return -1;
	      }
	    depth++;
	  };

	if(depth >= MAX_DEPTH )
	  {
	    std::cerr << __func__ << ": error - Max directory depth exceeded, limit =  " << MAX_DEPTH << std::endl;
	    free(pathlocal);
	    free(dirstring);
	    return -1;
	  }
      }

    // Clean Up
    free(pathlocal);
    free(dirstring);

    return 0;
#endif
  }


int GRVY_CheckDir(const char *dirname)
{
  struct stat st;

  if(stat(dirname,&st) != 0)
    {
      if( mkdir(dirname,0700) != 0 )
	{
	  std::cerr << __func__ << ": error - unable to create directory " << dirname << std::endl;
	  return -1;
	}
    }
  else if (!S_ISDIR(st.st_mode))
    {
      std::cerr << __func__ << ": error - entry exists, but is not a directory " << dirname << std::endl;
      return -1;
    }

  return 0;
}

template <class T>
bool
MiscCheckForSameValueInAllNodes(T&                    inputValue, // Yes, 'not' const
                                  double                acceptableTreshold,
                                  const MpiComm& comm,
                                  const char*           whereString)
{
  // Filter out those nodes that should not participate
  if (comm.MyPID() < 0) return true;

  double localValue = (double) inputValue;
  double sumValue = 0.;
  comm.Allreduce<double>(&localValue, &sumValue, (int) 1, RawValue_MPI_SUM,
                 whereString,
                 "failed MPI on 'sumValue' inside MiscCheckForSameValueInAllNodes()");

  double totalNumNodes = (double) comm.NumProc();
  double testValue = fabs(1. - localValue/(sumValue/totalNumNodes));
  unsigned int boolSum = 0;
#if 1
  unsigned int boolResult = 0;
  if (testValue > acceptableTreshold) boolResult = 1;
  comm.Allreduce<unsigned int>(&boolResult, &boolSum, (int) 1, RawValue_MPI_SUM,
                 whereString,
                 "failed MPI on 'boolSum' inside MiscCheckForSameValueInAllNodes()");

  if (boolSum > 0) {
    comm.Barrier();
    for (int i = 0; i < comm.NumProc(); ++i) {
      if (i == comm.MyPID()) {
        std::cerr << "WARNING, "
                  << whereString
                  << ", inside MiscCheckForSameValueInAllNodes()"
                  << ", rank (in this communicator) = " << i
                  << ": boolSum = "       << boolSum
                  << ", localValue = "    << localValue
                  << ", sumValue = "      << sumValue
                  << ", totalNumNodes = " << totalNumNodes
                  << ", avgValue = "      << (sumValue/totalNumNodes)
                  << ", relativeTest = "  << testValue
                  << std::endl;
      }
      comm.Barrier();
    }
    comm.Barrier();

    comm.Bcast((void *) &localValue, (int) 1, RawValue_MPI_DOUBLE, 0,
               whereString,
               "failed MPI on 'boolSum' inside MiscCheckForSameValueInAllNodes()");
    inputValue = localValue; // IMPORTANT
  }
#else
  queso_require_less_equal_msg(testValue, acceptableTreshold, "not all nodes have the same value inside MiscCheckForSameValueInAllNodes()");
#endif

  return (boolSum == 0);
}

template <class V>
void
MiscComputePositionsBetweenMinMax(V                minValues,
                                    V                maxValues,
                                    std::vector<V*>& positions)
{
  double factor = 0.5;
  switch (positions.size()) {
    case 0:
      // Do nothing
    break;

    case 1:
      positions[0] = new V((1. - factor) * minValues + factor * maxValues);
    break;

    default:
      for (unsigned int i = 0; i < positions.size(); ++i) {
        factor = ((double) i)/(((double) positions.size()) - 1.);
        positions[i] = new V((1. - factor) * minValues + factor * maxValues);
      }
    break;
  }

  return;
}

template <class V1,class V2>
void
MiscCheckTheParallelEnvironment(const V1& vec1, const V2& vec2)
{
  const BaseEnvironment& env = vec1.env();

  if (env.numSubEnvironments() == (unsigned int) env.fullComm().NumProc()) {
    queso_require_equal_to_msg(env.subRank(), 0, "there should exist only one processor per sub environment");
    queso_require_equal_to_msg(vec1.numOfProcsForStorage(), 1,
      "only 1 processor (per sub environment) should be necessary for the storage of a parameter vector");
    queso_require_equal_to_msg(vec2.numOfProcsForStorage(), 1,
      "only 1 processor (per sub environment) should be necessary for the storage of a parameter vector");
  }
  else if (env.numSubEnvironments() < (unsigned int) env.fullComm().NumProc()) {
    queso_require_equal_to_msg(env.fullComm().NumProc()%env.numSubEnvironments(), 0, "total number of processors should be a multiple of the number of sub environments");
    unsigned int numProcsPerSubEnvironment = env.fullComm().NumProc()/env.numSubEnvironments();
    queso_require_equal_to_msg(env.subComm().NumProc(), (int) numProcsPerSubEnvironment, "inconsistent number of processors per sub environment");
    if ((vec1.numOfProcsForStorage() == 1) &&
        (vec2.numOfProcsForStorage() == 1)) {
      // Ok
    }
    else if ((vec1.numOfProcsForStorage() == numProcsPerSubEnvironment) &&
             (vec2.numOfProcsForStorage() == numProcsPerSubEnvironment)) {
      queso_error_msg("parallel vectors are not supported yet");
    }
    else {
      queso_error_msg("number of processors required for a vector storage should be equal to either 1 or to the number of processors in the sub environment");
    }
  }
  else {
    queso_error_msg("number of processors per sub environment is less than 1!");
  }

  return;
}

}  // End namespace QUESO

template void QUESO::MiscCheckTheParallelEnvironment<QUESO::GslVector, QUESO::GslVector>(QUESO::GslVector const&, QUESO::GslVector const&);
template bool QUESO::MiscCheckForSameValueInAllNodes<bool>(bool&, double, QUESO::MpiComm const&, char const*);
template bool QUESO::MiscCheckForSameValueInAllNodes<double>(double&, double, QUESO::MpiComm const&, char const*);
