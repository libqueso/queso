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

#include <uqDefines.h>
#include <uqMiscellaneous.h>
#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <libgen.h>
#include <sys/stat.h>
#include <cmath>

void
uqMiscReadDoublesFromString(
  const std::string&         inputString,
        std::vector<double>& outputDoubles)
{
  //std::cout << "In uqMiscReadDoublesFromString()"
  //          << ": inputString = " << inputString
  //          << std::endl;
  outputDoubles.clear();
  bool aDoubleIsBeingRead = false;
  std::string::size_type positionOfFirstChar = 0;
  std::string::size_type numberOfChars = 0;
  for (std::string::size_type i = 0; i < inputString.size(); ++i) {
    UQ_FATAL_TEST_MACRO((inputString[i] == '\0'),
                        UQ_UNAVAILABLE_RANK,
                        "uqMiscReadDoublesFromString()",
                        "character '\0' should not be found!");
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
uqMiscReadWordsFromString(
  const std::string&        inputString,
  std::vector<std::string>& outputWords)
{
  //std::cout << "In uqMiscReadWordsFromString()"
  //          << ": inputString = " << inputString
  //          << std::endl;
  outputWords.clear();
  bool aWordIsBeingRead = false;
  std::string::size_type positionOfFirstChar = 0;
  std::string::size_type numberOfChars = 0;
  for (std::string::size_type i = 0; i < inputString.size(); ++i) {
    UQ_FATAL_TEST_MACRO((inputString[i] == '\0'),
                        UQ_UNAVAILABLE_RANK,
                        "uqMiscReadWordsFromString()",
                        "character '\0' should not be found!");
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
//uqMiscExtractDoubleFromString(
//  std::string& inputString,
//  double&      outputDouble)
//{
//  return;
//}

//void
//uqMiscExtractWordFromString(
//  std::string& inputString,
//  std::string& outputWord)
//{
//  return;
//}

int
uqMiscReadStringAndDoubleFromFile(
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
uqMiscReadCharsAndDoubleFromFile(
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
uqMiscGammar(
  double                a,
  double                b,
  const uqRngBaseClass* rngObject)
{
  double result = 0.;
  if (a < 1.) {
    result = uqMiscGammar(1.+a,b,rngObject)*std::pow( rngObject->uniformSample(),1./a );
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
uqMiscGetEllapsedSeconds(struct timeval *timeval0)
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

double uqMiscHammingWindow(unsigned int N, unsigned int j)
{
  double angle = 2.*M_PI*((double) j)/((double) N);
  double result = 0.53836 - 0.46164*cos(angle);

  return result;
}

double uqMiscGaussianDensity(double x, double mu, double sigma)
{
  double sigma2 = sigma*sigma;
  double diff   = x-mu;

  return (1./std::sqrt(2*M_PI*sigma2))*std::exp(-.5*diff*diff/sigma2);
}

unsigned int uqMiscUintDebugMessage(
  unsigned int value,
  const char*  message)
{
  if (message) {
    std::cout << "Passing in uqMiscUintDebugMessage(), value = " << value << ", message = " << message << std::endl;
  }
  return value;
}

int uqMiscIntDebugMessage(
  int         value,
  const char* message)
{
  if (message) {
    std::cout << "Passing in uqMiscIntDebugMessage(), value = " << value << ", message = " << message << std::endl;
  }
  return value;
}

double uqMiscDoubleDebugMessage(
  double     value,
  const char* message)
{
  if (message) {
    std::cout << "Passing in uqMiscDoubleDebugMessage(), value = " << value << ", message = " << message << std::endl;
  }
  return value;
}

///int uqCheckFilePath(const char *path) 
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

int uqCheckFilePath(const char *pathname)
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

	if ( uqGRVY_CheckDir(path_to_check.c_str()) )
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

	    if(uqGRVY_CheckDir(strcat(dirstring,token)))
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


int uqGRVY_CheckDir(const char *dirname)
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


