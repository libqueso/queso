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

#include <uqMiscellaneous.h>
#include <gsl/gsl_randist.h>
#include <sys/time.h>
#include <math.h>
#include <uqDefines.h>
#include <iostream>
#include <fstream>

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
  double         a,
  double         b,
  const gsl_rng* rng)
{
  double result = 0.;
  if (a < 1.) {
    result = uqMiscGammar(1.+a,b,rng)*std::pow( gsl_rng_uniform(rng),1./a );
  }
  else {
    double d = a-1./3.;
    double c = 1./std::sqrt(9.*d);
    double x = 0.;
    double w = 0.;
    while (1) {
      while (1) {
        x = gsl_ran_gaussian(rng,1.);
        w = 1.+c*x;
        if (w > 0.) break;
      }
      w = std::pow(w,3.);
      double u = gsl_rng_uniform(rng);
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
