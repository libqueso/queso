#include <uqMiscellaneous.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <uqDefines.h>

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

void
uqMiscExtractDoubleFromString(
  std::string& inputString,
  double&      outputDouble)
{
  return;
}

void
uqMiscExtractWordFromString(
  std::string& inputString,
  std::string& outputWord)
{
  return;
}

double
uqMiscGammar(
  double   a,
  double   b,
  gsl_rng* rng)
{
  double result = 0.;
  if (a < 1.) {
    result = uqMiscGammar(1.+a,b,rng)*pow( gsl_rng_uniform(rng),1./a );
  }
  else {
    double d = a-1./3.;
    double c = 1./sqrt(9.*d);
    double x = 0.;
    double w = 0.;
    while (1) {
      while (1) {
        x = gsl_ran_gaussian(rng,1.);
        w = 1.+c*x;
        if (w > 0.) break;
      }
      w = pow(w,3.);
      double u = gsl_rng_uniform(rng);
      double compValue = 1.-0.0331*pow(x,4.);
      if (u < compValue) break;
      compValue = 0.5*pow(x,2.)+d*(1.-w+log(w));
      if (log(u) < compValue) break;
    }
    result = b*d*w;
  }

  return result;
}
