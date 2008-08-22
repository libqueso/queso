/* uq/libs/mcmc/inc/uqObservableSpace.h
 *
 * Copyright (C) 2008 The PECOS Team, http://www.ices.utexas.edu/centers/pecos
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_OBSERVABLE_SPACE_H__
#define __UQ_OBSERVABLE_SPACE_H__

#include <uqFinDimLinearSpace.h>
#include <uqMiscellaneous.h>
#include <uqObservable.h>
#include <vector>

template <class V, class M>
class uqObservableSpaceClass : public uqFinDimLinearSpaceClass<V,M>
{
public:
           uqObservableSpaceClass();
           uqObservableSpaceClass(const uqEnvironmentClass& env,
                              const char*               prefix);
  virtual ~uqObservableSpaceClass();

          bool                     shouldVariancesBeUpdated   () const;
          unsigned int             dim                        () const;
          int                      setObservable              (unsigned int       observableId,
                                                               const std::string& name,
                                                               unsigned int       numberOfObservations,
                                                               double             priorVariance = 1.,
                                                               double             varianceAccuracy = 1.);
          const uqObservableClass& observable                 (unsigned int observableId) const;
          const V&                 numbersOfObservations      () const;
          const V&                 priorVariances             () const;
          const V&                 varianceAccuracies         () const;
  virtual void                     print                      (std::ostream& os) const;
          void                     printObservableNames       (std::ostream& os, bool printHorizontally) const; // See template specialization

protected:
          void                     defineMyOptions            (po::options_description& optionsDesc) const;
          void                     getMyOptionValues          (po::options_description& optionsDesc);
          void                     readObservablesFromSpecFile(std::string& specFileName);
          void                     resetValues                ();
          void                     createNumbersOfObservations() const; // See template specialization
          void                     createPriorVariances       () const; // See template specialization
          void                     createVarianceAccuracies   () const; // See template specialization

  po::options_description* m_optionsDesc;
  std::string              m_option_help;
  std::string              m_option_updateVariances;
  std::string              m_option_numSubSpaces;
  std::string              m_option_specificationFile;

  std::vector<uqObservableClass*> m_observables; // FIXME: will need to be a parallel vector in case of a very large number of observables
  uqObservableClass               m_dummyObservable;
  bool                            m_updateVariances;
  mutable V*                      m_numbersOfObservations;
  mutable V*                      m_priorVariances;
  mutable V*                      m_varianceAccuracies;

  using uqFinDimLinearSpaceClass<V,M>::m_env;
  using uqFinDimLinearSpaceClass<V,M>::m_dim;
};

template <class V, class M>
uqObservableSpaceClass<V,M>::uqObservableSpaceClass()
  :
  uqFinDimLinearSpaceClass<V,M>()
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqObservableSpaceClass<V,M>::constructor(), default",
                      "should not be used by user");
}

template <class V, class M>
uqObservableSpaceClass<V,M>::uqObservableSpaceClass(
  const uqEnvironmentClass& env,
  const char*               prefix)
  :
  uqFinDimLinearSpaceClass<V,M>(env,prefix),
  m_optionsDesc                (NULL),
  m_observables                (0),//,NULL),
  m_dummyObservable            ("NonExistentYet",0),
  m_updateVariances            (false),
  m_numbersOfObservations      (0),
  m_priorVariances             (NULL),
  m_varianceAccuracies         (NULL)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqObservableSpaceClass<V,M>::constructor()"
              << std::endl;
  }

  m_option_help              = uqFinDimLinearSpaceClass<V,M>::m_prefix + "observableSpace_help";
  m_option_updateVariances   = uqFinDimLinearSpaceClass<V,M>::m_prefix + "observableSpace_updateVariances";
  m_option_numSubSpaces      = uqFinDimLinearSpaceClass<V,M>::m_prefix + "observableSpace_numSubSpaces";
  m_option_specificationFile = uqFinDimLinearSpaceClass<V,M>::m_prefix + "observableSpace_specificationFile";

  m_optionsDesc = new po::options_description("Observable space options");
  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  // Now that 'm_dim' has been set, construct Trilinos map
  this->constructMap();

  if (m_env.rank() == 0) std::cout << "After getting values of options with prefix '" << uqFinDimLinearSpaceClass<V,M>::m_prefix
                                   << "', state of uqObservableSpaceClass object is:"
                                   << "\n" << *this
                                   << std::endl;

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqObservableSpaceClass<V,M>::constructor()"
              << std::endl;
  }
}

template <class V, class M>
uqObservableSpaceClass<V,M>::~uqObservableSpaceClass()
{
  //std::cout << "Entering uqObservableSpaceClass<V,M>::destructor()"
  //          << std::endl;

  if (m_varianceAccuracies   ) delete m_varianceAccuracies;
  if (m_priorVariances       ) delete m_priorVariances;
  if (m_numbersOfObservations) delete m_numbersOfObservations;

  for (unsigned int i = 0; i < m_observables.size(); ++i) {
    if (m_observables[i]) delete m_observables[i];
  }

  if (m_optionsDesc) delete m_optionsDesc;

  //std::cout << "Leaving uqObservableSpaceClass<V,M>::destructor()"
  //          << std::endl;
}

template <class V, class M>
bool
uqObservableSpaceClass<V,M>::shouldVariancesBeUpdated() const
{
  return m_updateVariances;
}

template <class V, class M>
unsigned int
uqObservableSpaceClass<V,M>::dim() const
{
  return m_dim;
}

template <class V, class M>
void
uqObservableSpaceClass<V,M>::defineMyOptions(
  po::options_description& optionsDesc) const
{
  m_optionsDesc->add_options()
    (m_option_help.c_str(),                                                            "produce help message for UQ observable space"                 )
    (m_option_updateVariances.c_str(),   po::value<bool        >()->default_value(0),  "update variance of likelihood results"                        )
    (m_option_numSubSpaces.c_str(),      po::value<unsigned int>()->default_value(0),  "number of vectors of observable quantities"                   )
    (m_option_specificationFile.c_str(), po::value<std::string >()->default_value(""), "File with the specification of all observables to be computed")
  ;

  return;
}

template <class V, class M>
void
uqObservableSpaceClass<V,M>::getMyOptionValues(po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    std::cout << optionsDesc
              << std::endl;
  }

  if (m_env.allOptionsMap().count(m_option_updateVariances.c_str())) {
    const po::variables_map& tmpMap = m_env.allOptionsMap();
    m_updateVariances = tmpMap[m_option_updateVariances.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_numSubSpaces.c_str())) {
    const po::variables_map& tmpMap = m_env.allOptionsMap();
    m_dim = tmpMap[m_option_numSubSpaces.c_str()].as<unsigned int>();
  }

  // Read observable specification file only if 0 dimension was passed to constructor
  if (m_observables.size() == 0) {
    std::string specFileName("");
    if (m_env.allOptionsMap().count(m_option_specificationFile.c_str())) {
      const po::variables_map& tmpMap = m_env.allOptionsMap();
      specFileName = tmpMap[m_option_specificationFile.c_str()].as<std::string>();
      readObservablesFromSpecFile(specFileName);
    }
  }

  return;
}

template <class V, class M>
void
uqObservableSpaceClass<V,M>::readObservablesFromSpecFile(std::string& specFileName)
{
  unsigned int maxCharsPerLine = 512;

  std::ifstream ifs(specFileName.c_str());

  // Determine number of lines
  unsigned int numLines = std::count(std::istreambuf_iterator<char>(ifs),
                                     std::istreambuf_iterator<char>(),
                                     '\n');

  // Determine number of observables
  int iRC;
  ifs.seekg(0,std::ios_base::beg);
  unsigned int lineId = 0;
  unsigned int numObservables = 0;
  std::string tempString;
  while ((lineId < numLines) && (ifs.eof() == false)) {
    iRC = uqMiscReadStringAndDoubleFromFile(ifs,tempString,NULL);
    UQ_FATAL_TEST_MACRO(iRC,
                        m_env.rank(),
                        "uqObservableSpaceClass<V,M>::constructor()",
                        "failed reading during the determination of the number of observables");
    //std::cout << "lineId = "           << lineId
    //          << ", numObservables = " << numObservables
    //          << ", tempString = "     << tempString
    //          << std::endl;
    if (tempString[0] != '#') numObservables++;
    lineId++;
    ifs.ignore(maxCharsPerLine,'\n');
  }
  UQ_FATAL_TEST_MACRO(lineId != numLines,
                      m_env.rank(),
                      "uqObservableSpaceClass<V,M>::constructor()",
                      "the first number of lines read is nonconsistent");
  if (m_dim != numObservables) {
    char errorExplanation[512];
    sprintf(errorExplanation,"number of observables (%d) in observable specification file does not match dimension (%d) passed in the main input file",numObservables,m_dim);
    UQ_FATAL_TEST_MACRO(true,
                        m_env.rank(),
                        "uqObservableSpaceClass<V,M>::constructor()",
                        errorExplanation);
  }

  std::cout << "Observable specification file '" << specFileName
            << "' has "                          << numLines
            << " lines and specifies "           << numObservables
            << " observables."
            << std::endl;
  m_observables.resize(numObservables,NULL);

  // Read file until End Of File character is reached
  ifs.seekg(0,std::ios_base::beg);
  lineId = 0;
  unsigned int observableId = 0;
  std::string  observableName            ("");
  std::string  numberOfObservationsString("");
  std::string  priorVarianceString       ("");
  std::string  varianceAccuracyString    ("");
  unsigned int numberOfObservations;
  double       priorVariance;
  double       varianceAccuracy;
  while ((lineId < numLines) && (ifs.eof() == false)) {
    //std::cout << "Beginning read of line (in observable specification file) of id = " << lineId << std::endl;
    bool endOfLineAchieved = false;

    iRC = uqMiscReadCharsAndDoubleFromFile(ifs, observableName, NULL, endOfLineAchieved);
    UQ_FATAL_TEST_MACRO(iRC,
                        m_env.rank(),
                        "uqObservableSpaceClass<V,M>::constructor()",
                        "failed reading a observable name during the observables reading loop");

    lineId++;
    if (observableName[0] == '#') {
      if (!endOfLineAchieved) ifs.ignore(maxCharsPerLine,'\n');
      continue;
    }
    UQ_FATAL_TEST_MACRO(endOfLineAchieved,
                        m_env.rank(),
                        "uqObservableSpaceClass<V,M>::constructor()",
                        "failed to provide information beyond observable name during the observables reading loop");

    // Check 'observableId' before setting one more observable
    if (observableId >= m_observables.size()) {
      char errorExplanation[512];
      sprintf(errorExplanation,"observableId (%d) got too large during reading of observable specification file",observableId);
      UQ_FATAL_TEST_MACRO(true,
                          m_env.rank(),
                          "uqObservableSpaceClass<V,M>::constructor()",
                          errorExplanation);
    }

    double tmpDouble;
    iRC = uqMiscReadCharsAndDoubleFromFile(ifs, numberOfObservationsString, &tmpDouble, endOfLineAchieved);
    UQ_FATAL_TEST_MACRO(iRC,
                        m_env.rank(),
                        "uqObservableSpaceClass<V,M>::constructor()",
                        "failed reading an initial value during the observables reading loop");
    numberOfObservations = (unsigned int) tmpDouble;

    priorVariance    = 1.;
    varianceAccuracy = 1.;

    if (!endOfLineAchieved) {
      iRC = uqMiscReadCharsAndDoubleFromFile(ifs, priorVarianceString, &priorVariance, endOfLineAchieved);
      UQ_FATAL_TEST_MACRO(iRC,
                          m_env.rank(),
                          "uqObservableSpaceClass<V,M>::constructor()",
                          "failed reading a minimal value during the observables reading loop");
    }

    if (!endOfLineAchieved) {
      iRC = uqMiscReadCharsAndDoubleFromFile(ifs, varianceAccuracyString, &varianceAccuracy, endOfLineAchieved);
      UQ_FATAL_TEST_MACRO(iRC,
                          m_env.rank(),
                          "uqObservableSpaceClass<V,M>::constructor()",
                          "failed reading a maximum value during the observables reading loop");
    }

    if (!endOfLineAchieved) ifs.ignore(maxCharsPerLine,'\n');
    //std::cout << "Just read, for observableId = " << observableId
    //          << ": observableName = "            << observableName
    //          << ", numberOfObservations = "      << numberOfObservations
    //          << ", priorVariance = "             << priorVariance
    //          << ", varianceAccuracy = "          << varianceAccuracy
    //          << std::endl;
    setObservable(observableId,
                  observableName,
                  numberOfObservations,
                  priorVariance,
                  varianceAccuracy);
    observableId++;
  }

  UQ_FATAL_TEST_MACRO(lineId != numLines,
                      m_env.rank(),
                      "uqObservableSpaceClass<V,M>::constructor()",
                      "the second number of lines read is nonconsistent");
  UQ_FATAL_TEST_MACRO(observableId != m_observables.size(),
                      m_env.rank(),
                      "uqObservableSpaceClass<V,M>::constructor()",
                      "the number of observables just read is nonconsistent");
  return;
}

template <class V, class M>
int
uqObservableSpaceClass<V,M>::setObservable(
  unsigned int       observableId,
  const std::string& name,
  unsigned int       numberOfObservations,
  double             priorVariance,
  double             varianceAccuracy)
{
  UQ_TEST_MACRO((observableId > m_observables.size()),
                m_env.rank(),
                "uqObservableSpaceClass<V,M>::setObservable()",
                "observableId is too big",
                UQ_INVALID_OBSERVABLE_SPEC_RC);

  if (m_observables[observableId] == NULL) {
    m_observables[observableId] = new uqObservableClass(name,
                                                        numberOfObservations,
                                                        priorVariance,
                                                        varianceAccuracy);
  }
  else {
    m_observables[observableId]->setName                (name);
    m_observables[observableId]->setNumberOfObservations(numberOfObservations);
    m_observables[observableId]->setPriorVariance       (priorVariance);
    m_observables[observableId]->setVarianceAccuracy    (varianceAccuracy);
  }

  // These values cannot be trusted anymore
  // They need to be updated
  // They will be updated the next time they are requested
  resetValues();

  return 0;
}

template <class V, class M>
const uqObservableClass&
uqObservableSpaceClass<V,M>::observable(unsigned int observableId) const
{
  if (observableId > m_observables.size()) return m_dummyObservable;
  if (m_observables[observableId] == NULL) return m_dummyObservable;
  return *(m_observables[observableId]);
}

template <class V, class M>
void
uqObservableSpaceClass<V,M>::resetValues()
{
  if (m_varianceAccuracies   ) delete m_varianceAccuracies;
  if (m_priorVariances       ) delete m_priorVariances;
  if (m_numbersOfObservations) delete m_numbersOfObservations;
  m_varianceAccuracies    = NULL;
  m_priorVariances        = NULL;
  m_numbersOfObservations = NULL;
}

template <class V, class M>
const V&
uqObservableSpaceClass<V,M>::numbersOfObservations() const
{
  if (m_numbersOfObservations == NULL) this->createNumbersOfObservations();
  return *m_numbersOfObservations;
}

template <class V, class M>
const V&
uqObservableSpaceClass<V,M>::priorVariances() const
{
  if (m_priorVariances == NULL) this->createPriorVariances();
  return *m_priorVariances;
}

template <class V, class M>
const V&
uqObservableSpaceClass<V,M>::varianceAccuracies() const
{
  if (m_varianceAccuracies == NULL) this->createVarianceAccuracies();
  return *m_varianceAccuracies;
}

template<class V, class M>
void
  uqObservableSpaceClass<V,M>::printObservableNames(std::ostream& os, bool printHorizontally) const
{
  if (printHorizontally) { 
    for (unsigned int i = 0; i < this->dim(); ++i) {
      os << "'" << this->observable(i).name() << "'"
         << " ";
    }
  }
  else {
    for (unsigned int i = 0; i < this->dim(); ++i) {
      os << "'" << this->observable(i).name() << "'"
         << std::endl;
    }
  }

  return;
}

template <class V, class M>
void
uqObservableSpaceClass<V,M>::print(std::ostream& os) const
{
  os <<         uqFinDimLinearSpaceClass<V,M>::m_prefix << "updateVariances = " << m_updateVariances
     << "\n" << uqFinDimLinearSpaceClass<V,M>::m_prefix << "numSubspaces = "    << m_dim
     << "\nObservables are:"
     << std::endl;
  for (unsigned int i = 0; i < this->dim(); ++i) {
    os << i << " ";
    if (m_observables[i]) {
      os << *(m_observables[i]);
    }
    else {
      os << "NULL";
    }
    os << std::endl;
  }

  return;
}

template<class V, class M>
std::ostream&
operator<<(std::ostream& os, const uqObservableSpaceClass<V,M>& space)
{
  space.print(os);

  return os;
}
#endif // __UQ_OBSERVABLE_SPACE_H__

