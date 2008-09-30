/* uq/libs/basic/inc/uqAsciiTable.h
 *
 * Copyright (C) 2008 The PECOS Team, http://queso.ices.utexas.edu
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

#ifndef __UQ_ASCII_TABLE_H__
#define __UQ_ASCII_TABLE_H__

#include <uqVectorSpace.h>
#include <uqMiscellaneous.h>
#include <uqEnvironment.h>
#include <EpetraExt_DistArray.h>

template <class V, class M>
class uqAsciiTableClass
{
public:
  uqAsciiTableClass(const uqEnvironmentClass&      env,
                    const uqVectorSpaceClass<V,M>& vectorSpace,
                          unsigned int             numExtraCols,
                    const std::vector<bool>*       extraColIsString,
                    const std::string&             fileName);
 ~uqAsciiTableClass();

  unsigned int numRows  ()                                             const;
  unsigned int numCols  ()                                             const;
  void         getColumn(unsigned int                       columnId,
                         EpetraExt::DistArray<std::string>* stringCol,
                         V*                                 doubleCol) const;
  void         print    (std::ostream& os)                             const;

private:
  const uqEnvironmentClass&                       m_env;
  const uqVectorSpaceClass<V,M>&                  m_vectorSpace;
  unsigned int                                    m_numCols;
  std::vector<bool>                               m_colIsString;
  std::vector<EpetraExt::DistArray<std::string>*> m_stringColumns;
  std::vector<V*>                                 m_doubleColumns;

  void readColumnsFromFile(std::string& fileName);
};

template <class V, class M>
uqAsciiTableClass<V,M>::uqAsciiTableClass(
  const uqEnvironmentClass&      env,
  const uqVectorSpaceClass<V,M>& vectorSpace,
        unsigned int             numExtraCols,
  const std::vector<bool>*       extraColIsString,
  const std::string&             fileName)
  :
  m_env          (env),
  m_vectorSpace  (vectorSpace),
  m_numCols      (1+numExtraCols),
  m_colIsString  (1,true),
  m_stringColumns(0),
  m_doubleColumns(0)
{
  if (m_numCols > 1) {
    m_colIsString.resize(m_numCols,false);
    if (extraColIsString == NULL) {
      // Nothing extra needs to be done
    }
    else {
      unsigned int maxJ = std::min(numExtraCols,extraColIsString->size());
      for (unsigned int j = 0; j < maxJ; ++j) {
        m_colIsString[1+j] = (*extraColIsString)[j];
      }
    }
  }
  m_stringColumns.resize(m_numCols,NULL);
  m_doubleColumns.resize(m_numCols,NULL);
}

template <class V, class M>
uqAsciiTableClass<V,M>::~uqAsciiTableClass()
{
  for (unsigned int j = 0; j < m_doubleColumns.size(); ++j) {
    if (m_doubleColumns[j]) delete m_doubleColumns[j];
  }
  for (unsigned int j = 0; j < m_stringColumns.size(); ++j) {
    if (m_stringColumns[j]) delete m_stringColumns[j];
  }
}

template <class V, class M>
void
uqAsciiTableClass<V,M>::readColumnsFromFile(std::string& fileName)
{
  unsigned int maxCharsPerLine = 512;

  std::ifstream ifs(fileName.c_str());

  // Determine number of lines
  unsigned int numLines = std::count(std::istreambuf_iterator<char>(ifs),
                                     std::istreambuf_iterator<char>(),
                                     '\n');

  // Determine number of valid lines
  int iRC;
  ifs.seekg(0,std::ios_base::beg);
  unsigned int lineId = 0;
  unsigned int numValid = 0;
  std::string tempString;
  while ((lineId < numLines) && (ifs.eof() == false)) {
    iRC = uqMiscReadStringAndDoubleFromFile(ifs,tempString,NULL);
    UQ_FATAL_TEST_MACRO(iRC,
                        m_env.rank(),
                        "uqAsciiTableClass<V,M>::constructor()",
                        "failed reading during the determination of the number of valid lines");
    //std::cout << "lineId = "           << lineId
    //          << ", numValid = " << numValid
    //          << ", tempString = "     << tempString
    //          << std::endl;
    if (tempString[0] != '#') numValid++;
    lineId++;
    ifs.ignore(maxCharsPerLine,'\n');
  }
  UQ_FATAL_TEST_MACRO(lineId != numLines,
                      m_env.rank(),
                      "uqAsciiTableClass<V,M>::constructor()",
                      "the first number of lines read is nonconsistent");
  if (m_vectorSpace.dim() != numValid) {
    char errorExplanation[512];
    sprintf(errorExplanation,"number of valid lines (%d) in ascii table file does not match dimension (%d) of vector space",numValid,m_vectorSpace.dim());
    UQ_FATAL_TEST_MACRO(true,
                        m_env.rank(),
                        "uqAsciiTableClass<V,M>::constructor()",
                        errorExplanation);
  }

  std::cout << "Ascii table file '"    << fileName
            << "' has "                << numLines
            << " lines and specifies " << numValid
            << " valid lines."
            << std::endl;

  for (unsigned int j=0; j < m_numCols; +j) {
    if (m_colIsString[j]) {
      m_stringColumns[j] = new EpetraExt::DistArray<std::string>(m_vectorSpace.map(),1);
    }
    else {
      m_doubleColumns[j] = new V(m_vectorSpace.zeroVector());
    }
  }

  // Read file until End Of File character is reached
  ifs.seekg(0,std::ios_base::beg);
  lineId = 0;
  unsigned int validId = 0;
  std::string  firstColumn("");
  while ((lineId < numLines) && (ifs.eof() == false)) {
    //std::cout << "Beginning read of line (in ascii table file) of id = " << lineId << std::endl;
    bool endOfLineAchieved = false;

    iRC = uqMiscReadCharsAndDoubleFromFile(ifs, firstColumn, NULL, endOfLineAchieved);
    UQ_FATAL_TEST_MACRO(iRC,
                        m_env.rank(),
                        "uqAsciiTableClass<V,M>::constructor()",
                        "failed reading a first column during the valid lines reading loop");

    lineId++;
    if (firstColumn[0] == '#') {
      if (!endOfLineAchieved) ifs.ignore(maxCharsPerLine,'\n');
      continue;
    }

    // Check 'validId' before setting one more valid line
    if (validId >= numValid) {
      char errorExplanation[512];
      sprintf(errorExplanation,"validId (%d) got too large during reading of ascii table file",validId);
      UQ_FATAL_TEST_MACRO(true,
                          m_env.rank(),
                          "uqAsciiTableClass<V,M>::constructor()",
                          errorExplanation);
    }

    std::cout << "Just read, for validId = " << validId
              << ": firstColumn = "          << firstColumn
              << std::endl;
    //setComponentName(validId,
    //                 firstColumn);
    validId++;
  }

  UQ_FATAL_TEST_MACRO(lineId != numLines,
                      m_env.rank(),
                      "uqAsciiTableClass<V,M>::constructor()",
                      "the second number of lines read is not consistent");
  UQ_FATAL_TEST_MACRO(validId != numValid,
                      m_env.rank(),
                      "uqAsciiTableClass<V,M>::constructor()",
                      "the number of valid lines just read is not consistent");

  return;
}

template <class V, class M>
unsigned int
uqAsciiTableClass<V,M>::numRows() const
{
  return m_vectorSpace.dim();
}

template <class V, class M>
unsigned int
uqAsciiTableClass<V,M>::numCols() const
{
  return m_numCols;
}

template <class V, class M>
void
uqAsciiTableClass<V,M>::getColumn(
  unsigned int                       columnId,
  EpetraExt::DistArray<std::string>* stringCol,
  V*                                 doubleCol) const
{
  return;
}

template <class V, class M>
void
uqAsciiTableClass<V,M>::print(std::ostream& os) const
{
  os << "\nTables contents are:"
     << std::endl;
  for (unsigned int i = 0; i < m_vectorSpace.dim(); ++i) {
    UQ_FATAL_TEST_MACRO((m_stringColumns[i] != NULL) && (m_doubleColumns[i] != NULL),
                        m_env.rank(),
                        "uqAsciiTableClass<V,M>::print()",
                        "column is not null on both possible ways");
    UQ_FATAL_TEST_MACRO((m_stringColumns[i] == NULL) && (m_doubleColumns[i] == NULL),
                        m_env.rank(),
                        "uqAsciiTableClass<V,M>::print()",
                        "column is null on both possible ways");

    os << i << " ";
    if (m_stringColumns[i] != NULL) {
      os << *(m_stringColumns[i]);
    }
    else {
      os << *(m_doubleColumns[i]);
    }
    os << std::endl;
  }

  return;
}

template<class V, class M>
std::ostream& operator<<(std::ostream& os, const uqAsciiTableClass<V,M>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_ASCII_TABLE_H__

