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

#include <uqEnvironment.h>
#include <EpetraExt_DistArray.h>

template <class V>
class uqAsciiTableClass
{
public:
  uqAsciiTableClass(const uqEnvironmentClass& env,
                          unsigned int        numRows,
                          unsigned int        numExtraCols,
                    const std::vector<bool>*  extraColIsString,
                    const std::string&        fileName);
 ~uqAsciiTableClass();

  unsigned int numRows  ()                                             const;
  unsigned int numCols  ()                                             const;
  void         getColumn(unsigned int                       columnId,
                         EpetraExt::DistArray<std::string>* stringCol,
                         V*                                 doubleCol) const;
  void         print    (std::ostream& os)                             const;

private:
  const uqEnvironmentClass&                       m_env;
  unsigned int                                    m_numRows;
  unsigned int                                    m_numCols;
  std::vector<bool>                               m_colIsString;
  std::vector<EpetraExt::DistArray<std::string>*> m_stringColumns;
  std::vector<V*>                                 m_doubleColumns;

  void readColumnsFromFile(std::string& fileName);
};

template <class V>
uqAsciiTableClass<V>::uqAsciiTableClass(
  const uqEnvironmentClass& env,
        unsigned int        numRows,
        unsigned int        numCols,
  const std::vector<bool>*  colIsString,
  const std::string&        fileName)
  :
  m_env          (env),
  m_numRows      (numRows),
  m_numCols      (numCols),
  m_colIsString  (1,true),
  m_stringColumns(1,NULL),
  m_doubleColumns(1,NULL)
{
}

template <class V>
uqAsciiTableClass<V>::~uqAsciiTableClass()
{
}

template <class V>
void
uqAsciiTableClass<V>::readColumnsFromFile(std::string& fileName)
{
  unsigned int maxCharsPerLine = 512;

  std::ifstream ifs(fileName.c_str());

  // Determine number of lines
  unsigned int numLines = std::count(std::istreambuf_iterator<char>(ifs),
                                     std::istreambuf_iterator<char>(),
                                     '\n');
#if 0
  // Determine number of components
  int iRC;
  ifs.seekg(0,std::ios_base::beg);
  unsigned int lineId = 0;
  unsigned int numComponents = 0;
  std::string tempString;
  while ((lineId < numLines) && (ifs.eof() == false)) {
    iRC = uqMiscReadStringAndDoubleFromFile(ifs,tempString,NULL);
    UQ_FATAL_TEST_MACRO(iRC,
                        m_env.rank(),
                        "uqAsciiTableClass<V>::constructor()",
                        "failed reading during the determination of the number of components");
    //std::cout << "lineId = "           << lineId
    //          << ", numComponents = " << numComponents
    //          << ", tempString = "     << tempString
    //          << std::endl;
    if (tempString[0] != '#') numComponents++;
    lineId++;
    ifs.ignore(maxCharsPerLine,'\n');
  }
  UQ_FATAL_TEST_MACRO(lineId != numLines,
                      m_env.rank(),
                      "uqAsciiTableClass<V>::constructor()",
                      "the first number of lines read is nonconsistent");
  if (m_dim != numComponents) {
    char errorExplanation[512];
    sprintf(errorExplanation,"number of components (%d) in space spec file does not match dimension (%d) passed in the main input file",numComponents,m_dim);
    UQ_FATAL_TEST_MACRO(true,
                        m_env.rank(),
                        "uqAsciiTableClass<V>::constructor()",
                        errorExplanation);
  }

  std::cout << "Space spec file '"     << fileName
            << "' has "                << numLines
            << " lines and specifies " << numComponents
            << " components."
            << std::endl;
  m_componentsNames->resize(numComponents,"");

  // Read file until End Of File character is reached
  ifs.seekg(0,std::ios_base::beg);
  lineId = 0;
  unsigned int componentId = 0;
  std::string  componentName("");
  while ((lineId < numLines) && (ifs.eof() == false)) {
    //std::cout << "Beginning read of line (in space spec file) of id = " << lineId << std::endl;
    bool endOfLineAchieved = false;

    iRC = uqMiscReadCharsAndDoubleFromFile(ifs, componentName, NULL, endOfLineAchieved);
    UQ_FATAL_TEST_MACRO(iRC,
                        m_env.rank(),
                        "uqAsciiTableClass<V>::constructor()",
                        "failed reading a component name during the components names reading loop");

    lineId++;
    if (componentName[0] == '#') {
      if (!endOfLineAchieved) ifs.ignore(maxCharsPerLine,'\n');
      continue;
    }

    // Check 'componentId' before setting one more component
    if (componentId >= m_componentsNames->size()) {
      char errorExplanation[512];
      sprintf(errorExplanation,"componentId (%d) got too large during reading of space spec file",componentId);
      UQ_FATAL_TEST_MACRO(true,
                          m_env.rank(),
                          "uqAsciiTableClass<V>::constructor()",
                          errorExplanation);
    }

    std::cout << "Just read, for componentId = " << componentId
              << ": componentName = "            << componentName
              << std::endl;
    setComponentName(componentId,
                     componentName);
    componentId++;
  }

  UQ_FATAL_TEST_MACRO(lineId != numLines,
                      m_env.rank(),
                      "uqAsciiTableClass<V>::constructor()",
                      "the second number of lines read is nonconsistent");
  UQ_FATAL_TEST_MACRO(componentId != m_componentsNames->size(),
                      m_env.rank(),
                      "uqAsciiTableClass<V>::constructor()",
                      "the number of components just read is nonconsistent");
#endif
  return;
}

template <class V>
unsigned int
uqAsciiTableClass<V>::numRows() const
{
  return m_numRows;
}

template <class V>
unsigned int
uqAsciiTableClass<V>::numCols() const
{
  return m_numCols;
}

template <class V>
void
uqAsciiTableClass<V>::getColumn(
  unsigned int                       columnId,
  EpetraExt::DistArray<std::string>* stringCol,
  V*                                 doubleCol) const
{
  return;
}

template <class V>
void
uqAsciiTableClass<V>::print(std::ostream& os) const
{
  os << "\nTables contents are:"
     << std::endl;
  for (unsigned int i = 0; i < m_numRows; ++i) {
    UQ_FATAL_TEST_MACRO((m_stringColumns[i] != NULL) && (m_doubleColumns[i] != NULL),
                        m_env.rank(),
                        "uqAsciiTableClass<V>::print()",
                        "column is not null on both possible ways");
    UQ_FATAL_TEST_MACRO((m_stringColumns[i] == NULL) && (m_doubleColumns[i] == NULL),
                        m_env.rank(),
                        "uqAsciiTableClass<V>::print()",
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

template<class V>
std::ostream& operator<<(std::ostream& os, const uqAsciiTableClass<V>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_ASCII_TABLE_H__

