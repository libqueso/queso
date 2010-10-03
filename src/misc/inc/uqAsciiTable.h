//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010 The PECOS Development Team
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

#ifndef __UQ_ASCII_TABLE_H__
#define __UQ_ASCII_TABLE_H__

#include <uqEnvironment.h>
#include <uqMiscellaneous.h>
#include <EpetraExt_DistArray.h>

template <class V, class M>
class uqAsciiTableClass
{
public:
  uqAsciiTableClass(const uqBaseEnvironmentClass& env,
                          unsigned int            numRows,
                          unsigned int            numExtraCols,
                    const std::vector<bool>*      extraColIsString,
                    const std::string&            fileName);
 ~uqAsciiTableClass();

  unsigned int                             numRows     ()                 const;
  unsigned int                             numCols     ()                 const;
  const EpetraExt::DistArray<std::string>& stringColumn(unsigned int j)   const;
  const V&                                 doubleColumn(unsigned int j)   const;
  void                                     print       (std::ostream& os) const;

private:
  Epetra_Map* newMap(); // See template specialization

  const uqBaseEnvironmentClass&                   m_env;
  unsigned int                                    m_numRows;
  unsigned int                                    m_numCols;
  std::vector<bool>                               m_colIsString;
  std::string                                     m_fileName;

  const Epetra_Map*                               m_map;
  std::vector<EpetraExt::DistArray<std::string>*> m_stringColumns;
  std::vector<V*>                                 m_doubleColumns;

  void readColumnsFromFile();
};

template <class V, class M>
uqAsciiTableClass<V,M>::uqAsciiTableClass(
  const uqBaseEnvironmentClass& env,
        unsigned int            numRows,
        unsigned int            numExtraCols,
  const std::vector<bool>*      extraColIsString,
  const std::string&            fileName)
  :
  m_env          (env),
  m_numRows      (numRows),
  m_numCols      (1+numExtraCols),
  m_colIsString  (1,true),
  m_fileName     (fileName),
  m_map          (newMap()),
  m_stringColumns(0),
  m_doubleColumns(0)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqAsciiTableClass<V,M>::constructor()..."
                            << std::endl;
  }

  if (m_numCols > 1) {
    m_colIsString.resize(m_numCols,false);
    if (extraColIsString == NULL) {
      // Nothing extra needs to be done
    }
    else {
      unsigned int maxJ = std::min(numExtraCols,(unsigned int) extraColIsString->size());
      for (unsigned int j = 0; j < maxJ; ++j) {
        m_colIsString[1+j] = (*extraColIsString)[j];
      }
    }
  }
  m_stringColumns.resize(m_numCols,NULL);
  m_doubleColumns.resize(m_numCols,NULL);
  readColumnsFromFile();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqAsciiTableClass<V,M>::constructor()"
                            << std::endl;
  }
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
uqAsciiTableClass<V,M>::readColumnsFromFile()
{
  unsigned int maxCharsPerLine = 512;

  std::ifstream ifs(m_fileName.c_str());
  UQ_FATAL_TEST_MACRO(ifs.is_open() == false,
                      m_env.worldRank(),
                      "uqAsciiTableClass<V,M>::readColumnsFromFile()",
                      "file was not found");

  // Determine number of lines
  unsigned int numLines = std::count(std::istreambuf_iterator<char>(ifs),
                                     std::istreambuf_iterator<char>(),
                                     '\n');

  // Determine number of valid lines
  int iRC;
  ifs.seekg(0,std::ios_base::beg);
  unsigned int lineId = 0;
  unsigned int numValidLines = 0;
  std::string tempString;
  while ((lineId < numLines) && (ifs.eof() == false)) {
    iRC = uqMiscReadStringAndDoubleFromFile(ifs,tempString,NULL);
    UQ_FATAL_TEST_MACRO(iRC,
                        m_env.worldRank(),
                        "uqAsciiTableClass<V,M>::readColumnsFromFile()",
                        "failed reading during the determination of the number of valid lines");
    //*m_env.subDisplayFile() << "lineId = "          << lineId
    //                        << ", numValidLines = " << numValidLines
    //                        << ", tempString = "    << tempString
    //                        << std::endl;
    if (tempString[0] != '#') numValidLines++;
    lineId++;
    ifs.ignore(maxCharsPerLine,'\n');
  }
  UQ_FATAL_TEST_MACRO(lineId != numLines,
                      m_env.worldRank(),
                      "uqAsciiTableClass<V,M>::readColumnsFromFile()",
                      "the first number of lines read is nonconsistent");
  if (m_numRows != numValidLines) {
    char errorExplanation[512];
    sprintf(errorExplanation,"number of valid lines (%u) in ascii table file does not match number of rows (%u)",numValidLines,m_numRows);
    UQ_FATAL_TEST_MACRO(true,
                        m_env.worldRank(),
                        "uqAsciiTableClass<V,M>::readColumnsFromFile()",
                        errorExplanation);
  }

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Ascii table file '"    << m_fileName
                            << "' has "                << numLines
                            << " lines and specifies " << numValidLines
                            << " valid lines."
                            << std::endl;
  }

  for (unsigned int j=0; j < m_numCols; ++j) {
    if (m_colIsString[j]) {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
        *m_env.subDisplayFile() << "Column j = " << j
                                << " is a columns of strings"
                                << std::endl;
      }
      m_stringColumns[j] = new EpetraExt::DistArray<std::string>(*m_map,1);
    }
    else {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
        *m_env.subDisplayFile() << "Column j = " << j
                                << " is a columns of doubles"
                                << std::endl;
      }
      m_doubleColumns[j] = new V(m_env,*m_map);
    }
  }

  // Read file until End Of File character is reached
  ifs.seekg(0,std::ios_base::beg);
  lineId = 0;
  unsigned int validLineId = 0;
  std::string tmpString;
  while ((lineId < numLines) && (ifs.eof() == false)) {
    //*m_env.subDisplayFile() << "Beginning read of line (in ascii table file) of id = " << lineId << std::endl;
    bool endOfLineAchieved = false;

    iRC = uqMiscReadCharsAndDoubleFromFile(ifs, tmpString, NULL, endOfLineAchieved);
    UQ_FATAL_TEST_MACRO(iRC,
                        m_env.worldRank(),
                        "uqAsciiTableClass<V,M>::readColumnsFromFile()",
                        "failed reading a first column during the valid lines reading loop");

    lineId++;
    if (tmpString[0] == '#') {
      if (!endOfLineAchieved) ifs.ignore(maxCharsPerLine,'\n');
      continue;
    }

    EpetraExt::DistArray<std::string>& firstColumn = *m_stringColumns[0];
    firstColumn(validLineId,0) = tmpString;

    // Check 'validLineId' before setting one more valid line
    if (validLineId >= numValidLines) {
      char errorExplanation[512];
      sprintf(errorExplanation,"validLineId (%u) got too large during reading of ascii table file",validLineId);
      UQ_FATAL_TEST_MACRO(true,
                          m_env.worldRank(),
                          "uqAsciiTableClass<V,M>::readColumnsFromFile()",
                          errorExplanation);
    }

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
      *m_env.subDisplayFile() << "Just read a string: table[" << validLineId
                              << ","                          << 0 // j=0
                              << "] = "                       << firstColumn(validLineId,0)
                              << std::endl;
    }

    for (unsigned int j=1; j < m_numCols; ++j) {
      UQ_FATAL_TEST_MACRO(endOfLineAchieved,
                          m_env.worldRank(),
                          "uqAsciiTableClass<V,M>::readColumnsFromFile()",
                          "failed reading all columns in a valid line");
      if (m_colIsString[j]) {
        EpetraExt::DistArray<std::string>& arrayOfStrings = *m_stringColumns[j];
        iRC = uqMiscReadCharsAndDoubleFromFile(ifs, arrayOfStrings(validLineId,0), NULL, endOfLineAchieved);
        UQ_FATAL_TEST_MACRO(iRC,
                            m_env.worldRank(),
                            "uqAsciiTableClass<V,M>::readColumnsFromFile()",
                            "failed reading a string column in a valid line");
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
          *m_env.subDisplayFile() << "Just read a string: table[" << validLineId
                                  << ","                          << j
                                  << "] = "                       << arrayOfStrings(validLineId,0)
                                  << std::endl;
        }
      }
      else {
        iRC = uqMiscReadCharsAndDoubleFromFile(ifs, tmpString, &(*m_doubleColumns[j])[validLineId], endOfLineAchieved);
        UQ_FATAL_TEST_MACRO(iRC,
                            m_env.worldRank(),
                            "uqAsciiTableClass<V,M>::readColumnsFromFile()",
                            "failed reading a double column in a valid line");
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
          *m_env.subDisplayFile() << "Just read a double: table[" << validLineId
                                  << ","                          << j
                                  << "] = "                       << (*m_doubleColumns[j])[validLineId]
                                  << std::endl;
        }
      }
    }
    if (!endOfLineAchieved) ifs.ignore(maxCharsPerLine,'\n');

    validLineId++;
  }

  UQ_FATAL_TEST_MACRO(lineId != numLines,
                      m_env.worldRank(),
                      "uqAsciiTableClass<V,M>::readColumnsFromFile()",
                      "the second number of lines read is not consistent");
  UQ_FATAL_TEST_MACRO(validLineId != numValidLines,
                      m_env.worldRank(),
                      "uqAsciiTableClass<V,M>::readColumnsFromFile()",
                      "the number of valid lines just read is not consistent");

  if (m_env.displayVerbosity() >= 5) {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "Finished reading table '" << m_fileName
                              << "'. Its contents per column are:"
                              << std::endl;
      *m_env.subDisplayFile() << *this; // FIX ME: output might need to be in parallel
      *m_env.subDisplayFile() << std::endl;
    }    
  }

  return;
}

template <class V, class M>
unsigned int
uqAsciiTableClass<V,M>::numRows() const
{
  return m_numRows;
}

template <class V, class M>
unsigned int
uqAsciiTableClass<V,M>::numCols() const
{
  return m_numCols;
}

template <class V, class M>
const EpetraExt::DistArray<std::string>&
uqAsciiTableClass<V,M>::stringColumn(unsigned int j) const
{
  UQ_FATAL_TEST_MACRO(j >= m_numCols,
                      m_env.worldRank(),
                      "uqAsciiTableClass<V,M>::stringColumn()",
                      "invalid j");

  UQ_FATAL_TEST_MACRO(m_stringColumns[j] == NULL,
                      m_env.worldRank(),
                      "uqAsciiTableClass<V,M>::stringColumn()",
                      "string column is not ready");

  return *m_stringColumns[j];
}

template <class V, class M>
const V&
uqAsciiTableClass<V,M>::doubleColumn(unsigned int j) const
{
  UQ_FATAL_TEST_MACRO(j >= m_numCols,
                      m_env.worldRank(),
                      "uqAsciiTableClass<V,M>::doubleColumn()",
                      "invalid j");

  UQ_FATAL_TEST_MACRO(m_doubleColumns[j] == NULL,
                      m_env.worldRank(),
                      "uqAsciiTableClass<V,M>::doubleColumn()",
                      "double column is not ready");

  return *m_doubleColumns[j];
}

template <class V, class M>
void
uqAsciiTableClass<V,M>::print(std::ostream& os) const
{
  for (unsigned int j = 0; j < m_numCols; ++j) {
    UQ_FATAL_TEST_MACRO((m_stringColumns[j] != NULL) && (m_doubleColumns[j] != NULL),
                        m_env.worldRank(),
                        "uqAsciiTableClass<V,M>::print()",
                        "column is not null on both possible ways");
    UQ_FATAL_TEST_MACRO((m_stringColumns[j] == NULL) && (m_doubleColumns[j] == NULL),
                        m_env.worldRank(),
                        "uqAsciiTableClass<V,M>::print()",
                        "column is null on both possible ways");

    os << "\nContents of table '" << m_fileName
       << "', column "            << j
       << ":"
       << std::endl;
    if (m_stringColumns[j] != NULL) {
      os << *m_stringColumns[j];
      //EpetraExt::DistArray<std::string>& arrayOfStrings = *m_stringColumns[j];
      //for (unsigned int i = 0; i < m_numRows; ++i) {
      //  os << arrayOfStrings(i,0)
      //     << std::endl;
      //}
    }
    else {
      os << *m_doubleColumns[j];
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

