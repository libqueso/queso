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
// $Id:$
//
//--------------------------------------------------------------------------

#ifndef __UQ_MULTI_LEVEL_SAMPLING2_H__
#define __UQ_MULTI_LEVEL_SAMPLING2_H__

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::checkpointML(
  double                                   currExponent,            // input
  double                                   currEta,                 // input
  const uqSequenceOfVectorsClass<P_V,P_M>& currChain,               // input
  const uqScalarSequenceClass<double>&     currLogLikelihoodValues, // input
  const uqScalarSequenceClass<double>&     currLogTargetValues)     // input
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "\n CHECKPOINTING initiating at level " << m_currLevel
                            << "\n" << std::endl;
  }

  //******************************************************************************
  // Write 'control' file without 'level' spefication in name
  //******************************************************************************
  unsigned int quantity1 = currChain.unifiedSequenceSize();
  unsigned int quantity2 = currLogLikelihoodValues.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
  unsigned int quantity3 = currLogTargetValues.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
  if (m_env.inter0Rank() >= 0) {
    UQ_FATAL_TEST_MACRO(m_logEvidenceFactors.size() != m_currLevel,
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::checkpointML()",
                        "number of evidence factors is not consistent");
    UQ_FATAL_TEST_MACRO(quantity1 != quantity2,
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::checkpointML()",
                        "quantity2 is not consistent");
    UQ_FATAL_TEST_MACRO(quantity1 != quantity3,
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::checkpointML()",
                        "quantity3 is not consistent");
  }

  if (m_env.fullRank() == 0) {
    std::ofstream* ofsVar = new std::ofstream((m_options.m_restartOutput_baseNameForFiles + "Control.txt").c_str(), 
                                              std::ofstream::out | std::ofstream::trunc);
    *ofsVar << m_currLevel               << std::endl  // 1
            << m_vectorSpace.dimGlobal() << std::endl  // 2
            << currExponent              << std::endl  // 3
            << currEta                   << std::endl  // 4
            << quantity1                 << std::endl; // 5
    unsigned int savedPrecision = ofsVar->precision();
    ofsVar->precision(16);
    for (unsigned int i = 0; i < m_logEvidenceFactors.size(); ++i) {
      *ofsVar << m_logEvidenceFactors[i] << std::endl;
    }
    ofsVar->precision(savedPrecision);
    *ofsVar << "COMPLETE"                << std::endl; // 6 = ML_CHECKPOINT_FIXED_AMOUNT_OF_DATA

    delete ofsVar;
  }
  m_env.fullComm().Barrier();

  //******************************************************************************
  // Write three 'data' files
  //******************************************************************************
  char levelSufix[256];
  sprintf(levelSufix,"%d",m_currLevel+LEVEL_REF_ID); // Yes, '+0'

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "\n CHECKPOINTING chain at level " << m_currLevel
                            << "\n" << std::endl;
  }
  currChain.unifiedWriteContents(m_options.m_restartOutput_baseNameForFiles + "Chain_l" + levelSufix,
                                 m_options.m_restartOutput_fileType);
  m_env.fullComm().Barrier();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "\n CHECKPOINTING like at level " << m_currLevel
                            << "\n" << std::endl;
  }
  currLogLikelihoodValues.unifiedWriteContents(m_options.m_restartOutput_baseNameForFiles + "LogLike_l" + levelSufix,
                                               m_options.m_restartOutput_fileType);
  m_env.fullComm().Barrier();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "\n CHECKPOINTING target at level " << m_currLevel
                            << "\n" << std::endl;
  }
  currLogTargetValues.unifiedWriteContents(m_options.m_restartOutput_baseNameForFiles + "LogTarget_l" + levelSufix,
                                           m_options.m_restartOutput_fileType);
  m_env.fullComm().Barrier();

  //******************************************************************************
  // Write 'control' file *with* 'level' spefication in name
  //******************************************************************************
  if (m_env.fullRank() == 0) {
    std::ofstream* ofsVar = new std::ofstream((m_options.m_restartOutput_baseNameForFiles + "Control_l" + levelSufix + ".txt").c_str(), 
                                              std::ofstream::out | std::ofstream::trunc);
    *ofsVar << m_currLevel               << std::endl  // 1
            << m_vectorSpace.dimGlobal() << std::endl  // 2
            << currExponent              << std::endl  // 3
            << currEta                   << std::endl  // 4
            << quantity1                 << std::endl; // 5
    unsigned int savedPrecision = ofsVar->precision();
    ofsVar->precision(16);
    for (unsigned int i = 0; i < m_logEvidenceFactors.size(); ++i) {
      *ofsVar << m_logEvidenceFactors[i] << std::endl;
    }
    ofsVar->precision(savedPrecision);
    *ofsVar << "COMPLETE"                << std::endl; // 6 = ML_CHECKPOINT_FIXED_AMOUNT_OF_DATA

    delete ofsVar;
  }
  m_env.fullComm().Barrier();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "\n CHECKPOINTING done at level " << m_currLevel
                            << "\n" << std::endl;
  }

  return;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::restartML(
  double&                            currExponent,            // output
  double&                            currEta,                 // output
  uqSequenceOfVectorsClass<P_V,P_M>& currChain,               // output
  uqScalarSequenceClass<double>&     currLogLikelihoodValues, // output
  uqScalarSequenceClass<double>&     currLogTargetValues)     // output
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "\n RESTARTING initiating at level " << m_currLevel
                            << "\n" << std::endl;
  }

  //******************************************************************************
  // Read 'control' file
  //******************************************************************************
  unsigned int vectorSpaceDim  = 0;
  unsigned int quantity1       = 0;
  std::string  checkingString("");
  if (m_env.fullRank() == 0) {
    std::ifstream* ifsVar = new std::ifstream((m_options.m_restartInput_baseNameForFiles + "Control.txt").c_str(), 
                                              std::ifstream::in);

    //******************************************************************************
    // Determine number of lines
    //******************************************************************************
    unsigned int numLines = std::count(std::istreambuf_iterator<char>(*ifsVar),
                                       std::istreambuf_iterator<char>(),
                                       '\n');
    ifsVar->seekg(0,std::ios_base::beg);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
      *m_env.subDisplayFile() << "Restart input file has " << numLines
                              << " lines"
                              << std::endl;
    }

    //******************************************************************************
    // Read all values
    //******************************************************************************
    *ifsVar >> m_currLevel; // 1
    UQ_FATAL_TEST_MACRO(numLines != (ML_CHECKPOINT_FIXED_AMOUNT_OF_DATA + m_currLevel),
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::restartML()",
                        "number of lines read is different than pre-established number of lines in control file");

    m_logEvidenceFactors.clear();
    m_logEvidenceFactors.resize(m_currLevel,0.);
    *ifsVar >> vectorSpaceDim  // 2
            >> currExponent    // 3
            >> currEta         // 4
            >> quantity1;      // 5
    for (unsigned int i = 0; i < m_logEvidenceFactors.size(); ++i) {
      *ifsVar >> m_logEvidenceFactors[i];
    }
    *ifsVar >> checkingString; // 6 = ML_CHECKPOINT_FIXED_AMOUNT_OF_DATA
    UQ_FATAL_TEST_MACRO(checkingString != "COMPLETE",
                        m_env.fullRank(),
                        "uqMLSamplingClass<P_V,P_M>::restartML()",
                        "control txt input file is not complete");

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
      *m_env.subDisplayFile() << "Restart input file has the following information:"
                              << "\n m_currLevel = "      << m_currLevel
                              << "\n vectorSpaceDim = "   << vectorSpaceDim
                              << "\n currExponent = "     << currExponent
                              << "\n currEta = "          << currEta
                              << "\n quantity1 = "        << quantity1;
      for (unsigned int i = 0; i < m_logEvidenceFactors.size(); ++i) {
        *m_env.subDisplayFile() << "\n [" << i << "] = " << m_logEvidenceFactors[i];
      }
      *m_env.subDisplayFile() << std::endl;
    }

#if 0 // For debug only
    std::string tmpString;
    for (unsigned int i = 0; i < 2; ++i) {
      *ifsVar >> tmpString;
      std::cout << "Just read '" << tmpString << "'" << std::endl;
    }
    while ((lineId < numLines) && (ifsVar->eof() == false)) {
    }
    ifsVar->ignore(maxCharsPerLine,'\n');
#endif

    delete ifsVar;
  } // if (m_env.fullRank() == 0)
  m_env.fullComm().Barrier();

  //******************************************************************************
  // MPI_Bcast 'm_currLevel'
  //******************************************************************************
  unsigned int tmpUint = (unsigned int) m_currLevel;
  m_env.fullComm().Bcast((void *) &tmpUint, (int) 1, uqRawValue_MPI_UNSIGNED, 0, // Yes, 'fullComm'
                         "uqMLSamplingClass<P_V,P_M>::restartML()",
                         "failed MPI.Bcast() for m_currLevel");
  if (m_env.fullRank() != 0) {
    m_currLevel = tmpUint;
  }

  //******************************************************************************
  // MPI_Bcast the rest of the information just read
  //******************************************************************************
  std::vector<double> tmpData(ML_CHECKPOINT_FIXED_AMOUNT_OF_DATA-1+m_currLevel,0.);
  if (m_env.fullRank() == 0) {
    tmpData[0] = vectorSpaceDim;
    tmpData[1] = currExponent;
    tmpData[2] = currEta;
    tmpData[3] = quantity1;
    for (unsigned int i = 0; i < m_logEvidenceFactors.size(); ++i) {
      tmpData[4+i] = m_logEvidenceFactors[i];
    }
  }
  else {
    m_logEvidenceFactors.clear();
    m_logEvidenceFactors.resize(m_currLevel,0.);
  }
  m_env.fullComm().Bcast((void *) &tmpData[0], (int) tmpData.size(), uqRawValue_MPI_DOUBLE, 0, // Yes, 'fullComm'
                         "uqMLSamplingClass<P_V,P_M>::restartML()",
                         "failed MPI.Bcast() for rest of information read from input file");
  if (m_env.fullRank() != 0) {
    vectorSpaceDim = tmpData[0];
    currExponent   = tmpData[1];
    currEta        = tmpData[2];
    quantity1      = tmpData[3];
    for (unsigned int i = 0; i < m_logEvidenceFactors.size(); ++i) {
      m_logEvidenceFactors[i] = tmpData[4+i];
    }
  }

  //******************************************************************************
  // Process read data in all MPI nodes now
  //******************************************************************************
  UQ_FATAL_TEST_MACRO(vectorSpaceDim != m_vectorSpace.dimGlobal(),
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::restartML()",
                      "read vector space dimension is not consistent");
  UQ_FATAL_TEST_MACRO((currExponent < 0.) || (currExponent > 1.),
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::restartML()",
                      "read currExponent is not consistent");
  UQ_FATAL_TEST_MACRO((quantity1 % m_env.numSubEnvironments()) != 0,
                      m_env.fullRank(),
                      "uqMLSamplingClass<P_V,P_M>::restartML()",
                      "read size of chain should be a multiple of the number of subenvironments");
  unsigned int subSequenceSize = 0;
  subSequenceSize = ((double) quantity1) / ((double) m_env.numSubEnvironments());

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Restart input file has the following information"
                            << ": subSequenceSize = " << subSequenceSize
                            << std::endl;
  }

  //******************************************************************************
  // Read three 'data' files
  //******************************************************************************
  char levelSufix[256];
  sprintf(levelSufix,"%d",m_currLevel+LEVEL_REF_ID); // Yes, '+0'

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "\n RESTARTING chain at level " << m_currLevel
                            << "\n" << std::endl;
  }
  currChain.unifiedReadContents(m_options.m_restartInput_baseNameForFiles + "Chain_l" + levelSufix,
                                m_options.m_restartInput_fileType,
                                subSequenceSize);
  m_env.fullComm().Barrier();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "\n RESTARTING like at level " << m_currLevel
                            << "\n" << std::endl;
  }
  currLogLikelihoodValues.unifiedReadContents(m_options.m_restartInput_baseNameForFiles + "LogLike_l" + levelSufix,
                                              m_options.m_restartInput_fileType,
                                              subSequenceSize);
  m_env.fullComm().Barrier();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "\n RESTARTING target at level " << m_currLevel
                            << "\n" << std::endl;
  }
  currLogTargetValues.unifiedReadContents(m_options.m_restartInput_baseNameForFiles + "LogTarget_l" + levelSufix,
                                          m_options.m_restartInput_fileType,
                                          subSequenceSize);
  m_env.fullComm().Barrier();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "\n RESTARTING done at level " << m_currLevel
                            << "\n" << std::endl;
  }

  return;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::generateSequence_Level0_all(
  const uqMLSamplingLevelOptionsClass& currOptions,                // input
  unsigned int&                        unifiedRequestedNumSamples, // output
  uqSequenceOfVectorsClass<P_V,P_M>&   currChain,                  // output
  uqScalarSequenceClass<double>&       currLogLikelihoodValues,    // output
  uqScalarSequenceClass<double>&       currLogTargetValues)        // output
{
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "KEY In uqMLSampling<P_V,P_M>::generateSequence()"
                              << ": beginning level "              << m_currLevel+LEVEL_REF_ID
                              << ", currOptions.m_rawChainSize = " << currOptions.m_rawChainSize // Ok to use rawChainSize
                              << std::endl;
    }

    int iRC = UQ_OK_RC;
    struct timeval timevalLevel;
    iRC = gettimeofday(&timevalLevel, NULL);
    if (iRC) {}; // just to remove compiler warning

    if (m_env.inter0Rank() >= 0) {
      unsigned int tmpSize = currOptions.m_rawChainSize;
      m_env.inter0Comm().Allreduce((void *) &tmpSize, (void *) &unifiedRequestedNumSamples, (int) 1, uqRawValue_MPI_UNSIGNED, uqRawValue_MPI_SUM,
                                   "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                                   "failed MPI.Allreduce() for requested num samples in level 0");
    }
    else {
      unifiedRequestedNumSamples = currOptions.m_rawChainSize;
    }

    currChain.setName              (currOptions.m_prefix + "rawChain"        );
    currLogLikelihoodValues.setName(currOptions.m_prefix + "rawLogLikelihood");
    currLogTargetValues.setName    (currOptions.m_prefix + "rawLogTarget"    );

    currChain.resizeSequence              (currOptions.m_rawChainSize); // Ok to use rawChainSize
    currLogLikelihoodValues.resizeSequence(currOptions.m_rawChainSize); // Ok to use rawChainSize
    currLogTargetValues.resizeSequence    (currOptions.m_rawChainSize); // Ok to use rawChainSize

    P_V auxVec(m_vectorSpace.zeroVector());
    uqScalarFunctionSynchronizerClass<P_V,P_M> likelihoodSynchronizer(m_likelihoodFunction,auxVec); // prudencio 2010-08-01
    for (unsigned int i = 0; i < currChain.subSequenceSize(); ++i) {
      //std::cout << "In QUESO: before prior realizer with i = " << i << std::endl;
      bool outOfSupport = true;
      do {

        m_priorRv.realizer().realization(auxVec);
        auxVec.mpiBcast(0, m_env.subComm()); // prudencio 2010-08-01

        outOfSupport = !(m_targetDomain->contains(auxVec));
      } while (outOfSupport); // prudenci 2011-Oct-04
      
      currChain.setPositionValues(i,auxVec);
      // KAUST: all nodes should call likelihood
#if 1 // prudencio 2010-08-01
      currLogLikelihoodValues[i] = likelihoodSynchronizer.callFunction(&auxVec,NULL,NULL,NULL,NULL,NULL,NULL); // likelihood is important
#else
      currLogLikelihoodValues[i] = m_likelihoodFunction.lnValue(auxVec,NULL,NULL,NULL,NULL); // likelihood is important
#endif
      currLogTargetValues[i]     = m_priorRv.pdf().lnValue(auxVec,NULL,NULL,NULL,NULL) + currLogLikelihoodValues[i];
      //std::cout << "In QUESO: currLogTargetValues[" << i << "] = " << currLogTargetValues[i] << std::endl;
    }

    if (m_env.inter0Rank() >= 0) { // KAUST
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
      if (currOptions.m_rawChainComputeStats) {
        uqFilePtrSetStruct filePtrSet;
        m_env.openOutputFile(currOptions.m_dataOutputFileName,
                             UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT, // Yes, always ".m"
                             currOptions.m_dataOutputAllowedSet,
                             false,
                             filePtrSet);

        //m_env.syncPrintDebugMsg("At level 0, calling computeStatistics for chain",1,10,m_env.inter0Comm()); // output debug
        currChain.computeStatistics(*currOptions.m_rawChainStatisticalOptionsObj,
                                    filePtrSet.ofsVar);

        m_env.closeFile(filePtrSet,UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT);
      }
      // Compute MLE and MAP
      // rr0
#endif
      if (currOptions.m_rawChainDataOutputFileName != UQ_MH_SG_FILENAME_FOR_NO_FILE) {
        currChain.unifiedWriteContents              (currOptions.m_rawChainDataOutputFileName,
                                                     currOptions.m_rawChainDataOutputFileType);
        currLogLikelihoodValues.unifiedWriteContents(currOptions.m_rawChainDataOutputFileName,
                                                     currOptions.m_rawChainDataOutputFileType);
        currLogTargetValues.unifiedWriteContents    (currOptions.m_rawChainDataOutputFileName,
                                                     currOptions.m_rawChainDataOutputFileType);
      }

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ": finished generating " << currChain.subSequenceSize()
                                << " chain positions"
                                << std::endl;

        //unsigned int numZeros = 0;
        //for (unsigned int i = 0; i < currTargetValues.subSequenceSize(); ++i) {
        //  *m_env.subDisplayFile() << "currTargetValues[" << i
        //                          << "] = " << currTargetValues[i]
        //                          << std::endl;
        //  if (currTargetValues[i] == 0.) numZeros++;
        //}
        //*m_env.subDisplayFile() << "Number of zeros in currTargetValues = " << numZeros
        //                        << std::endl;
      }

      if (currOptions.m_filteredChainGenerate) {
        // todo
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
        if (currOptions.m_filteredChainComputeStats) {
          // todo

          //currChain.computeStatistics(*currOptions.m_filteredChainStatisticalOptionsObj,
          //                            filePtrSet.ofsVar);
        }
        // Compute MLE and MAP
        // rr0
#endif
      }

    } // KAUST

    UQ_FATAL_TEST_MACRO((currChain.subSequenceSize() != currOptions.m_rawChainSize), // Ok to use rawChainSize
                        m_env.worldRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                        "currChain (first one) has been generated with invalid size");

    double levelRunTime = uqMiscGetEllapsedSeconds(&timevalLevel);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                              << ": ending level " << m_currLevel+LEVEL_REF_ID
                              << ", total level time = " << levelRunTime << " seconds"
                              << std::endl;
    }

  return;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::generateSequence_Step01_inter0(
  const uqMLSamplingLevelOptionsClass* currOptions,                // input
  unsigned int&                        unifiedRequestedNumSamples) // output
{
  int iRC = UQ_OK_RC;
  struct timeval timevalStep;
  iRC = gettimeofday(&timevalStep, NULL);
  if (iRC) {}; // just to remove compiler warning

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": beginning step 1 of 11"
                                << std::endl;
      }

      unsigned int tmpSize = currOptions->m_rawChainSize;
      // This computed 'unifiedRequestedNumSamples' needs to be recomputed only at the last
      // level, when 'currOptions' is replaced by 'lastLevelOptions' (see step 3 of 11)
      m_env.inter0Comm().Allreduce((void *) &tmpSize, (void *) &unifiedRequestedNumSamples, (int) 1, uqRawValue_MPI_UNSIGNED, uqRawValue_MPI_SUM,
                                   "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                                   "failed MPI.Allreduce() for requested num samples in step 1");

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "KEY In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ", currOptions->m_rawChainSize = " << currOptions->m_rawChainSize // Ok to use rawChainSize
                                << std::endl;
      }

  double stepRunTime = uqMiscGetEllapsedSeconds(&timevalStep);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving uqMLSampling<P_V,P_M>::generateSequence_Step()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", after " << stepRunTime << " seconds"
                            << std::endl;
  }

  return;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::generateSequence_Step02_inter0(
  const uqMLSamplingLevelOptionsClass* currOptions,             // input
  uqSequenceOfVectorsClass<P_V,P_M>&   currChain,               // input/output
  uqScalarSequenceClass<double>&       currLogLikelihoodValues, // input/output
  uqScalarSequenceClass<double>&       currLogTargetValues,     // input/output
  uqSequenceOfVectorsClass<P_V,P_M>&   prevChain,               // output
  uqScalarSequenceClass<double>&       prevLogLikelihoodValues, // output
  uqScalarSequenceClass<double>&       prevLogTargetValues,     // output
  unsigned int&                        indexOfFirstWeight,      // output
  unsigned int&                        indexOfLastWeight)       // output
{
  int iRC = UQ_OK_RC;
  struct timeval timevalStep;
  iRC = gettimeofday(&timevalStep, NULL);
  if (iRC) {}; // just to remove compiler warning

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": beginning step 2 of 11"
                                << std::endl;
      }

      prevChain = currChain;
      currChain.clear();
      currChain.setName(currOptions->m_prefix + "rawChain");

      prevLogLikelihoodValues = currLogLikelihoodValues; // likelihood is important
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence_Step()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ", prevLogLikelihoodValues[0] = " << prevLogLikelihoodValues[0]
                                << std::endl;
      }
      prevLogTargetValues     = currLogTargetValues;

      currLogLikelihoodValues.clear();
      currLogLikelihoodValues.setName(currOptions->m_prefix + "rawLogLikelihood");

      currLogTargetValues.clear();
      currLogTargetValues.setName(currOptions->m_prefix + "rawLogTarget");

#if 0 // For debug only
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        P_V prevPosition(m_vectorSpace.zeroVector());
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ":"
                                << std::endl;
        for (unsigned int i = 0; i < prevChain.subSequenceSize(); ++i) {
          prevChain.getPositionValues(i,prevPosition);
          *m_env.subDisplayFile() << "  prevChain[" << i
                                  << "] = " << prevPosition
                                  << ", prevLogLikelihoodValues[" << i
                                  << "] = " << prevLogLikelihoodValues[i]
                                  << ", prevLogTargetValues[" << i
                                  << "] = " << prevLogTargetValues[i]
                                  << std::endl;
        }
      }
#endif

      unsigned int quantity1 = prevChain.unifiedSequenceSize();
      unsigned int quantity2 = currChain.unifiedSequenceSize();
      unsigned int quantity3 = prevLogLikelihoodValues.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
      unsigned int quantity4 = currLogLikelihoodValues.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
      unsigned int quantity5 = prevLogTargetValues.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
      unsigned int quantity6 = currLogTargetValues.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": prevChain.unifiedSequenceSize() = " << quantity1
                                << ", currChain.unifiedSequenceSize() = " << quantity2
                                << ", prevLogLikelihoodValues.unifiedSequenceSize() = " << quantity3
                                << ", currLogLikelihoodValues.unifiedSequenceSize() = " << quantity4
                                << ", prevLogTargetValues.unifiedSequenceSize() = " << quantity5
                                << ", currLogTargetValues.unifiedSequenceSize() = " << quantity6
                                << std::endl;
      }

      UQ_FATAL_TEST_MACRO((prevChain.subSequenceSize() != prevLogLikelihoodValues.subSequenceSize()),
                          m_env.worldRank(),
                          "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                          "different sizes between previous chain and previous sequence of likelihood values");

      UQ_FATAL_TEST_MACRO((prevChain.subSequenceSize() != prevLogTargetValues.subSequenceSize()),
                          m_env.worldRank(),
                          "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                          "different sizes between previous chain and previous sequence of target values");

      // Set 'indexOfFirstWeight' and 'indexOfLastWeight' // KAUST
      indexOfFirstWeight = 0;
      indexOfLastWeight  = indexOfFirstWeight + prevChain.subSequenceSize()-1;
      {
        //std::cout << "m_env.inter0Comm().NumProc() = " << m_env.inter0Comm().NumProc() << std::endl;
        int r = m_env.inter0Rank();
        //std::cout << "r = " << r << std::endl;
        m_env.inter0Comm().Barrier();
        unsigned int auxUint = 0;
        if (r > 0) {
          uqRawType_MPI_Status status;
	  //std::cout << "Rank " << r << " is entering MPI_Recv()" << std::endl;
          m_env.inter0Comm().Recv((void*) &auxUint, 1, uqRawValue_MPI_UNSIGNED, r-1, r-1, &status,
                                  "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                                  "failed MPI.Recv()");
	  //std::cout << "Rank " << r << " received auxUint = " << auxUint << std::endl;
          indexOfFirstWeight = auxUint;
          indexOfLastWeight = indexOfFirstWeight + prevChain.subSequenceSize()-1;
        }
        if (r < (m_env.inter0Comm().NumProc()-1)) {
          auxUint = indexOfLastWeight + 1;
	  //std::cout << "Rank " << r << " is sending auxUint = " << auxUint << std::endl;
          m_env.inter0Comm().Send((void*) &auxUint, 1, uqRawValue_MPI_UNSIGNED, r+1, r,
                                  "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                                  "failed MPI.Send()");
	  //std::cout << "Rank " << r << " sent auxUint = " << auxUint << std::endl;
        }
        m_env.inter0Comm().Barrier();
      }

  double stepRunTime = uqMiscGetEllapsedSeconds(&timevalStep);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving uqMLSampling<P_V,P_M>::generateSequence_Step()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", after " << stepRunTime << " seconds"
                            << std::endl;
  }

  return;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::generateSequence_Step03_inter0(
  const uqMLSamplingLevelOptionsClass* currOptions,             // input
  const uqScalarSequenceClass<double>& prevLogLikelihoodValues, // input
  double                               prevExponent,            // input
  double                               failedExponent,          // input // gpmsa
  double&                              currExponent,            // output
  uqScalarSequenceClass<double>&       weightSequence)          // output
{
  int iRC = UQ_OK_RC;
  struct timeval timevalStep;
  iRC = gettimeofday(&timevalStep, NULL);
  if (iRC) {}; // just to remove compiler warning

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ", failedExponent = " << failedExponent // gpmsa
                                << ": beginning step 3 of 11"
                                << std::endl;
      }

      std::vector<double> exponents(2,0.);
      exponents[0] = prevExponent;
      exponents[1] = 1.;

      double nowExponent = 1.; // Try '1.' right away
      double nowEffectiveSizeRatio = 0.; // To be computed

      unsigned int nowAttempt = 0;
      bool testResult = false;
      double meanEffectiveSizeRatio = .5*(currOptions->m_minEffectiveSizeRatio + currOptions->m_maxEffectiveSizeRatio);
      uqScalarSequenceClass<double> omegaLnDiffSequence(m_env,prevLogLikelihoodValues.subSequenceSize(),"");

      double nowUnifiedEvidenceLnFactor = 0.;
      do {
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                  << ", level " << m_currLevel+LEVEL_REF_ID
                                  << ", step "  << m_currStep
                                  << ", failedExponent = " << failedExponent // gpmsa
                                  << ": entering loop for computing next exponent"
                                  << ", with nowAttempt = " << nowAttempt
                                  << std::endl;
        }

        if (failedExponent > 0.) { // gpmsa
          nowExponent = .5*(prevExponent+failedExponent);
        }
        else {
          if (nowAttempt > 0) {
            if (nowEffectiveSizeRatio > meanEffectiveSizeRatio) {
              exponents[0] = nowExponent;
            }
            else {
              exponents[1] = nowExponent;
            }
            nowExponent = .5*(exponents[0] + exponents[1]);
          }
        }
        double auxExponent = nowExponent;
        if (prevExponent != 0.) {
          auxExponent /= prevExponent;
          auxExponent -= 1.;
        }
        double subWeightRatioSum     = 0.;
        double unifiedWeightRatioSum = 0.;

        for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
          omegaLnDiffSequence[i] = prevLogLikelihoodValues[i]*auxExponent; // likelihood is important
        }

#if 1 // prudenci-2012-07-06
      //double unifiedOmegaLnMin = omegaLnDiffSequence.unifiedMinPlain(m_vectorSpace.numOfProcsForStorage() == 1);
        double unifiedOmegaLnMax = omegaLnDiffSequence.unifiedMaxPlain(m_vectorSpace.numOfProcsForStorage() == 1);
#else
        double unifiedOmegaLnMin = 0.;
        double unifiedOmegaLnMax = 0.;
        omegaLnDiffSequence.unifiedMinMaxExtra(m_vectorSpace.numOfProcsForStorage() == 1, // KAUST3
                                               0,
                                               omegaLnDiffSequence.subSequenceSize(),
                                               unifiedOmegaLnMin,
                                               unifiedOmegaLnMax);
#endif
        for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
          omegaLnDiffSequence[i] -= unifiedOmegaLnMax;
          weightSequence[i] = exp(omegaLnDiffSequence[i]);
          subWeightRatioSum += weightSequence[i];
#if 0 // For debug only
          if ((m_currLevel == 1) && (nowAttempt == 6))  {
            if (m_env.subDisplayFile() && (m_env.displayVerbosity() >= 99)) {
              *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                      << ", level "                        << m_currLevel+LEVEL_REF_ID
                                      << ", step "                         << m_currStep
                                      << ", i = "                          << i
                                      << ", prevLogLikelihoodValues[i] = " << prevLogLikelihoodValues[i]
                                      << ", omegaLnDiffSequence[i] = "     << omegaLnDiffSequence[i]
                                      << ", weightSequence[i] = "          << weightSequence[i]
		//<< ", subWeightRatioSum = "          << subWeightRatioSum
                                      << std::endl;
            }
          }
#endif
        }
        m_env.inter0Comm().Allreduce((void *) &subWeightRatioSum, (void *) &unifiedWeightRatioSum, (int) 1, uqRawValue_MPI_DOUBLE, uqRawValue_MPI_SUM,
                                     "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                                     "failed MPI.Allreduce() for weight ratio sum");

        unsigned int auxQuantity = weightSequence.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
        nowUnifiedEvidenceLnFactor = log(unifiedWeightRatioSum) + unifiedOmegaLnMax - log(auxQuantity);

        double effectiveSampleSize = 0.;
        for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
          weightSequence[i] /= unifiedWeightRatioSum;
          effectiveSampleSize += weightSequence[i]*weightSequence[i];
          //if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          //  *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
          //                          << ", level "                 << m_currLevel+LEVEL_REF_ID
          //                          << ", step "                  << m_currStep
          //                          << ": i = "                   << i
          //                          << ", effectiveSampleSize = " << effectiveSampleSize
          //                          << std::endl;
          //}
        }

        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                  << ", level "                                  << m_currLevel+LEVEL_REF_ID
                                  << ", step "                                   << m_currStep
                                  << ": nowAttempt = "                           << nowAttempt
                                  << ", prevExponent = "                         << prevExponent
                                  << ", exponents[0] = "                         << exponents[0]
                                  << ", nowExponent = "                          << nowExponent
                                  << ", exponents[1] = "                         << exponents[1]
                                  << ", subWeightRatioSum = "                    << subWeightRatioSum
                                  << ", unifiedWeightRatioSum = "                << unifiedWeightRatioSum
                                  << ", unifiedOmegaLnMax = "                    << unifiedOmegaLnMax
                                  << ", weightSequence.unifiedSequenceSize() = " << auxQuantity
                                  << ", nowUnifiedEvidenceLnFactor = "           << nowUnifiedEvidenceLnFactor
                                  << ", effectiveSampleSize = "                  << effectiveSampleSize
                                  << std::endl;
        }

#if 0 // For debug only
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                  << ", level " << m_currLevel+LEVEL_REF_ID
                                  << ", step "  << m_currStep
                                  << ":"
                                  << std::endl;
        }
        for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
          if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
            *m_env.subDisplayFile() << "  weightSequence[" << i
                                    << "] = "              << weightSequence[i]
                                    << std::endl;
          }
        }
#endif

        double subQuantity = effectiveSampleSize;
        effectiveSampleSize = 0.;
        m_env.inter0Comm().Allreduce((void *) &subQuantity, (void *) &effectiveSampleSize, (int) 1, uqRawValue_MPI_DOUBLE, uqRawValue_MPI_SUM,
                                     "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                                     "failed MPI.Allreduce() for effective sample size");

        effectiveSampleSize = 1./effectiveSampleSize;
        nowEffectiveSizeRatio = effectiveSampleSize/((double) weightSequence.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1));
        UQ_FATAL_TEST_MACRO((nowEffectiveSizeRatio > (1.+1.e-8)),
                            m_env.worldRank(),
                            "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                            "effective sample size ratio cannot be > 1");
        //UQ_FATAL_TEST_MACRO((nowEffectiveSizeRatio < (1.-1.e-8)),
        //                    m_env.worldRank(),
        //                    "uqMLSamplingClass<P_V,P_M>::generateSequence()",
        //                    "effective sample size ratio cannot be < 1");

        if (failedExponent > 0.) { // gpmsa
          testResult = true;
        }
        else {
          //bool aux1 = (nowEffectiveSizeRatio == meanEffectiveSizeRatio);
          bool aux2 = (nowExponent == 1.                             )
                      &&
                      (nowEffectiveSizeRatio > meanEffectiveSizeRatio);
          bool aux3 = (nowEffectiveSizeRatio >= currOptions->m_minEffectiveSizeRatio)
                      &&
                      (nowEffectiveSizeRatio <= currOptions->m_maxEffectiveSizeRatio);
          testResult = aux2 || aux3;
        }

        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                  << ", level "                   << m_currLevel+LEVEL_REF_ID
                                  << ", step "                    << m_currStep
                                  << ": nowAttempt = "            << nowAttempt
                                  << ", prevExponent = "          << prevExponent
                                  << ", failedExponent = "        << failedExponent // gpmsa
                                  << ", exponents[0] = "          << exponents[0]
                                  << ", nowExponent = "           << nowExponent
                                  << ", exponents[1] = "          << exponents[1]
                                  << ", effectiveSampleSize = "   << effectiveSampleSize
                                  << ", weightSequenceSize = "    << weightSequence.subSequenceSize()
                                  << ", minEffectiveSizeRatio = " << currOptions->m_minEffectiveSizeRatio
                                  << ", nowEffectiveSizeRatio = " << nowEffectiveSizeRatio
                                  << ", maxEffectiveSizeRatio = " << currOptions->m_maxEffectiveSizeRatio
	    //<< ", aux2 = "                  << aux2
	    //<< ", aux3 = "                  << aux3
                                  << ", testResult = "            << testResult
                                  << std::endl;
        }
        nowAttempt++;

        // Make sure all nodes in 'inter0Comm' have the same value of 'nowExponent'
        if (uqMiscCheckForSameValueInAllNodes(nowExponent,
                                              0., // kept 'zero' on 2010/03/05
                                              m_env.inter0Comm(),
                                              "uqMLSamplingClass<P_V,P_M>::generateSequence(), step 3, nowExponent") == false) {
          if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
            *m_env.subDisplayFile() << "WARNING, In uqMLSampling<P_V,P_M>::generateSequence()"
                                    << ", level "        << m_currLevel+LEVEL_REF_ID
                                    << ", step "         << m_currStep
                                    << ": nowAttempt = " << nowAttempt
                                    << ", uqMiscCheck for 'nowExponent' detected a problem"
                                    << std::endl;
          }
        }

        // Make sure all nodes in 'inter0Comm' have the same value of 'testResult'
        if (uqMiscCheckForSameValueInAllNodes(testResult,
                                              0., // kept 'zero' on 2010/03/05
                                              m_env.inter0Comm(),
                                              "uqMLSamplingClass<P_V,P_M>::generateSequence(), step 3, testResult") == false) {
          if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
            *m_env.subDisplayFile() << "WARNING, In uqMLSampling<P_V,P_M>::generateSequence()"
                                    << ", level "        << m_currLevel+LEVEL_REF_ID
                                    << ", step "         << m_currStep
                                    << ": nowAttempt = " << nowAttempt
                                    << ", uqMiscCheck for 'testResult' detected a problem"
                                    << std::endl;
          }
        }
      } while (testResult == false);
      currExponent = nowExponent;
      if (failedExponent > 0.) { // gpmsa
        m_logEvidenceFactors[m_logEvidenceFactors.size()-1] = nowUnifiedEvidenceLnFactor;
      }
      else {
        m_logEvidenceFactors.push_back(nowUnifiedEvidenceLnFactor); // restart
      }

      unsigned int quantity1 = weightSequence.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level "                                  << m_currLevel+LEVEL_REF_ID
                                << ", step "                                   << m_currStep
                                << ": weightSequence.subSequenceSize() = "     << weightSequence.subSequenceSize()
                                << ", weightSequence.unifiedSequenceSize() = " << quantity1
                                << ", failedExponent = "                       << failedExponent // gpmsa
                                << ", currExponent = "                         << currExponent
                                << ", effective ratio = "                      << nowEffectiveSizeRatio
                                << ", log(evidence factor) = "                 << m_logEvidenceFactors[m_logEvidenceFactors.size()-1]
                                << ", evidence factor = "                      << exp(m_logEvidenceFactors[m_logEvidenceFactors.size()-1])
                                << std::endl;

        //unsigned int numZeros = 0;
        //for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
        //  *m_env.subDisplayFile() << "weightSequence[" << i
        //                          << "] = " << weightSequence[i]
        //                         << std::endl;
        //  if (weightSequence[i] == 0.) numZeros++;
        //}
        //*m_env.subDisplayFile() << "Number of zeros in weightSequence = " << numZeros
        //                        << std::endl;
      }

      // Make sure all nodes in 'inter0Comm' have the same value of 'logEvidenceFactor'
      if (uqMiscCheckForSameValueInAllNodes(m_logEvidenceFactors[m_logEvidenceFactors.size()-1],
                                            3.0e-16, // changed from 'zero' on 2010/03/03
                                            m_env.inter0Comm(),
                                            "uqMLSamplingClass<P_V,P_M>::generateSequence(), step 3, logEvidenceFactor") == false) {
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "WARNING, In uqMLSampling<P_V,P_M>::generateSequence()"
                                  << ", level "        << m_currLevel+LEVEL_REF_ID
                                  << ", step "         << m_currStep
                                  << ", failedExponent = " << failedExponent // gpmsa
                                  << ": nowAttempt = " << nowAttempt
                                  << ", uqMiscCheck for 'logEvidenceFactor' detected a problem"
                                  << std::endl;
        }
      }

  double stepRunTime = uqMiscGetEllapsedSeconds(&timevalStep);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving uqMLSampling<P_V,P_M>::generateSequence_Step()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", failedExponent = " << failedExponent // gpmsa
                            << ", after " << stepRunTime << " seconds"
                            << std::endl;
  }

  return;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::generateSequence_Step04_inter0(
  const uqSequenceOfVectorsClass<P_V,P_M>& prevChain,        // input
  const uqScalarSequenceClass<double>&     weightSequence,   // input
  P_M&                                     unifiedCovMatrix) // output
{
  int iRC = UQ_OK_RC;
  struct timeval timevalStep;
  iRC = gettimeofday(&timevalStep, NULL);
  if (iRC) {}; // just to remove compiler warning

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": beginning step 4 of 11"
                                << std::endl;
      }

      P_V auxVec(m_vectorSpace.zeroVector());
      P_V subWeightedMeanVec(m_vectorSpace.zeroVector());
      for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
        prevChain.getPositionValues(i,auxVec);
        subWeightedMeanVec += weightSequence[i]*auxVec;
      }

      // Todd Oliver 2010-09-07: compute weighted mean over all processors
      P_V unifiedWeightedMeanVec(m_vectorSpace.zeroVector());
      if (m_env.inter0Rank() >= 0) {
        subWeightedMeanVec.mpiAllReduce(uqRawValue_MPI_SUM,m_env.inter0Comm(),unifiedWeightedMeanVec);
      }
      else {
        unifiedWeightedMeanVec = subWeightedMeanVec;
      }

      P_V diffVec(m_vectorSpace.zeroVector());
      P_M subCovMatrix(m_vectorSpace.zeroVector());
      for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
        prevChain.getPositionValues(i,auxVec);
        diffVec = auxVec - unifiedWeightedMeanVec;
        subCovMatrix += weightSequence[i]*matrixProduct(diffVec,diffVec);
      }

      for (unsigned int i = 0; i < unifiedCovMatrix.numRowsLocal(); ++i) { // KAUST5
        for (unsigned int j = 0; j < unifiedCovMatrix.numCols(); ++j) {
          double localValue = subCovMatrix(i,j);
          double sumValue = 0.;
          if (m_env.inter0Rank() >= 0) {
            m_env.inter0Comm().Allreduce((void *) &localValue, (void *) &sumValue, (int) 1, uqRawValue_MPI_DOUBLE, uqRawValue_MPI_SUM,
                                         "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                                         "failed MPI.Allreduce() for cov matrix");
          }
          else {
            sumValue = localValue;
          }
          unifiedCovMatrix(i,j) = sumValue;
        }
      }

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level "              << m_currLevel+LEVEL_REF_ID
                                << ", step "               << m_currStep
                                << ": unifiedCovMatrix = " << unifiedCovMatrix
                                << std::endl;
      }

  double stepRunTime = uqMiscGetEllapsedSeconds(&timevalStep);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving uqMLSampling<P_V,P_M>::generateSequence_Step()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", after " << stepRunTime << " seconds"
                            << std::endl;
  }

  return;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::generateSequence_Step05_inter0(
  unsigned int                         unifiedRequestedNumSamples,        // input
  const uqScalarSequenceClass<double>& weightSequence,                    // input
  std::vector<unsigned int>&           unifiedIndexCountersAtProc0Only,   // output
  std::vector<double>&                 unifiedWeightStdVectorAtProc0Only) // output
{
  int iRC = UQ_OK_RC;
  struct timeval timevalStep;
  iRC = gettimeofday(&timevalStep, NULL);
  if (iRC) {}; // just to remove compiler warning

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": beginning step 5 of 11"
                                << std::endl;
      }

#if 0 // For debug only
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ", before weightSequence.getUnifiedContentsAtProc0Only()"
                                << ":"
                                << std::endl;
      }
      for (unsigned int i = 0; i < weightSequence.subSequenceSize(); ++i) {
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << ", weightSequence[" << i
                                  << "] = "              << weightSequence[i]
                                  << std::endl;
        }
      }
#endif

      weightSequence.getUnifiedContentsAtProc0Only(m_vectorSpace.numOfProcsForStorage() == 1,
                                                   unifiedWeightStdVectorAtProc0Only);

#if 0 // For debug only
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ", after weightSequence.getUnifiedContentsAtProc0Only()"
                                << ":"
                                << std::endl;
      }
      for (unsigned int i = 0; i < unifiedWeightStdVectorAtProc0Only.size(); ++i) {
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "  unifiedWeightStdVectorAtProc0Only[" << i
                                  << "] = "                                 << unifiedWeightStdVectorAtProc0Only[i]
                                  << std::endl;
        }
      }
#endif
      sampleIndexes_proc0(unifiedRequestedNumSamples,        // input
                          unifiedWeightStdVectorAtProc0Only, // input
                          unifiedIndexCountersAtProc0Only);  // output

      unsigned int auxUnifiedSize = weightSequence.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
      if (m_env.inter0Rank() == 0) {
        UQ_FATAL_TEST_MACRO(unifiedIndexCountersAtProc0Only.size() != auxUnifiedSize,
                            m_env.worldRank(),
                            "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                            "wrong output from sampleIndexesAtProc0() in step 5");
      }

  double stepRunTime = uqMiscGetEllapsedSeconds(&timevalStep);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving uqMLSampling<P_V,P_M>::generateSequence_Step()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", after " << stepRunTime << " seconds"
                            << std::endl;
  }

  return;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::generateSequence_Step06_all(
  const uqMLSamplingLevelOptionsClass* currOptions,                     // input
  unsigned int                         indexOfFirstWeight,              // input
  unsigned int                         indexOfLastWeight,               // input
  const std::vector<unsigned int>&     unifiedIndexCountersAtProc0Only, // input
  bool&                                useBalancedChains,               // output
  std::vector<uqExchangeInfoStruct>&   exchangeStdVec)                  // output
{
  int iRC = UQ_OK_RC;
  struct timeval timevalStep;
  iRC = gettimeofday(&timevalStep, NULL);
  if (iRC) {}; // just to remove compiler warning

  useBalancedChains = decideOnBalancedChains_all(currOptions,                     // input
                                                 indexOfFirstWeight,              // input
                                                 indexOfLastWeight,               // input
                                                 unifiedIndexCountersAtProc0Only, // input
                                                 exchangeStdVec);                 // output

  double stepRunTime = uqMiscGetEllapsedSeconds(&timevalStep);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving uqMLSampling<P_V,P_M>::generateSequence_Step()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", after " << stepRunTime << " seconds"
                            << std::endl;
  }

  return;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::generateSequence_Step07_inter0(
  bool                                      useBalancedChains,               // input
  unsigned int                              indexOfFirstWeight,              // input
  unsigned int                              indexOfLastWeight,               // input
  const std::vector<unsigned int>&          unifiedIndexCountersAtProc0Only, // input
  uqUnbalancedLinkedChainsPerNodeStruct&    unbalancedLinkControl,           // (possible) output
  const uqMLSamplingLevelOptionsClass*      currOptions,                     // input
  const uqSequenceOfVectorsClass<P_V,P_M>&  prevChain,                       // input
  std::vector<uqExchangeInfoStruct>&        exchangeStdVec,                  // (possible) input/output
  uqBalancedLinkedChainsPerNodeStruct<P_V>& balancedLinkControl)             // (possible) output
{
  int iRC = UQ_OK_RC;
  struct timeval timevalStep;
  iRC = gettimeofday(&timevalStep, NULL);
  if (iRC) {}; // just to remove compiler warning

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": beginning step 7 of 11"
                                << std::endl;
      }

      if (useBalancedChains) {
        prepareBalLinkedChains_inter0(currOptions,                     // input
                                      prevChain,                       // input
                                      exchangeStdVec,                  // input/output
                                      balancedLinkControl);            // output
      }
      else {
        prepareUnbLinkedChains_inter0(indexOfFirstWeight,              // input
                                      indexOfLastWeight,               // input
                                      unifiedIndexCountersAtProc0Only, // input
                                      unbalancedLinkControl);          // output
      }

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": balancedLinkControl.balLinkedChains.size() = "   << balancedLinkControl.balLinkedChains.size()
                                << ", unbalancedLinkControl.unbLinkedChains.size() = " << unbalancedLinkControl.unbLinkedChains.size()
                                << std::endl;
      }

  double stepRunTime = uqMiscGetEllapsedSeconds(&timevalStep);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving uqMLSampling<P_V,P_M>::generateSequence_Step()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", after " << stepRunTime << " seconds"
                            << std::endl;
  }

  return;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::generateSequence_Step08_all(
  uqBayesianJointPdfClass<P_V,P_M>& currPdf, // input/output
  uqGenericVectorRVClass<P_V,P_M>&  currRv)  // output
{
  int iRC = UQ_OK_RC;
  struct timeval timevalStep;
  iRC = gettimeofday(&timevalStep, NULL);
  if (iRC) {}; // just to remove compiler warning

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": beginning step 8 of 11"
                                << std::endl;
      }

      currRv.setPdf(currPdf);

  double stepRunTime = uqMiscGetEllapsedSeconds(&timevalStep);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving uqMLSampling<P_V,P_M>::generateSequence_Step()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", after " << stepRunTime << " seconds"
                            << std::endl;
  }

  return;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::generateSequence_Step09_all(
  const uqSequenceOfVectorsClass<P_V,P_M>& prevChain,                         // input
  unsigned int                             indexOfFirstWeight,                // input
  unsigned int                             indexOfLastWeight,                 // input
  const std::vector<double>&               unifiedWeightStdVectorAtProc0Only, // input
  const uqScalarSequenceClass<double>&     weightSequence,                    // input
  double                                   prevEta,                           // input
  const uqGenericVectorRVClass<P_V,P_M>&   currRv,                            // input
  uqMLSamplingLevelOptionsClass*           currOptions,                       // input (changed temporarily internally)
  P_M&                                     unifiedCovMatrix,                  // input/output
  double&                                  currEta)                           // output
{
  int iRC = UQ_OK_RC;
  struct timeval timevalStep;
  iRC = gettimeofday(&timevalStep, NULL);
  if (iRC) {}; // just to remove compiler warning

    if (currOptions->m_scaleCovMatrix == false) {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": skipping step 9 of 11"
                                << std::endl;
      }
    }
    else {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence_Step09_all()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": beginning step 9 of 11"
                                << std::endl;
      }

      double beforeEta           = prevEta;
      double beforeRejectionRate = 0.;               // To be updated
      bool   beforeRejectionRateIsBelowRange = true; // To be updated

      double nowEta           = prevEta;
      double nowRejectionRate = 0.;               // To be computed
      bool   nowRejectionRateIsBelowRange = true; // To be computed

      std::vector<double> etas(2,0.);
      etas[0] = beforeEta;
      etas[1] = 1.;

      std::vector<double> rejs(2,0.);
      rejs[0] = 0.; // To be computed
      rejs[1] = 0.; // To be computed

      unsigned int nowAttempt = 0;
      bool testResult = false;
      double meanRejectionRate = .5*(currOptions->m_minRejectionRate + currOptions->m_maxRejectionRate);
      bool useMiddlePointLogicForEta = false;
      P_M nowCovMatrix(unifiedCovMatrix);
#if 0 // KAUST, to check
      std::vector<double> unifiedWeightStdVectorAtProc0Only(0);
      weightSequence.getUnifiedContentsAtProc0Only(m_vectorSpace.numOfProcsForStorage() == 1,
                                                   unifiedWeightStdVectorAtProc0Only);
#endif
      do {
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence_Step09_all()"
                                  << ", level " << m_currLevel+LEVEL_REF_ID
                                  << ", step "  << m_currStep
                                  << ": entering loop for assessing rejection rate"
                                  << ", with nowAttempt = "  << nowAttempt
                                  << ", nowRejectionRate = " << nowRejectionRate
                                  << std::endl;
        }
        nowCovMatrix = unifiedCovMatrix;

        if (nowRejectionRate < currOptions->m_minRejectionRate) {
          nowRejectionRateIsBelowRange = true;
        }
        else if (nowRejectionRate > currOptions->m_maxRejectionRate) {
          nowRejectionRateIsBelowRange = false;
        }
        else {
          UQ_FATAL_TEST_MACRO(true,
                              m_env.worldRank(),
                              "uqMLSamplingClass<P_V,P_M>::generateSequence_Step09_all()",
                              "nowRejectionRate should be out of the requested range at this point of the logic");
        }

        if (m_env.inter0Rank() >= 0) { // KAUST
          if (nowAttempt > 0) {
            if (useMiddlePointLogicForEta == false) {
              if (nowAttempt == 1) {
                // Ok, keep useMiddlePointLogicForEta = false
              }
              else if ((beforeRejectionRateIsBelowRange == true) &&
                       (nowRejectionRateIsBelowRange    == true)) {
                // Ok
              }
              else if ((beforeRejectionRateIsBelowRange == false) &&
                       (nowRejectionRateIsBelowRange    == false)) {
                // Ok
              }
              else if ((beforeRejectionRateIsBelowRange == true ) &&
                       (nowRejectionRateIsBelowRange    == false)) {
                useMiddlePointLogicForEta = true;

                // This is the first time the middle point logic will be used below
                etas[0] = std::min(beforeEta,nowEta);
                etas[1] = std::max(beforeEta,nowEta);

                if (etas[0] == beforeEta) {
                  rejs[0] = beforeRejectionRate;
                  rejs[1] = nowRejectionRate;
                }
                else {
                  rejs[0] = nowRejectionRate;
                  rejs[1] = beforeRejectionRate;
                }
              }
              else if ((beforeRejectionRateIsBelowRange == false) &&
                       (nowRejectionRateIsBelowRange    == true )) {
                useMiddlePointLogicForEta = true;

                // This is the first time the middle point logic will be used below
                etas[0] = std::min(beforeEta,nowEta);
                etas[1] = std::max(beforeEta,nowEta);
              }
              else {
                UQ_FATAL_TEST_MACRO(true,
                                    m_env.worldRank(),
                                    "uqMLSamplingClass<P_V,P_M>::generateSequence_Step09_all()",
                                    "before and now range flags are inconsistent");
              }
            } // if (useMiddlePointLogicForEta == false)

            beforeEta                       = nowEta;
            beforeRejectionRate             = nowRejectionRate;
            beforeRejectionRateIsBelowRange = nowRejectionRateIsBelowRange;
            if (useMiddlePointLogicForEta == false) {
              if (beforeRejectionRateIsBelowRange) nowEta *= 4.;
              else                                 nowEta /= 4.;
              if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
                *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence_Step09_all()"
                                        << ", level " << m_currLevel+LEVEL_REF_ID
                                        << ", step "  << m_currStep
                                        << ": in loop for assessing rejection rate"
                                        << ", with nowAttempt = "  << nowAttempt
                                        << ", useMiddlePointLogicForEta = false"
                                        << ", nowEta just updated to value (to be tested) " << nowEta
                                        << std::endl;
              }
            }
            else {
              if (nowRejectionRate > meanRejectionRate) {
                if (rejs[0] > meanRejectionRate) {
                  etas[0] = nowEta;
                  etas[1] = etas[1];
                }
                else {
                  etas[0] = etas[0];
                  etas[1] = nowEta;
                }
              }
              else {
                if (rejs[0] < meanRejectionRate) {
                  etas[0] = nowEta;
                  etas[1] = etas[1];
                }
                else {
                  etas[0] = etas[0];
                  etas[1] = nowEta;
                }
              }
              nowEta = .5*(etas[0] + etas[1]);
              if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
                *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence_Step09_all()"
                                        << ", level " << m_currLevel+LEVEL_REF_ID
                                        << ", step "  << m_currStep
                                        << ": in loop for assessing rejection rate"
                                        << ", with nowAttempt = " << nowAttempt
                                        << ", useMiddlePointLogicForEta = true"
                                        << ", nowEta just updated to value (to be tested) " << nowEta
                                        << ", etas[0] = " << etas[0]
                                        << ", etas[1] = " << etas[1]
                                        << std::endl;
              }
            }
          } // if (nowAttempt > 0)
        } // if (m_env.inter0Rank() >= 0) // KAUST

        nowCovMatrix *= nowEta;

        // prudencio 2010-12-09: logic 'originalSubNumSamples += 1' added because of the difference of results between GNU and INTEL compiled codes
        double       doubSubNumSamples     = (1.-meanRejectionRate)/meanRejectionRate/currOptions->m_covRejectionRate/currOptions->m_covRejectionRate; // e.g. 19.99...; or 20.0; or 20.1; or 20.9
        unsigned int originalSubNumSamples = 1 + (unsigned int) (doubSubNumSamples); // e.g. 20; or 21; or 21; or 21
        double       auxDouble             = (double) originalSubNumSamples; // e.g. 20.0; or 21.0; or 21.0; or 21.0
        if ((auxDouble - doubSubNumSamples) < 1.e-8) { // e.g. 0.00...01; or 1.0; or 0.9; or 0.1
          originalSubNumSamples += 1;
        }

        if (m_env.inter0Rank() >= 0) { // KAUST
          if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
            *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence_Step09_all()"
                                    << ", level " << m_currLevel+LEVEL_REF_ID
                                    << ", step "  << m_currStep
                                    << ": in loop for assessing rejection rate"
                                    << ", about to sample "     << originalSubNumSamples << " indexes"
                                    << ", meanRejectionRate = " << meanRejectionRate
                                    << ", covRejectionRate = "  << currOptions->m_covRejectionRate
                                    << std::endl;
          }
        } // KAUST

        std::vector<unsigned int> nowUnifiedIndexCountersAtProc0Only(0); // It will be resized by 'sampleIndexes_proc0()' below
        if (m_env.inter0Rank() >= 0) { // KAUST
          unsigned int tmpUnifiedNumSamples = originalSubNumSamples*m_env.inter0Comm().NumProc();
          sampleIndexes_proc0(tmpUnifiedNumSamples,                // input
                              unifiedWeightStdVectorAtProc0Only,   // input
                              nowUnifiedIndexCountersAtProc0Only); // output

          unsigned int auxUnifiedSize = weightSequence.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
          if (m_env.inter0Rank() == 0) {
            UQ_FATAL_TEST_MACRO(nowUnifiedIndexCountersAtProc0Only.size() != auxUnifiedSize,
                                m_env.worldRank(),
                                "uqMLSamplingClass<P_V,P_M>::generateSequence_Step09_all()",
                                "wrong output from sampleIndexesAtProc0() in step 9");
          }

          if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
            *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence_Step09_all()"
                                    << ", level " << m_currLevel+LEVEL_REF_ID
                                    << ", step "  << m_currStep
                                    << ": in loop for assessing rejection rate"
                                    << ", about to distribute sampled assessment indexes"
                                    << std::endl;
          }
        } // KAUST

        std::vector<uqExchangeInfoStruct>        exchangeStdVec(0);
        uqBalancedLinkedChainsPerNodeStruct<P_V> nowBalLinkControl;
        uqUnbalancedLinkedChainsPerNodeStruct    nowUnbLinkControl; // KAUST

        // All processors should call this routine in order to have the same decision value
        bool useBalancedChains = decideOnBalancedChains_all(currOptions,                        // input
                                                            indexOfFirstWeight,                 // input
                                                            indexOfLastWeight,                  // input
                                                            nowUnifiedIndexCountersAtProc0Only, // input
                                                            exchangeStdVec);                    // output

        if (m_env.inter0Rank() >= 0) { // KAUST
          if (useBalancedChains) {
            prepareBalLinkedChains_inter0(currOptions,                        // input
                                          prevChain,                          // input
                                          exchangeStdVec,                     // input/output
                                          nowBalLinkControl);                 // output
          }
          else {
            prepareUnbLinkedChains_inter0(indexOfFirstWeight,                 // input
                                          indexOfLastWeight,                  // input
                                          nowUnifiedIndexCountersAtProc0Only, // input
                                          nowUnbLinkControl);                 // output
          }
        } // KAUST

        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence_Step09_all()"
                                  << ", level " << m_currLevel+LEVEL_REF_ID
                                  << ", step "  << m_currStep
                                  << ": in loop for assessing rejection rate"
                                  << ", about to generate assessment chain"
                                  << std::endl;
        }

        uqSequenceOfVectorsClass<P_V,P_M> nowChain(m_vectorSpace,
                                                   0,
                                                   m_options.m_prefix+"now_chain");
        double       nowRunTime    = 0.;
        unsigned int nowRejections = 0;

        // KAUST: all nodes should call here
        bool         savedTotallyMute           = currOptions->m_totallyMute; // HERE - ENHANCEMENT
        unsigned int savedRawChainSize          = currOptions->m_rawChainSize; // Ok to use rawChainSize
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
        bool         savedRawChainComputeStats  = currOptions->m_rawChainComputeStats;
#endif
        bool         savedFilteredChainGenerate = currOptions->m_filteredChainGenerate;
        unsigned int savedDrMaxNumExtraStages   = currOptions->m_drMaxNumExtraStages;
        unsigned int savedAmAdaptInterval       = currOptions->m_amAdaptInterval;

        currOptions->m_totallyMute = true;
        if (m_env.displayVerbosity() >= 999999) {
          currOptions->m_totallyMute = false;
        }
        currOptions->m_rawChainSize          = 0; // will be set inside generateXYZLinkedChains()
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
        currOptions->m_rawChainComputeStats  = false;
#endif
        currOptions->m_filteredChainGenerate = false;
        currOptions->m_drMaxNumExtraStages   = 0;
        currOptions->m_amAdaptInterval       = 0;

        // KAUST: all nodes in 'subComm' should call here, important
        if (useBalancedChains) {
          generateBalLinkedChains_all(*currOptions,       // input, only m_rawChainSize changes
                                      nowCovMatrix,       // input
                                      currRv,             // input
                                      nowBalLinkControl,  // input // Round Rock
                                      nowChain,           // output 
                                      nowRunTime,         // output
                                      nowRejections,      // output
                                      NULL,               // output
                                      NULL);              // output
        }
        else {
          generateUnbLinkedChains_all(*currOptions,       // input, only m_rawChainSize changes
                                      nowCovMatrix,       // input
                                      currRv,             // input
                                      nowUnbLinkControl,  // input // Round Rock
                                      indexOfFirstWeight, // input // Round Rock
                                      prevChain,          // input // Round Rock
                                      nowChain,           // output 
                                      nowRunTime,         // output
                                      nowRejections,      // output
                                      NULL,               // output
                                      NULL);              // output
        }

        // KAUST: all nodes should call here
        currOptions->m_totallyMute           = savedTotallyMute;
        currOptions->m_rawChainSize          = savedRawChainSize;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
        currOptions->m_rawChainComputeStats  = savedRawChainComputeStats;
#endif
        currOptions->m_filteredChainGenerate = savedFilteredChainGenerate; // FIX ME
        currOptions->m_drMaxNumExtraStages   = savedDrMaxNumExtraStages;
        currOptions->m_amAdaptInterval       = savedAmAdaptInterval;

        for (unsigned int i = 0; i < nowBalLinkControl.balLinkedChains.size(); ++i) {
          UQ_FATAL_TEST_MACRO(nowBalLinkControl.balLinkedChains[i].initialPosition == NULL,
                              m_env.worldRank(),
                              "uqMLSamplingClass<P_V,P_M>::generateSequence_Step09_all()",
                              "Initial position pointer in step 9 should not be NULL");
          delete nowBalLinkControl.balLinkedChains[i].initialPosition;
          nowBalLinkControl.balLinkedChains[i].initialPosition = NULL;
        }
        nowBalLinkControl.balLinkedChains.clear();

        if (m_env.inter0Rank() >= 0) { // KAUST
          // If only one cov matrix is used, then the rejection should be assessed among all inter0Comm nodes // KAUST3
          unsigned int nowUnifiedRejections = 0;
          m_env.inter0Comm().Allreduce((void *) &nowRejections, (void *) &nowUnifiedRejections, (int) 1, uqRawValue_MPI_UNSIGNED, uqRawValue_MPI_SUM,
                                       "uqMLSamplingClass<P_V,P_M>::generateSequence_Step09_all()",
                                       "failed MPI.Allreduce() for now rejections");

#if 0 // Round Rock 2009 12 29
          unsigned int tmpUnifiedNumSamples = 0;
          m_env.inter0Comm().Allreduce((void *) &tmpSubNumSamples, (void *) &tmpUnifiedNumSamples, (int) 1, uqRawValue_MPI_UNSIGNED, uqRawValue_MPI_SUM,
                                       "uqMLSamplingClass<P_V,P_M>::generateSequence_Step09_all()",
                                       "failed MPI.Allreduce() for num samples in step 9");
#endif

          unsigned int tmpUnifiedNumSamples = originalSubNumSamples*m_env.inter0Comm().NumProc();
          nowRejectionRate = ((double) nowUnifiedRejections) / ((double) tmpUnifiedNumSamples);

          //bool aux1 = (nowRejectionRate == meanRejectionRate);
          bool aux2 = (nowRejectionRate >= currOptions->m_minRejectionRate)
                      &&
                      (nowRejectionRate <= currOptions->m_maxRejectionRate);
          testResult = aux2;

          // Make sure all nodes in 'inter0Comm' have the same value of 'testResult'
          if (uqMiscCheckForSameValueInAllNodes(testResult,
                                                0., // kept 'zero' on 2010/03/03
                                                m_env.inter0Comm(),
                                                "uqMLSamplingClass<P_V,P_M>::generateSequence_Step09_all(), step 9, testResult") == false) {
            if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
              *m_env.subDisplayFile() << "WARNING, In uqMLSampling<P_V,P_M>::generateSequence()"
                                      << ", level "        << m_currLevel+LEVEL_REF_ID
                                      << ", step "         << m_currStep
                                      << ": nowAttempt = " << nowAttempt
                                      << ", uqMiscCheck for 'testResult' detected a problem"
                                      << std::endl;
            }
	  }
        } // if (m_env.inter0Rank() >= 0) { // KAUST

        // KAUST: all nodes in 'subComm' should have the same 'testResult'
        unsigned int tmpUint = (unsigned int) testResult;
        m_env.subComm().Bcast((void *) &tmpUint, (int) 1, uqRawValue_MPI_UNSIGNED, 0, // Yes, 'subComm', important
                              "uqMLSamplingClass<P_V,P_M>::generateSequence_Step09_all()",
                              "failed MPI.Bcast() for testResult");
        testResult = (bool) tmpUint;

        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
          *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence_Step09_all()"
                                  << ", level "              << m_currLevel+LEVEL_REF_ID
                                  << ", step "               << m_currStep
                                  << ": in loop for assessing rejection rate"
                                  << ", nowAttempt = "       << nowAttempt
                                  << ", beforeEta = "        << beforeEta
                                  << ", etas[0] = "          << etas[0]
                                  << ", nowEta = "           << nowEta
                                  << ", etas[1] = "          << etas[1]
                                  << ", minRejectionRate = " << currOptions->m_minRejectionRate
                                  << ", nowRejectionRate = " << nowRejectionRate
                                  << ", maxRejectionRate = " << currOptions->m_maxRejectionRate
                                  << std::endl;
        }
        nowAttempt++;

        if (m_env.inter0Rank() >= 0) { // KAUST
          // Make sure all nodes in 'inter0Comm' have the same value of 'nowEta'
          if (uqMiscCheckForSameValueInAllNodes(nowEta,
                                                1.0e-16, // changed from 'zero' on 2009/11/dd
                                                m_env.inter0Comm(),
                                                "uqMLSamplingClass<P_V,P_M>::generateSequence_Step09_all(), step 9, nowEta") == false) {
            if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
              *m_env.subDisplayFile() << "WARNING, In uqMLSampling<P_V,P_M>::generateSequence()"
                                      << ", level "        << m_currLevel+LEVEL_REF_ID
                                      << ", step "         << m_currStep
                                      << ": nowAttempt = " << nowAttempt
                                      << ", uqMiscCheck for 'nowEta' detected a problem"
                                      << std::endl;
            }
          }
        }
      } while (testResult == false);
      currEta = nowEta;
      if (currEta != 1.) {
        unifiedCovMatrix *= currEta;
      }

      unsigned int quantity1 = weightSequence.unifiedSequenceSize(m_vectorSpace.numOfProcsForStorage() == 1);
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence_Step09_all()"
                                << ", level "                                  << m_currLevel+LEVEL_REF_ID
                                << ", step "                                   << m_currStep
                                << ": weightSequence.subSequenceSize() = "     << weightSequence.subSequenceSize()
                                << ", weightSequence.unifiedSequenceSize() = " << quantity1
                                << ", currEta = "                              << currEta
                                << ", assessed rejection rate = "              << nowRejectionRate
                                << std::endl;
      }
    }

  double stepRunTime = uqMiscGetEllapsedSeconds(&timevalStep);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving uqMLSampling<P_V,P_M>::generateSequence_Step()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", after " << stepRunTime << " seconds"
                            << std::endl;
  }

  return;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::generateSequence_Step10_all(
  uqMLSamplingLevelOptionsClass&                  currOptions,                  // input (changed temporarily internally)
  const P_M&                                      unifiedCovMatrix,             // input
  const uqGenericVectorRVClass  <P_V,P_M>&        currRv,                       // input
  bool                                            useBalancedChains,            // input
  const uqUnbalancedLinkedChainsPerNodeStruct&    unbalancedLinkControl,        // input // Round Rock
  unsigned int                                    indexOfFirstWeight,           // input // Round Rock
  const uqSequenceOfVectorsClass<P_V,P_M>&        prevChain,                    // input // Round Rock
  const uqBalancedLinkedChainsPerNodeStruct<P_V>& balancedLinkControl,          // input // Round Rock
  uqSequenceOfVectorsClass      <P_V,P_M>&        currChain,                    // output
  double&                                         cumulativeRawChainRunTime,    // output
  unsigned int&                                   cumulativeRawChainRejections, // output
  uqScalarSequenceClass         <double>*         currLogLikelihoodValues,      // output
  uqScalarSequenceClass         <double>*         currLogTargetValues)          // output
{
  int iRC = UQ_OK_RC;
  struct timeval timevalStep;
  iRC = gettimeofday(&timevalStep, NULL);
  if (iRC) {}; // just to remove compiler warning

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ": beginning step 10 of 11"
                                << ", currLogLikelihoodValues = " << currLogLikelihoodValues
                                << std::endl;
      }

      // All nodes should call here
      bool         savedTotallyMute           = currOptions.m_totallyMute; // HERE - ENHANCEMENT
      unsigned int savedRawChainSize          = currOptions.m_rawChainSize; // Ok to use rawChainSize
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
      bool         savedRawChainComputeStats  = currOptions.m_rawChainComputeStats;
#endif
      bool         savedFilteredChainGenerate = currOptions.m_filteredChainGenerate;

      currOptions.m_totallyMute = true;
      if (m_env.displayVerbosity() >= 999999) {
        currOptions.m_totallyMute = false;
      }
      currOptions.m_rawChainSize          = 0; // will be set inside generateXYZLinkedChains()
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
      currOptions.m_rawChainComputeStats  = false;
#endif
      currOptions.m_filteredChainGenerate = false;

      // All nodes should call here
      if (useBalancedChains) {
        generateBalLinkedChains_all(currOptions,                  // input, only m_rawChainSize changes
                                    unifiedCovMatrix,             // input
                                    currRv,                       // input
                                    balancedLinkControl,          // input // Round Rock
                                    currChain,                    // output
                                    cumulativeRawChainRunTime,    // output
                                    cumulativeRawChainRejections, // output
                                    currLogLikelihoodValues,      // output // likelihood is important
                                    currLogTargetValues);         // output
      }
      else {
        generateUnbLinkedChains_all(currOptions,                  // input, only m_rawChainSize changes
                                    unifiedCovMatrix,             // input
                                    currRv,                       // input
                                    unbalancedLinkControl,        // input // Round Rock
                                    indexOfFirstWeight,           // input // Round Rock
                                    prevChain,                    // input // Round Rock
                                    currChain,                    // output
                                    cumulativeRawChainRunTime,    // output
                                    cumulativeRawChainRejections, // output
                                    currLogLikelihoodValues,      // output // likelihood is important
                                    currLogTargetValues);         // output
      }

      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
        double tmpValue = INFINITY;
        if (currLogLikelihoodValues) tmpValue = (*currLogLikelihoodValues)[0];
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence_Step()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ", after chain generatrion"
                                << ", currLogLikelihoodValues[0] = " << tmpValue
                                << std::endl;
      }

      // All nodes should call here
      currOptions.m_totallyMute           = savedTotallyMute;
      currOptions.m_rawChainSize          = savedRawChainSize;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
      currOptions.m_rawChainComputeStats  = savedRawChainComputeStats;
#endif
      currOptions.m_filteredChainGenerate = savedFilteredChainGenerate; // FIX ME

  double stepRunTime = uqMiscGetEllapsedSeconds(&timevalStep);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving uqMLSampling<P_V,P_M>::generateSequence_Step()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", after " << stepRunTime << " seconds"
                            << std::endl;
  }

  return;
}

template <class P_V,class P_M>
void
uqMLSamplingClass<P_V,P_M>::generateSequence_Step11_inter0(
  const uqMLSamplingLevelOptionsClass* currOptions,                  // input
  unsigned int                         unifiedRequestedNumSamples,   // input
  unsigned int                         cumulativeRawChainRejections, // input
  uqSequenceOfVectorsClass<P_V,P_M>&   currChain,                    // input/output
  uqScalarSequenceClass<double>&       currLogLikelihoodValues,      // input/output
  uqScalarSequenceClass<double>&       currLogTargetValues,          // input/output
  unsigned int&                        unifiedNumberOfRejections)    // output
{
  int iRC = UQ_OK_RC;
  struct timeval timevalStep;
  iRC = gettimeofday(&timevalStep, NULL);
  if (iRC) {}; // just to remove compiler warning

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ": beginning step 11 of 11"
                            << std::endl;
  }

  //if (m_env.subComm().MyPID() == 0) std::cout << "Aqui 000" << std::endl;

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (currOptions->m_rawChainComputeStats) {
    uqFilePtrSetStruct filePtrSet;
    m_env.openOutputFile(currOptions->m_dataOutputFileName,
                         UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT, // Yes, always ".m"
                         currOptions->m_dataOutputAllowedSet,
                         false,
                         filePtrSet);

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 10)) { // output debug
      *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence_Step()"
                              << ", level " << m_currLevel+LEVEL_REF_ID
                              << ", step "  << m_currStep
                              << ", calling computeStatistics for raw chain"
                              << ". Ofstream pointer value = " << filePtrSet.ofsVar
                              << ", statistical options are"
                              << "\n" << *currOptions->m_rawChainStatisticalOptionsObj
                              << std::endl;
    }
    //m_env.syncPrintDebugMsg("At step 11, calling computeStatistics for raw chain",1,10,m_env.inter0Comm()); // output debug
    currChain.computeStatistics(*currOptions->m_rawChainStatisticalOptionsObj,
                                filePtrSet.ofsVar);

    m_env.closeFile(filePtrSet,UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT);
  }
  // Compute MLE and MAP
  // rr0
#endif
  if (currOptions->m_rawChainDataOutputFileName != UQ_MH_SG_FILENAME_FOR_NO_FILE) {
    currChain.unifiedWriteContents(currOptions->m_rawChainDataOutputFileName,
                                   currOptions->m_rawChainDataOutputFileType); // KAUST5
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence_Step()"
                              << ", level " << m_currLevel+LEVEL_REF_ID
                              << ", step "  << m_currStep
                              << ", before calling currLogLikelihoodValues.unifiedWriteContents()"
                              << ", currLogLikelihoodValues[0] = " << currLogLikelihoodValues[0]
                              << std::endl;
    }
    currLogLikelihoodValues.unifiedWriteContents(currOptions->m_rawChainDataOutputFileName,
                                                 currOptions->m_rawChainDataOutputFileType);
    currLogTargetValues.unifiedWriteContents    (currOptions->m_rawChainDataOutputFileName,
                                                 currOptions->m_rawChainDataOutputFileType);
  }

  if (currOptions->m_filteredChainGenerate) {
    uqFilePtrSetStruct filePtrSet;
    m_env.openOutputFile(currOptions->m_dataOutputFileName,
                         UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT, // Yes, always ".m"
                         currOptions->m_dataOutputAllowedSet,
                         false,
                         filePtrSet);

    unsigned int filterInitialPos = (unsigned int) (currOptions->m_filteredChainDiscardedPortion * (double) currChain.subSequenceSize());
    unsigned int filterSpacing    = currOptions->m_filteredChainLag;
    if (filterSpacing == 0) {
      currChain.computeFilterParams(filePtrSet.ofsVar,
                                    filterInitialPos,
                                    filterSpacing);
    }

    // Filter positions from the converged portion of the chain
    currChain.filter(filterInitialPos,
                     filterSpacing);
    currChain.setName(currOptions->m_prefix + "filtChain");

    currLogLikelihoodValues.filter(filterInitialPos,
                                   filterSpacing);
    currLogLikelihoodValues.setName(currOptions->m_prefix + "filtLogLikelihood");

    currLogTargetValues.filter(filterInitialPos,
                               filterSpacing);
    currLogTargetValues.setName(currOptions->m_prefix + "filtLogTarget");

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    if (currOptions->m_filteredChainComputeStats) {
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 10)) { // output debug
        *m_env.subDisplayFile() << "In uqMLSampling<P_V,P_M>::generateSequence_Step()"
                                << ", level " << m_currLevel+LEVEL_REF_ID
                                << ", step "  << m_currStep
                                << ", calling computeStatistics for filtered chain"
                                << ". Ofstream pointer value = " << filePtrSet.ofsVar
                                << ", statistical options are"
                                << "\n" << *currOptions->m_rawChainStatisticalOptionsObj
                                << std::endl;
      }

      //m_env.syncPrintDebugMsg("At step 11, calling computeStatistics for filtered chain",1,10,m_env.inter0Comm()); // output debug
      currChain.computeStatistics(*currOptions->m_filteredChainStatisticalOptionsObj,
                                  filePtrSet.ofsVar);

    }
#endif
    // Compute MLE and MAP
    // rr0
    m_env.closeFile(filePtrSet,UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT);

    if (currOptions->m_filteredChainDataOutputFileName != UQ_MH_SG_FILENAME_FOR_NO_FILE) {
      currChain.unifiedWriteContents              (currOptions->m_filteredChainDataOutputFileName,
                                                   currOptions->m_filteredChainDataOutputFileType);
      currLogLikelihoodValues.unifiedWriteContents(currOptions->m_filteredChainDataOutputFileName,
                                                   currOptions->m_filteredChainDataOutputFileType);
      currLogTargetValues.unifiedWriteContents    (currOptions->m_filteredChainDataOutputFileName,
                                                   currOptions->m_filteredChainDataOutputFileType);
    }
  } // if (currOptions->m_filteredChainGenerate)

  if (currOptions->m_filteredChainGenerate) {
    // Do not check
  }
  else {
    // Check if unified size of generated chain matches the unified requested size // KAUST
    unsigned int tmpSize = currChain.subSequenceSize();
    unsigned int unifiedGeneratedNumSamples = 0;
    m_env.inter0Comm().Allreduce((void *) &tmpSize, (void *) &unifiedGeneratedNumSamples, (int) 1, uqRawValue_MPI_UNSIGNED, uqRawValue_MPI_SUM,
                                 "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                                 "failed MPI.Allreduce() for generated num samples in step 11");
    //std::cout << "unifiedGeneratedNumSamples = "   << unifiedGeneratedNumSamples
    //          << ", unifiedRequestedNumSamples = " << unifiedRequestedNumSamples
    //          << std::endl;
    UQ_FATAL_TEST_MACRO(unifiedGeneratedNumSamples != unifiedRequestedNumSamples,
                        m_env.worldRank(),
                        "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                        "currChain (linked one) has been generated with invalid size");
  }

  // Compute unified number of rejections
  m_env.inter0Comm().Allreduce((void *) &cumulativeRawChainRejections, (void *) &unifiedNumberOfRejections, (int) 1, uqRawValue_MPI_UNSIGNED, uqRawValue_MPI_SUM,
                               "uqMLSamplingClass<P_V,P_M>::generateSequence()",
                               "failed MPI.Allreduce() for number of rejections");

  double stepRunTime = uqMiscGetEllapsedSeconds(&timevalStep);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "Leaving uqMLSampling<P_V,P_M>::generateSequence_Step()"
                            << ", level " << m_currLevel+LEVEL_REF_ID
                            << ", step "  << m_currStep
                            << ", after " << stepRunTime << " seconds"
                            << std::endl;
  }

  return;
}
#endif // __UQ_MULTI_LEVEL_SAMPLING2_H__


