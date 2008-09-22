/* uq/libs/mcmc/inc/uqMonteCarloDC.h
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

#ifndef __UQ_MCDC_H__
#define __UQ_MCDC_H__

#define UQ_MCDC_FILENAME_FOR_NO_OUTPUT_FILE "."

// _ODV = option default value
#define UQ_MCDC_NUM_SAMPLES_ODV       100
#define UQ_MCDC_OUTPUT_FILE_NAME_ODV  UQ_MCDC_FILENAME_FOR_NO_OUTPUT_FILE
#define UQ_MCDC_USE2_ODV              0
#define UQ_MCDC_DISPLAY_PERIOD_ODV    500
#define UQ_MCDC_MEASURE_RUN_TIMES_ODV 0
#define UQ_MCDC_WRITE_ODV             0
#define UQ_MCDC_COMPUTE_STATS_ODV     0

/*! A templated class that implements a Monte Carlo Distribution Calculator
 */
template <class P_V,class P_M,class Q_V,class Q_M>
class uqMonteCarloDCClass
{
public:
  uqMonteCarloDCClass(const uqEnvironmentClass&                           env,               /*! The QUESO toolkit environment. */
                      const char*                                         prefix,            /*! Prefix.                        */
                      const uqParamSpaceClass          <P_V,P_M>&         paramSpace,        /*! The parameter space.           */
                      const uqQoISpaceClass            <Q_V,Q_M>&         qoiSpace,          /*! The QoI space.                 */
                      const uqProbDensity_BaseClass    <P_V,P_M>*         paramDensityObj,   /*! The parameter density.         */
                      const uqSampleGenerator_BaseClass<P_V,P_M>*         paramGeneratorObj, /*! The parameter generator.       */
                      const uqQoIFunction_BaseClass    <P_V,P_M,Q_V,Q_M>& qoiFunctionObj);   /*! The QoI function.              */
 ~uqMonteCarloDCClass();

  void calculateDistributions();
  void calculateDistributions(const uqSampleGenerator_BaseClass<P_V,P_M>& paramGeneratorObj);

  void print                 (std::ostream&            os) const;

private:
  void defineMyOptions       (po::options_description& optionsDesc);
  void getMyOptionValues     (po::options_description& optionsDesc);

  void calculateDistributions(const uqSampleGenerator_BaseClass<P_V,P_M>& paramGeneratorObj,
                                    uqChainBaseClass<P_V>&                workingSeq);

  void generateSequence      (      unsigned int                          seqSize,
                              const uqSampleGenerator_BaseClass<P_V,P_M>& paramGeneratorObj,
                                    uqChainBaseClass<P_V>&                workingSeq,
                              const std::string&                          seqName);

  const uqEnvironmentClass&                           m_env;
        std::string                                   m_prefix;
  const uqParamSpaceClass          <P_V,P_M>&         m_paramSpace;
  const uqQoISpaceClass            <Q_V,Q_M>&         m_qoiSpace;
  const uqProbDensity_BaseClass    <P_V,P_M>*         m_paramDensityObj;
  const uqSampleGenerator_BaseClass<P_V,P_M>*         m_paramGeneratorObj;
  const uqQoIFunction_BaseClass    <P_V,P_M,Q_V,Q_M>& m_qoiFunctionObj;

  po::options_description*        m_optionsDesc;
  std::string                     m_option_help;
  std::string                     m_option_numSamples;
  std::string                     m_option_outputFileName;
  std::string                     m_option_use2;
  std::string                     m_option_displayPeriod;
  std::string                     m_option_measureRunTimes;
  std::string                     m_option_write;
  std::string                     m_option_computeStats;

  unsigned int                    m_numSamples;
  std::string                     m_outputFileName;
  bool                            m_use2;
  unsigned int                    m_displayPeriod;
  bool                            m_measureRunTimes;
  bool                            m_write;
  bool                            m_computeStats;
  uqChainStatisticalOptionsClass* m_statisticalOptions;

  uqSequenceOfVectorsClass<P_V>   m_seq1;
  uqArrayOfSequencesClass<P_V>    m_seq2;
};

template<class P_V,class P_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const uqMonteCarloDCClass<P_V,P_M,Q_V,Q_M>& obj);

template <class P_V,class P_M,class Q_V,class Q_M>
uqMonteCarloDCClass<P_V,P_M,Q_V,Q_M>::uqMonteCarloDCClass(
  const uqEnvironmentClass&                           env,
  const char*                                         prefix,
  const uqParamSpaceClass          <P_V,P_M>&         paramSpace,
  const uqQoISpaceClass            <Q_V,Q_M>&         qoiSpace,
  const uqProbDensity_BaseClass    <P_V,P_M>*         paramDensityObj,
  const uqSampleGenerator_BaseClass<P_V,P_M>*         paramGeneratorObj,
  const uqQoIFunction_BaseClass    <P_V,P_M,Q_V,Q_M>& qoiFunctionObj)
  :
  m_env                   (env),
  m_prefix                ((std::string)prefix + "mcdc_"),
  m_paramSpace            (paramSpace),
  m_qoiSpace              (qoiSpace),
  m_paramDensityObj       (paramDensityObj),
  m_paramGeneratorObj     (paramGeneratorObj),
  m_qoiFunctionObj        (qoiFunctionObj),
  m_optionsDesc           (new po::options_description("Monte Carlo options")),
  m_option_help           (m_prefix + "help"           ),
  m_option_numSamples     (m_prefix + "numSamples"     ),
  m_option_outputFileName (m_prefix + "outputFileName" ),
  m_option_use2           (m_prefix + "use2"           ),
  m_option_displayPeriod  (m_prefix + "displayPeriod"  ),
  m_option_measureRunTimes(m_prefix + "measureRunTimes"),
  m_option_write          (m_prefix + "write"          ),
  m_option_computeStats   (m_prefix + "computeStats"   ),
  m_numSamples            (UQ_MCDC_NUM_SAMPLES_ODV      ),
  m_outputFileName        (UQ_MCDC_OUTPUT_FILE_NAME_ODV ),
  m_use2                  (UQ_MCDC_USE2_ODV             ),
  m_displayPeriod         (UQ_MCDC_DISPLAY_PERIOD_ODV   ),
  m_measureRunTimes       (UQ_MCDC_MEASURE_RUN_TIMES_ODV),
  m_write                 (UQ_MCDC_WRITE_ODV            ),
  m_computeStats          (UQ_MCDC_COMPUTE_STATS_ODV    ),
  m_statisticalOptions    (NULL),
  m_seq1                  (0,m_qoiSpace.zeroVector()),
  m_seq2                  (0,m_qoiSpace.zeroVector())
{
  if (m_env.rank() == 0) std::cout << "Entering uqMonteCarloDCClass<P_V,P_M,Q_V,Q_M>::constructor()"
                                   << std::endl;

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.rank() == 0) std::cout << "In uqMonteCarloDCClass<P_V,P_M,Q_V,Q_M>::constructor()"
                                   << ": after getting values of options with prefix '" << m_prefix
                                   << "', state of  object is:"
                                   << "\n" << *this
                                   << std::endl;

  if (m_computeStats) m_statisticalOptions = new uqChainStatisticalOptionsClass(m_env,m_prefix + "seq_");

  if (m_env.rank() == 0) std::cout << "Leaving uqMonteCarloDCClass<P_V,P_M,Q_V,Q_M>::constructor()"
                                   << std::endl;
}

template <class P_V,class P_M,class Q_V,class Q_M>
uqMonteCarloDCClass<P_V,P_M,Q_V,Q_M>::~uqMonteCarloDCClass()
{
  if (m_statisticalOptions) delete m_statisticalOptions;
}

template<class P_V,class P_M,class Q_V,class Q_M>
void
uqMonteCarloDCClass<P_V,P_M,Q_V,Q_M>::defineMyOptions(
  po::options_description& optionsDesc)
{
  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                                      "produce help message for Monte Carlo distribution calculator")
    (m_option_numSamples.c_str(),      po::value<unsigned int>()->default_value(UQ_MCDC_NUM_SAMPLES_ODV      ), "number of samples"                                           )
    (m_option_outputFileName.c_str(),  po::value<std::string >()->default_value(UQ_MCDC_OUTPUT_FILE_NAME_ODV ), "name of output file"                                         )
    (m_option_use2.c_str(),            po::value<bool        >()->default_value(UQ_MCDC_USE2_ODV             ), "use seq2"                                                    )
    (m_option_displayPeriod.c_str(),   po::value<unsigned int>()->default_value(UQ_MCDC_DISPLAY_PERIOD_ODV   ), "period of message display during sequence generation"        )
    (m_option_measureRunTimes.c_str(), po::value<bool        >()->default_value(UQ_MCDC_MEASURE_RUN_TIMES_ODV), "measure run times"                                           )
    (m_option_write.c_str(),           po::value<bool        >()->default_value(UQ_MCDC_WRITE_ODV            ), "write sequence values to the output file"                    )
    (m_option_computeStats.c_str(),    po::value<bool        >()->default_value(UQ_MCDC_COMPUTE_STATS_ODV    ), "compute statistics on sequence of qoi"                       )
  ;

  return;
}

template<class P_V,class P_M,class Q_V,class Q_M>
void
uqMonteCarloDCClass<P_V,P_M,Q_V,Q_M>::getMyOptionValues(
  po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    std::cout << optionsDesc
              << std::endl;
  }

  if (m_env.allOptionsMap().count(m_option_numSamples.c_str())) {
    m_numSamples = m_env.allOptionsMap()[m_option_numSamples.c_str()].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_outputFileName.c_str())) {
    m_outputFileName = m_env.allOptionsMap()[m_option_outputFileName.c_str()].as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_use2.c_str())) {
    m_use2 = m_env.allOptionsMap()[m_option_use2.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_displayPeriod.c_str())) {
    m_displayPeriod = m_env.allOptionsMap()[m_option_displayPeriod.c_str()].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_measureRunTimes.c_str())) {
    m_measureRunTimes = m_env.allOptionsMap()[m_option_measureRunTimes.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_write.c_str())) {
    m_write = m_env.allOptionsMap()[m_option_write.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_computeStats.c_str())) {
    m_computeStats = m_env.allOptionsMap()[m_option_computeStats.c_str()].as<bool>();
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqMonteCarloDCClass<P_V,P_M,Q_V,Q_M>::calculateDistributions()
{
  UQ_FATAL_TEST_MACRO(m_paramGeneratorObj == NULL,
                      m_env.rank(),
                      "uqMonteCarloDCClass<P_V,P_M,Q_V,Q_M>::calculateDistributions()",
                      "m_paramGeneratorObj is NULL");

  if (m_use2) {
    calculateDistributions(*m_paramGeneratorObj,
                           m_seq2);
  }
  else {
    calculateDistributions(*m_paramGeneratorObj,
                           m_seq1);
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqMonteCarloDCClass<P_V,P_M,Q_V,Q_M>::calculateDistributions(
  const uqSampleGenerator_BaseClass<P_V,P_M>& paramGeneratorObj)
{
  if (m_use2) {
    calculateDistributions(paramGeneratorObj,
                           m_seq2);
  }
  else {
    calculateDistributions(paramGeneratorObj,
                           m_seq1);
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqMonteCarloDCClass<P_V,P_M,Q_V,Q_M>::calculateDistributions(
  const uqSampleGenerator_BaseClass<P_V,P_M>& paramGeneratorObj,
        uqChainBaseClass<P_V>&                workingSeq)
{
  std::string prefixName = m_prefix;
  std::string seqName    = prefixName + "seq";

  //****************************************************
  // Generate sequence of qoi values
  //****************************************************
  unsigned int actualNumSamples = std::min(m_numSamples,paramGeneratorObj.period());
  generateSequence(actualNumSamples,
                   paramGeneratorObj,
                   workingSeq,
                   seqName);

  //****************************************************
  // Open file      
  //****************************************************
  std::ofstream* ofs = NULL;
  if (m_outputFileName == UQ_MCDC_FILENAME_FOR_NO_OUTPUT_FILE) {
    if (m_env.rank() == 0) {
      std::cout << "No output file opened for qoi sequence " << seqName 
                << std::endl;
    }
  }
  else {
    if (m_env.rank() == 0) {
      std::cout << "Opening output file '" << m_outputFileName
                << "' for qoi sequence "   << seqName
                << std::endl;
    }

    // Open file
    ofs = new std::ofstream(m_outputFileName.c_str(), std::ofstream::out | std::ofstream::in | std::ofstream::ate);
    if ((ofs            == NULL ) ||
        (ofs->is_open() == false)) {
      delete ofs;
      ofs = new std::ofstream(m_outputFileName.c_str(), std::ofstream::out | std::ofstream::trunc);
    }

    UQ_FATAL_TEST_MACRO((ofs && ofs->is_open()) == false,
                        m_env.rank(),
                        "uqMonteCarloDCClass<P_V,P_M,Q_V,Q_M>::calculateDistributions()",
                        "failed to open file");
  }
  
  //****************************************************
  // Eventually:
  // --> write sequence
  // --> compute statistics on it
  //****************************************************
  if (m_write && ofs) {
    workingSeq.write(seqName,*ofs);
  }

  if (m_computeStats) {
    workingSeq.computeStatistics(*m_statisticalOptions,
                                 seqName,
                                 m_qoiSpace.componentsNames(),
                                 ofs);
  }

  //****************************************************
  // Close file      
  //****************************************************
  if (ofs) {
    // Close file
    ofs->close();

    if (m_env.rank() == 0) {
      std::cout << "Closed output file '" << m_outputFileName
                << "' for qoi sequence "  << seqName
                << std::endl;
    }
  }
  if (m_env.rank() == 0) {
    std::cout << std::endl;
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqMonteCarloDCClass<P_V,P_M,Q_V,Q_M>::generateSequence(
        unsigned int                          seqSize,
  const uqSampleGenerator_BaseClass<P_V,P_M>& paramGeneratorObj,
        uqChainBaseClass<P_V>&                workingSeq,
  const std::string&                          seqName)
{
  if (m_env.rank() == 0) {
    std::cout << "Starting the generation of qoi sequence " << seqName
              << ", with "                                  << seqSize
              << " samples..."
              << std::endl;
  }

  int iRC = UQ_OK_RC;
  struct timeval timevalSeq;
  struct timeval timevalQoIFunction;

  double seqRunTime         = 0;
  double qoiFunctionRunTime = 0;

  iRC = gettimeofday(&timevalSeq, NULL);

  workingSeq.resizeSequence(seqSize);
  P_V tmpV(m_paramSpace.zeroVector());
  Q_V tmpQ(m_qoiSpace.zeroVector());
  for (unsigned int i = 0; i < seqSize; ++i) {
    paramGeneratorObj.nextSample(tmpV);

    if (m_measureRunTimes) iRC = gettimeofday(&timevalQoIFunction, NULL);
    m_qoiFunctionObj.computeQoIs(tmpV,tmpQ);
    if (m_measureRunTimes) qoiFunctionRunTime += uqMiscGetEllapsedSeconds(&timevalQoIFunction);

    workingSeq.setPositionValues(i,tmpQ);

    if ((m_displayPeriod            > 0) && 
        (((i+1) % m_displayPeriod) == 0)) {
      if (m_env.rank() == 0) {
        std::cout << "Finished generating " << i+1
                  << " qoi samples"
                  << std::endl;
      }
    }
  }

  seqRunTime = uqMiscGetEllapsedSeconds(&timevalSeq);

  if (m_env.rank() == 0) {
    std::cout << "Finished the generation of qoi sequence " << seqName
              << ", with "                                  << workingSeq.sequenceSize()
              << " samples";
  }
  std::cout << "\nSome information about this sequence:"
            << "\n  Sequence run time = " << seqRunTime
            << " seconds";
  if (m_measureRunTimes) {
    std::cout << "\n\n Breaking of the seq run time:\n";
    std::cout << "\n  QoI function run time   = " << qoiFunctionRunTime
              << " seconds ("                     << 100.*qoiFunctionRunTime/seqRunTime
              << "%)";
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqMonteCarloDCClass<P_V,P_M,Q_V,Q_M>::print(std::ostream& os) const
{
  os <<         m_option_numSamples      << " = " << m_numSamples
     << "\n" << m_option_outputFileName  << " = " << m_outputFileName
     << "\n" << m_option_use2            << " = " << m_use2
     << "\n" << m_option_displayPeriod   << " = " << m_displayPeriod
     << "\n" << m_option_measureRunTimes << " = " << m_measureRunTimes
     << "\n" << m_option_write           << " = " << m_write
     << "\n" << m_option_computeStats    << " = " << m_computeStats;

  return;
}

template<class P_V,class P_M,class Q_V,class Q_M> 
std::ostream& operator<<(std::ostream& os, const uqMonteCarloDCClass<P_V,P_M,Q_V,Q_M>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_MCDC_H__
