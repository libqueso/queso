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

#ifndef __UQ_MOC_SG_H__
#define __UQ_MOC_SG_H__

#include <uqVectorRV.h>
#include <uqVectorFunction.h>
#include <uqVectorFunctionSynchronizer.h>

#define UQ_MOC_SG_FILENAME_FOR_NO_OUTPUT_FILE "."

// _ODV = option default value
#define UQ_MOC_SG_OUTPUT_FILE_NAME_ODV       UQ_MOC_SG_FILENAME_FOR_NO_OUTPUT_FILE
#define UQ_MOC_SG_OUTPUT_ALLOW_ODV           ""

#define UQ_MOC_SG_PSEQ_OUTPUT_FILE_NAME_ODV  UQ_MOC_SG_FILENAME_FOR_NO_OUTPUT_FILE
#define UQ_MOC_SG_PSEQ_OUTPUT_ALLOW_ODV      ""
#define UQ_MOC_SG_PSEQ_COMPUTE_STATS_ODV     0

#define UQ_MOC_SG_QSEQ_SIZE_ODV              100
#define UQ_MOC_SG_QSEQ_DISPLAY_PERIOD_ODV    500
#define UQ_MOC_SG_QSEQ_MEASURE_RUN_TIMES_ODV 0
#define UQ_MOC_SG_QSEQ_OUTPUT_FILE_NAME_ODV  UQ_MOC_SG_FILENAME_FOR_NO_OUTPUT_FILE
#define UQ_MOC_SG_QSEQ_OUTPUT_ALLOW_ODV      ""
#define UQ_MOC_SG_QSEQ_COMPUTE_STATS_ODV     0

/*! A templated class that implements a Monte Carlo Distribution Calculator
 */
template <class P_V,class P_M,class Q_V,class Q_M>
class uqMonteCarloSGClass
{
public:
  uqMonteCarloSGClass(const char*                                       prefix,         /*! Prefix.           */
                      const uqBaseVectorRVClass      <P_V,P_M>&         paramRv,        /*! The parameter rv. */
                      const uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>& qoiFunctionObj, /*! The qoi function. */
                            uqBaseVectorRVClass      <Q_V,Q_M>&         qoiRv);         /*! The qoi rv.       */
 ~uqMonteCarloSGClass();

  void generateSequence           (uqBaseVectorSequenceClass<P_V,P_M>& workingPSeq,
                                   uqBaseVectorSequenceClass<Q_V,Q_M>& workingQSeq);
  void checkTheParallelEnvironment();

  void print                      (std::ostream& os) const;

private:
  void defineMyOptions       (po::options_description& optionsDesc);
  void getMyOptionValues     (po::options_description& optionsDesc);

  void internGenerateSequence(const uqBaseVectorRVClass      <P_V,P_M>& paramRv,
                                    uqBaseVectorSequenceClass<P_V,P_M>& workingPSeq,
                                    uqBaseVectorSequenceClass<Q_V,Q_M>& workingQSeq);

  void actualGenerateSequence(const uqBaseVectorRVClass      <P_V,P_M>& paramRv,
                                    uqBaseVectorSequenceClass<P_V,P_M>& workingPSeq,
                                    uqBaseVectorSequenceClass<Q_V,Q_M>& workingQSeq,
                                    unsigned int                        seqSize);

  const uqBaseEnvironmentClass&                             m_env;
        std::string                                         m_prefix;
  const uqBaseVectorRVClass              <P_V,P_M>&         m_paramRv;
  const uqBaseVectorFunctionClass        <P_V,P_M,Q_V,Q_M>& m_qoiFunctionObj;
  const uqBaseVectorRVClass              <Q_V,Q_M>&         m_qoiRv;
  const uqVectorSpaceClass               <P_V,P_M>&         m_paramSpace;
  const uqVectorSpaceClass               <Q_V,Q_M>&         m_qoiSpace;
  const uqVectorFunctionSynchronizerClass<P_V,P_M,Q_V,Q_M>* m_qoiFunctionSynchronizer;

  po::options_description*        m_optionsDesc;
  std::string                     m_option_help;
  std::string                     m_option_outputFileName;
  std::string                     m_option_outputAllow;

  std::string                     m_option_pseq_outputFileName;
  std::string                     m_option_pseq_outputAllow;
  std::string                     m_option_pseq_computeStats;

  std::string                     m_option_qseq_size;
  std::string                     m_option_qseq_displayPeriod;
  std::string                     m_option_qseq_measureRunTimes;
  std::string                     m_option_qseq_outputFileName;
  std::string                     m_option_qseq_outputAllow;
  std::string                     m_option_qseq_computeStats;

  std::string                     m_outputFileName;
  std::set<unsigned int>          m_outputAllow;

  std::string                     m_pseqOutputFileName;
  std::set<unsigned int>          m_pseqOutputAllow;
  bool                            m_pseqComputeStats;
  uqChainStatisticalOptionsClass* m_pseqStatisticalOptions;

  unsigned int                    m_qseqSize;
  unsigned int                    m_qseqDisplayPeriod;
  bool                            m_qseqMeasureRunTimes;
  std::string                     m_qseqOutputFileName;
  std::set<unsigned int>          m_qseqOutputAllow;
  bool                            m_qseqComputeStats;
  uqChainStatisticalOptionsClass* m_qseqStatisticalOptions;
};

template<class P_V,class P_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>& obj);

template <class P_V,class P_M,class Q_V,class Q_M>
uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::uqMonteCarloSGClass(
  const char*                                       prefix,
  const uqBaseVectorRVClass      <P_V,P_M>&         paramRv,
  const uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>& qoiFunctionObj,
        uqBaseVectorRVClass      <Q_V,Q_M>&         qoiRv)
  :
  m_env                        (paramRv.env()),
  m_prefix                     ((std::string)(prefix) + "mc_"),
  m_paramRv                    (paramRv),
  m_qoiFunctionObj             (qoiFunctionObj),
  m_qoiRv                      (qoiRv),
  m_paramSpace                 (m_paramRv.imageSet().vectorSpace()),
  m_qoiSpace                   (m_qoiRv.imageSet().vectorSpace()),
  m_qoiFunctionSynchronizer    (new uqVectorFunctionSynchronizerClass<P_V,P_M,Q_V,Q_M>(m_qoiFunctionObj,m_paramRv.imageSet().vectorSpace().zeroVector(),m_qoiRv.imageSet().vectorSpace().zeroVector())),
  m_optionsDesc                (new po::options_description("Monte Carlo options")),
  m_option_help                (m_prefix + "help"                ),
  m_option_outputFileName      (m_prefix + "outputFileName"      ),
  m_option_outputAllow         (m_prefix + "outputAllow"         ),
  m_option_pseq_outputFileName (m_prefix + "pseq_outputFileName" ),
  m_option_pseq_outputAllow    (m_prefix + "pseq_outputAllow"    ),
  m_option_pseq_computeStats   (m_prefix + "pseq_computeStats"   ),
  m_option_qseq_size           (m_prefix + "qseq_size"           ),
  m_option_qseq_displayPeriod  (m_prefix + "qseq_displayPeriod"  ),
  m_option_qseq_measureRunTimes(m_prefix + "qseq_measureRunTimes"),
  m_option_qseq_outputFileName (m_prefix + "qseq_outputFileName" ),
  m_option_qseq_outputAllow    (m_prefix + "qseq_outputAllow"    ),
  m_option_qseq_computeStats   (m_prefix + "qseq_computeStats"   ),
  m_outputFileName             (UQ_MOC_SG_OUTPUT_FILE_NAME_ODV      ),
//m_outputAllow                (),
  m_pseqOutputFileName         (UQ_MOC_SG_PSEQ_OUTPUT_FILE_NAME_ODV ),
//m_pseqOutputAllow            (),
  m_pseqComputeStats           (UQ_MOC_SG_PSEQ_COMPUTE_STATS_ODV    ),
  m_pseqStatisticalOptions     (NULL),
  m_qseqSize                   (UQ_MOC_SG_QSEQ_SIZE_ODV             ),
  m_qseqDisplayPeriod          (UQ_MOC_SG_QSEQ_DISPLAY_PERIOD_ODV   ),
  m_qseqMeasureRunTimes        (UQ_MOC_SG_QSEQ_MEASURE_RUN_TIMES_ODV),
  m_qseqOutputFileName         (UQ_MOC_SG_QSEQ_OUTPUT_FILE_NAME_ODV ),
//m_qseqOutputAllow            (),
  m_qseqComputeStats           (UQ_MOC_SG_QSEQ_COMPUTE_STATS_ODV    ),
  m_qseqStatisticalOptions     (NULL)
{
  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "Entering uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::constructor()"
                           << std::endl;
  }

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::constructor()"
                           << ": after getting values of options with prefix '" << m_prefix
                           << "', state of  object is:"
                           << "\n" << *this
                           << std::endl;
  }

  if (m_pseqComputeStats) m_pseqStatisticalOptions = new uqChainStatisticalOptionsClass(m_env,m_prefix + "pseq_");
  if (m_qseqComputeStats) m_qseqStatisticalOptions = new uqChainStatisticalOptionsClass(m_env,m_prefix + "qseq_");

  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "Leaving uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::constructor()"
                           << std::endl;
  }
}

template <class P_V,class P_M,class Q_V,class Q_M>
uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::~uqMonteCarloSGClass()
{
  if (m_qseqStatisticalOptions ) delete m_qseqStatisticalOptions;
  if (m_pseqStatisticalOptions ) delete m_pseqStatisticalOptions;
  if (m_qoiFunctionSynchronizer) delete m_qoiFunctionSynchronizer;
}

template<class P_V,class P_M,class Q_V,class Q_M>
void
uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::defineMyOptions(
  po::options_description& optionsDesc)
{
  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                                                 "produce help message for Monte Carlo distribution calculator")
    (m_option_outputFileName.c_str(),       po::value<std::string >()->default_value(UQ_MOC_SG_OUTPUT_FILE_NAME_ODV      ), "name of generic output file"                                 )
    (m_option_outputAllow.c_str(),          po::value<std::string >()->default_value(UQ_MOC_SG_OUTPUT_ALLOW_ODV          ), "subEnvs that will write to generic output file"              )
    (m_option_pseq_outputFileName.c_str(),  po::value<std::string >()->default_value(UQ_MOC_SG_PSEQ_OUTPUT_FILE_NAME_ODV ), "name of output file for parameters"                          )
    (m_option_pseq_outputAllow.c_str(),     po::value<std::string >()->default_value(UQ_MOC_SG_PSEQ_OUTPUT_ALLOW_ODV     ), "subEnvs that will write to output file for parameters"       )
    (m_option_pseq_computeStats.c_str(),    po::value<bool        >()->default_value(UQ_MOC_SG_PSEQ_COMPUTE_STATS_ODV    ), "compute statistics on sequence of parameter"                 )
    (m_option_qseq_size.c_str(),            po::value<unsigned int>()->default_value(UQ_MOC_SG_QSEQ_SIZE_ODV             ), "size of qoi sequence"                                        )
    (m_option_qseq_displayPeriod.c_str(),   po::value<unsigned int>()->default_value(UQ_MOC_SG_QSEQ_DISPLAY_PERIOD_ODV   ), "period of message display during qoi sequence generation"    )
    (m_option_qseq_measureRunTimes.c_str(), po::value<bool        >()->default_value(UQ_MOC_SG_QSEQ_MEASURE_RUN_TIMES_ODV), "measure run times"                                           )
    (m_option_qseq_outputFileName.c_str(),  po::value<std::string >()->default_value(UQ_MOC_SG_QSEQ_OUTPUT_FILE_NAME_ODV ), "name of output file for qois"                                )
    (m_option_qseq_outputAllow.c_str(),     po::value<std::string >()->default_value(UQ_MOC_SG_QSEQ_OUTPUT_ALLOW_ODV     ), "subEnvs that will write to output file for qois"             )
    (m_option_qseq_computeStats.c_str(),    po::value<bool        >()->default_value(UQ_MOC_SG_QSEQ_COMPUTE_STATS_ODV    ), "compute statistics on sequence of qoi"                       )
  ;

  return;
}

template<class P_V,class P_M,class Q_V,class Q_M>
void
uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::getMyOptionValues(
  po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    if (m_env.subScreenFile()) {
      *m_env.subScreenFile() << optionsDesc
                             << std::endl;
    }
  }

  if (m_env.allOptionsMap().count(m_option_outputFileName.c_str())) {
    m_outputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_outputFileName.c_str()]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_outputAllow.c_str())) {
    m_outputAllow.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_outputAllow.c_str()].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpAllow);

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_outputAllow.insert((unsigned int) tmpAllow[i]);
      }
    }
  }

  if (m_env.allOptionsMap().count(m_option_pseq_outputFileName.c_str())) {
    m_pseqOutputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_pseq_outputFileName.c_str()]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_pseq_outputAllow.c_str())) {
    m_pseqOutputAllow.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_pseq_outputAllow.c_str()].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpAllow);

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_pseqOutputAllow.insert((unsigned int) tmpAllow[i]);
      }
    }
  }

  if (m_env.allOptionsMap().count(m_option_pseq_computeStats.c_str())) {
    m_pseqComputeStats = ((const po::variable_value&) m_env.allOptionsMap()[m_option_pseq_computeStats.c_str()]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_qseq_size.c_str())) {
    m_qseqSize = ((const po::variable_value&) m_env.allOptionsMap()[m_option_qseq_size.c_str()]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_qseq_displayPeriod.c_str())) {
    m_qseqDisplayPeriod = ((const po::variable_value&) m_env.allOptionsMap()[m_option_qseq_displayPeriod.c_str()]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_qseq_measureRunTimes.c_str())) {
    m_qseqMeasureRunTimes = ((const po::variable_value&) m_env.allOptionsMap()[m_option_qseq_measureRunTimes.c_str()]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_qseq_outputFileName.c_str())) {
    m_qseqOutputFileName = ((const po::variable_value&) m_env.allOptionsMap()[m_option_qseq_outputFileName.c_str()]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_qseq_outputAllow.c_str())) {
    m_qseqOutputAllow.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_qseq_outputAllow.c_str()].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpAllow);

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_qseqOutputAllow.insert((unsigned int) tmpAllow[i]);
      }
    }
  }

  if (m_env.allOptionsMap().count(m_option_qseq_computeStats.c_str())) {
    m_qseqComputeStats = ((const po::variable_value&) m_env.allOptionsMap()[m_option_qseq_computeStats.c_str()]).as<bool>();
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::generateSequence(
  uqBaseVectorSequenceClass<P_V,P_M>& workingPSeq,
  uqBaseVectorSequenceClass<Q_V,Q_M>& workingQSeq)
{
  checkTheParallelEnvironment();
  internGenerateSequence(m_paramRv,workingPSeq,workingQSeq);

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence(
  const uqBaseVectorRVClass      <P_V,P_M>& paramRv,
        uqBaseVectorSequenceClass<P_V,P_M>& workingPSeq,
        uqBaseVectorSequenceClass<Q_V,Q_M>& workingQSeq)
{
  workingPSeq.setName(m_prefix+"ParamSeq");
  workingQSeq.setName(m_prefix+"QoiSeq");

  //****************************************************
  // Generate sequence of qoi values
  //****************************************************
  unsigned int actualSize = std::min(m_qseqSize,paramRv.realizer().period());
  if ((m_env.subScreenFile()) && (m_env.verbosity() >= 0)) {
    *m_env.subScreenFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                           << ": m_qseqSize = "                                             << m_qseqSize
                           << ", paramRv.realizer().period() = "                            << paramRv.realizer().period()
                           << ", about to call actualGenerateSequence() with actualSize = " << actualSize
                           << std::endl;
  }
  actualGenerateSequence(paramRv,
                         workingPSeq,
                         workingQSeq,
                         actualSize);
  if ((m_env.subScreenFile()) && (m_env.verbosity() >= 0)) {
    *m_env.subScreenFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                           << ": returned from call to actualGenerateSequence() with actualSize = " << actualSize
                           << std::endl;
  }

  //****************************************************
  // Open generic output file      
  //****************************************************
  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                           << ", prefix = "                                                        << m_prefix
                           << ": checking necessity of opening generic output file (qseq name is " << workingQSeq.name()
                           << ") ..."
                           << std::endl;
  }
  std::ofstream* genericOfsVar = NULL;
  m_env.openOutputFile(m_outputFileName,
                       "m",
                       m_outputAllow,
                       false,
                       genericOfsVar);

  //****************************************************
  // Eventually:
  // --> write parameter sequence
  // --> compute statistics on it
  //****************************************************
  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                           << ", prefix = "                                            << m_prefix
                           << ": checking necessity of opening output files for pseq " << workingPSeq.name()
                           << "..."
                           << std::endl;
  }

  // Take "sub" care of pseq
  std::ofstream* pseqOfsVar = NULL;
  m_env.openOutputFile(m_pseqOutputFileName,
                       "hdf",
                       m_pseqOutputAllow,
                       true,
                       pseqOfsVar);

  if (pseqOfsVar) {
    workingPSeq.printContents(*pseqOfsVar);
  }

  if (pseqOfsVar) {
    pseqOfsVar->close();
    if (m_env.subScreenFile()) {
      *m_env.subScreenFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                             << ", prefix = "            << m_prefix
                             << ": closed output file '" << m_pseqOutputFileName
                             << "' for pseq "            << workingPSeq.name()
                             << std::endl;
    }
  }

  // Take "unified" care of pseq
#if 0
  std::ofstream* unifiedPSeqOfsVar = NULL;
  m_env.openUnifiedOutputFile(m_pseqOutputFileName,
                              "hdf",
                              true,
                              unifiedPSeqOfsVar);

  if (unifiedPSeqOfsVar) {
    workingPSeq.printUnifiedContents(*unifiedPSeqOfsVar);
  }

  if (unifiedPSeqOfsVar) {
    unifiedPSeqOfsVar->close();
    if (m_env.subScreenFile()) {
      *m_env.subScreenFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                             << ", prefix = "                    << m_prefix
                             << ": closed unified output file '" << m_pseqOutputFileName
                             << "' for pseq "                    << workingPSeq.name()
                             << std::endl;
    }
  }
#else
  if (m_pseqOutputFileName != ".") {
    workingPSeq.printUnifiedContents(m_pseqOutputFileName);
    if (m_env.subScreenFile()) {
      *m_env.subScreenFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                             << ", prefix = "                    << m_prefix
                             << ": closed unified output file '" << m_pseqOutputFileName
                             << "' for pseq "                    << workingPSeq.name()
                             << std::endl;
    }
  }
#endif

  // Take case of other aspects of pseq
  if (m_pseqComputeStats) {
    if ((m_env.subScreenFile()) && (m_env.verbosity() >= 0)) {
      *m_env.subScreenFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                             << ": about to call 'workingPSeq.computeStatistics()'"
                             << std::endl;
    }
    workingPSeq.computeStatistics(*m_pseqStatisticalOptions,
                                  genericOfsVar);
    if ((m_env.subScreenFile()) && (m_env.verbosity() >= 0)) {
      *m_env.subScreenFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                             << ": returned from call to 'workingPSeq.computeStatistics()'"
                             << std::endl;
    }
  }

  //****************************************************
  // Eventually:
  // --> write qoi sequence
  // --> compute statistics on it
  //****************************************************
  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                           << ", prefix = "                                            << m_prefix
                           << ": checking necessity of opening output files for qseq " << workingQSeq.name()
                           << "..."
                           << std::endl;
  }

  // Take "sub" care of qseq
  std::ofstream* qseqOfsVar = NULL;
  m_env.openOutputFile(m_qseqOutputFileName,
                       "hdf",
                       m_qseqOutputAllow,
                       true,
                       qseqOfsVar);

  if (qseqOfsVar) {
    workingQSeq.printContents(*qseqOfsVar);
  }

  if (qseqOfsVar) {
    qseqOfsVar->close();
    if (m_env.subScreenFile()) {
      *m_env.subScreenFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                             << ", prefix = "            << m_prefix
                             << ": closed output file '" << m_qseqOutputFileName
                             << "' for qseq "            << workingQSeq.name()
                             << std::endl;
    }
  }

  // Take "unified" care of qseq
#if 0
  std::ofstream* unifiedQSeqOfsVar = NULL;
  m_env.openUnifiedOutputFile(m_qseqOutputFileName,
                              "hdf",
                              true,
                              unifiedQSeqOfsVar);

  if (unifiedQSeqOfsVar) {
    workingQSeq.printUnifiedContents(*unifiedQSeqOfsVar);
  }

  if (unifiedQSeqOfsVar) {
    unifiedQSeqOfsVar->close();
    if (m_env.subScreenFile()) {
      *m_env.subScreenFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                             << ", prefix = "                    << m_prefix
                             << ": closed unified output file '" << m_qseqOutputFileName
                             << "' for qseq "                    << workingQSeq.name()
                             << std::endl;
    }
  }
#else
  if (m_qseqOutputFileName != ".") {
    workingQSeq.printUnifiedContents(m_qseqOutputFileName);
    if (m_env.subScreenFile()) {
      *m_env.subScreenFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                             << ", prefix = "                    << m_prefix
                             << ": closed unified output file '" << m_qseqOutputFileName
                             << "' for qseq "                    << workingQSeq.name()
                             << std::endl;
    }
  }
#endif

  // Take case of other aspects of qseq
  if (m_qseqComputeStats) {
    if ((m_env.subScreenFile()) && (m_env.verbosity() >= 0)) {
      *m_env.subScreenFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                             << ": about to call 'workingQSeq.computeStatistics()'"
                             << std::endl;
    }
    workingQSeq.computeStatistics(*m_qseqStatisticalOptions,
                                  genericOfsVar);
    if ((m_env.subScreenFile()) && (m_env.verbosity() >= 0)) {
      *m_env.subScreenFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                             << ": returned from call to 'workingQSeq.computeStatistics()'"
                             << std::endl;
    }
  }

  //****************************************************
  // Close generic output file      
  //****************************************************
  if (genericOfsVar) {
    genericOfsVar->close();
    if (m_env.subScreenFile()) {
      *m_env.subScreenFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                             << ", prefix = "                    << m_prefix
                             << ": closed generic output file '" << m_outputFileName
                             << "' for qoi sequence "            << workingQSeq.name()
                             << std::endl;
    }
  }
  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << std::endl;
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::actualGenerateSequence(
  const uqBaseVectorRVClass      <P_V,P_M>& paramRv,
        uqBaseVectorSequenceClass<P_V,P_M>& workingPSeq,
        uqBaseVectorSequenceClass<Q_V,Q_M>& workingQSeq,
        unsigned int                        seqSize)
{
  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "Starting the generation of qoi sequence " << workingQSeq.name()
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

  workingPSeq.resizeSequence(seqSize);
  workingQSeq.resizeSequence(seqSize);
  P_V tmpP(m_paramSpace.zeroVector());
  Q_V tmpQ(m_qoiSpace.zeroVector());

  for (unsigned int i = 0; i < seqSize; ++i) {
    paramRv.realizer().realization(tmpP);

    workingPSeq.setPositionValues(i,tmpP);

    if (m_qseqMeasureRunTimes) iRC = gettimeofday(&timevalQoIFunction, NULL);
    m_qoiFunctionSynchronizer->callFunction(&tmpP,NULL,&tmpQ,NULL,NULL,NULL); // Might demand parallel environment
    if (m_qseqMeasureRunTimes) qoiFunctionRunTime += uqMiscGetEllapsedSeconds(&timevalQoIFunction);

    workingQSeq.setPositionValues(i,tmpQ);

    if ((m_qseqDisplayPeriod            > 0) && 
        (((i+1) % m_qseqDisplayPeriod) == 0)) {
      if (m_env.subScreenFile()) {
        *m_env.subScreenFile() << "Finished generating " << i+1
                               << " qoi samples"
                               << std::endl;
      }
    }
  }

  seqRunTime = uqMiscGetEllapsedSeconds(&timevalSeq);

  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "Finished the generation of qoi sequence " << workingQSeq.name()
                           << ", with "                                  << workingQSeq.sequenceSize()
                           << " samples";
    *m_env.subScreenFile() << "\nSome information about this sequence:"
                           << "\n  Sequence run time = " << seqRunTime
                           << " seconds";
    *m_env.subScreenFile() << "\n\n Breaking of the seq run time:\n";
    *m_env.subScreenFile() << "\n  QoI function run time   = " << qoiFunctionRunTime
                           << " seconds ("                     << 100.*qoiFunctionRunTime/seqRunTime
                           << "%)";
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::checkTheParallelEnvironment()
{
  if (m_env.numSubEnvironments() == (unsigned int) m_env.fullComm().NumProc()) {
    UQ_FATAL_TEST_MACRO(m_env.subRank() != 0,
                        m_env.rank(),
                        "uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::checkTheParallelEnvironment()",
                        "there should exist only one processor per sub environment");
    UQ_FATAL_TEST_MACRO(m_paramRv.imageSet().vectorSpace().zeroVector().numberOfProcessorsRequiredForStorage() != 1,
                        m_env.rank(),
                        "uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::checkTheParallelEnvironment()",
                        "only 1 processor (per sub environment) should be necessary for the storage of a parameter vector");
  }
  else if (m_env.numSubEnvironments() < (unsigned int) m_env.fullComm().NumProc()) {
    UQ_FATAL_TEST_MACRO(m_env.fullComm().NumProc()%m_env.numSubEnvironments() != 0,
                        m_env.rank(),
                        "uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::checkTheParallelEnvironment()",
                        "total number of processors should be a multiple of the number of sub environments");
    unsigned int numProcsPerSubEnvironment = m_env.fullComm().NumProc()/m_env.numSubEnvironments();
    UQ_FATAL_TEST_MACRO(m_env.subComm().NumProc() != (int) numProcsPerSubEnvironment,
                        m_env.rank(),
                        "uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::checkTheParallelEnvironment()",
                        "inconsistent number of processors per sub environment");
    if ((m_paramRv.imageSet().vectorSpace().zeroVector().numberOfProcessorsRequiredForStorage() == 1) &&
        (m_qoiRv.imageSet().vectorSpace().zeroVector().numberOfProcessorsRequiredForStorage()   == 1)) {
      // Ok
    }
    else if ((m_paramRv.imageSet().vectorSpace().zeroVector().numberOfProcessorsRequiredForStorage() == numProcsPerSubEnvironment) &&
             (m_qoiRv.imageSet().vectorSpace().zeroVector().numberOfProcessorsRequiredForStorage()   == numProcsPerSubEnvironment)) {
      UQ_FATAL_TEST_MACRO(true,
                          m_env.rank(),
                          "uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::checkTheParallelEnvironment()",
                          "parallel vectors are not supported yet");
    }
    else {
      UQ_FATAL_TEST_MACRO(true,
                          m_env.rank(),
                          "uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::checkTheParallelEnvironment()",
                          "number of processors required for a vector storage should be equal to the number of processors in the sub environment");
    }
  }
  else {
    UQ_FATAL_TEST_MACRO(true,
                        m_env.rank(),
                        "uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::checkTheParallelEnvironment()",
                        "number of processors per sub environment is too large");
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::print(std::ostream& os) const
{
  os <<         m_option_outputFileName << " = " << m_outputFileName
     << "\n" << m_option_outputAllow    << " = ";
  for (std::set<unsigned int>::iterator setIt = m_outputAllow.begin(); setIt != m_outputAllow.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_pseq_outputFileName  << " = " << m_pseqOutputFileName
     << "\n" << m_option_pseq_outputAllow     << " = ";
  for (std::set<unsigned int>::iterator setIt = m_pseqOutputAllow.begin(); setIt != m_pseqOutputAllow.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_pseq_computeStats    << " = " << m_pseqComputeStats
     << "\n" << m_option_qseq_size            << " = " << m_qseqSize
     << "\n" << m_option_qseq_displayPeriod   << " = " << m_qseqDisplayPeriod
     << "\n" << m_option_qseq_measureRunTimes << " = " << m_qseqMeasureRunTimes
     << "\n" << m_option_qseq_outputFileName  << " = " << m_qseqOutputFileName
     << "\n" << m_option_qseq_outputAllow     << " = ";
  for (std::set<unsigned int>::iterator setIt = m_qseqOutputAllow.begin(); setIt != m_qseqOutputAllow.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_qseq_computeStats << " = " << m_qseqComputeStats;

  return;
}

template<class P_V,class P_M,class Q_V,class Q_M> 
std::ostream& operator<<(std::ostream& os, const uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_MOC_SG_H__
