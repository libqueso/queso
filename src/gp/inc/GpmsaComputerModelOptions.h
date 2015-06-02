//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
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

#ifndef UQ_GCM_OPTIONS_H
#define UQ_GCM_OPTIONS_H

#include <queso/Environment.h>
#include <queso/BoostInputOptionsParser.h>
#include <queso/SequenceStatisticalOptions.h>

#define UQ_GCM_FILENAME_FOR_NO_FILE "."

// _ODV = option default value
#define UQ_GCM_COMPUTE_SOLUTION_ODV                      1
#define UQ_GCM_CHECK_AGAINST_PREVIOUS_SAMPLE_ODV         1
#define UQ_GCM_DATA_OUTPUT_FILE_NAME_ODV                 UQ_GCM_FILENAME_FOR_NO_FILE
#define UQ_GCM_DATA_OUTPUT_ALLOW_ALL_ODV                 0
#define UQ_GCM_DATA_OUTPUT_ALLOWED_SET_ODV               ""
#define UQ_GCM_PRIOR_SEQ_NUM_SAMPLES_ODV                 0
#define UQ_GCM_PRIOR_SEQ_DATA_OUTPUT_FILE_NAME_ODV       UQ_GCM_FILENAME_FOR_NO_FILE
#define UQ_GCM_PRIOR_SEQ_DATA_OUTPUT_FILE_TYPE_ODV       "m"
#define UQ_GCM_PRIOR_SEQ_DATA_OUTPUT_ALLOW_ALL_ODV       0
#define UQ_GCM_PRIOR_SEQ_DATA_OUTPUT_ALLOWED_SET_ODV     ""
#define UQ_GCM_NUGGET_VALUE_FOR_BT_WY_B_ODV              0.
#define UQ_GCM_NUGGET_VALUE_FOR_BT_WY_B_INV_ODV          0.
#define UQ_GCM_FORM_C_MATRIX_ODV                         0
#define UQ_GCM_USE_TILDE_LOGIC_FOR_RANK_DEFFICIENT_C_ODV 0
#define UQ_GCM_PRED_LAG_ODV                              1
#define UQ_GCM_PRED_VUS_BY_SAMPLING_RVS_ODV              0
#define UQ_GCM_PRED_VUS_BY_SUMMING_RVS_ODV               1
#define UQ_GCM_PRED_VUS_AT_KEY_POINTS_ODV                0
#define UQ_GCM_PRED_WS_BY_SAMPLING_RVS_ODV               0
#define UQ_GCM_PRED_WS_BY_SUMMING_RVS_ODV                1
#define UQ_GCM_PRED_WS_AT_KEY_POINTS_ODV                 0

namespace boost {
  namespace program_options {
    class options_description;
  }
}

namespace QUESO {

class GcmOptionsValues
{
public:
  GcmOptionsValues            ();
  GcmOptionsValues(const BaseEnvironment * env, const char * prefix);
  GcmOptionsValues            (const GcmOptionsValues& src);
  GcmOptionsValues& operator= (const GcmOptionsValues& rhs);
  virtual ~GcmOptionsValues            ();

  std::string m_prefix;

  bool                   m_checkAgainstPreviousSample;
  std::string            m_dataOutputFileName;
  bool                   m_dataOutputAllowAll;
  std::set<unsigned int> m_dataOutputAllowedSet;
  unsigned int           m_priorSeqNumSamples;
  std::string            m_priorSeqDataOutputFileName;
  std::string            m_priorSeqDataOutputFileType;
  bool                   m_priorSeqDataOutputAllowAll;
  std::set<unsigned int> m_priorSeqDataOutputAllowedSet;
  double                 m_nuggetValueForBtWyB;
  double                 m_nuggetValueForBtWyBInv;
  double                 m_formCMatrix;
  bool                   m_useTildeLogicForRankDefficientC;
  unsigned int           m_predLag;
  bool                   m_predVUsBySamplingRVs;
  bool                   m_predVUsBySummingRVs;
  bool                   m_predVUsAtKeyPoints;
  bool                   m_predWsBySamplingRVs;
  bool                   m_predWsBySummingRVs;
  bool                   m_predWsAtKeyPoints;

  //MhOptionsValues m_mhOptionsValues;

private:
  BoostInputOptionsParser * m_parser;

  std::string                   m_option_help;
  std::string                   m_option_checkAgainstPreviousSample;
  std::string                   m_option_dataOutputFileName;
  std::string                   m_option_dataOutputAllowAll;
  std::string                   m_option_dataOutputAllowedSet;
  std::string                   m_option_priorSeqNumSamples;
  std::string                   m_option_priorSeqDataOutputFileName;
  std::string                   m_option_priorSeqDataOutputFileType;
  std::string                   m_option_priorSeqDataOutputAllowAll;
  std::string                   m_option_priorSeqDataOutputAllowedSet;
  std::string                   m_option_nuggetValueForBtWyB;
  std::string                   m_option_nuggetValueForBtWyBInv;
  std::string                   m_option_formCMatrix;
  std::string                   m_option_useTildeLogicForRankDefficientC;
  std::string                   m_option_predLag;
  std::string                   m_option_predVUsBySamplingRVs;
  std::string                   m_option_predVUsBySummingRVs;
  std::string                   m_option_predVUsAtKeyPoints;
  std::string                   m_option_predWsBySamplingRVs;
  std::string                   m_option_predWsBySummingRVs;
  std::string                   m_option_predWsAtKeyPoints;

  void copy(const GcmOptionsValues& src);

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  friend class GpmsaComputerModelOptions;
  SsOptionsValues m_alternativePriorSeqSsOptionsValues;
#endif
};

class GpmsaComputerModelOptions
{
public:
  GpmsaComputerModelOptions(const BaseEnvironment& env, const char* prefix);
  GpmsaComputerModelOptions(const BaseEnvironment& env, const char* prefix, const GcmOptionsValues& alternativeOptionsValues);
 ~GpmsaComputerModelOptions();

  void scanOptionsValues();
  void print            (std::ostream& os) const;

  GcmOptionsValues            m_ov;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  SequenceStatisticalOptions* m_priorSeqStatisticalOptionsObj;
  bool                               m_priorSeqStatOptsInstantiated;
#endif
  std::string                        m_prefix;

private:
  void   defineMyOptions  (boost::program_options::options_description& optionsDesc) const;
  void   getMyOptionValues(boost::program_options::options_description& optionsDesc);

  const BaseEnvironment& m_env;

  boost::program_options::options_description*      m_optionsDesc;
  std::string                   m_option_help;
  std::string                   m_option_checkAgainstPreviousSample;
  std::string                   m_option_dataOutputFileName;
  std::string                   m_option_dataOutputAllowAll;
  std::string                   m_option_dataOutputAllowedSet;
  std::string                   m_option_priorSeqNumSamples;
  std::string                   m_option_priorSeqDataOutputFileName;
  std::string                   m_option_priorSeqDataOutputFileType;
  std::string                   m_option_priorSeqDataOutputAllowAll;
  std::string                   m_option_priorSeqDataOutputAllowedSet;
  std::string                   m_option_nuggetValueForBtWyB;
  std::string                   m_option_nuggetValueForBtWyBInv;
  std::string                   m_option_formCMatrix;
  std::string                   m_option_useTildeLogicForRankDefficientC;
  std::string                   m_option_predLag;
  std::string                   m_option_predVUsBySamplingRVs;
  std::string                   m_option_predVUsBySummingRVs;
  std::string                   m_option_predVUsAtKeyPoints;
  std::string                   m_option_predWsBySamplingRVs;
  std::string                   m_option_predWsBySummingRVs;
  std::string                   m_option_predWsAtKeyPoints;
};

std::ostream& operator<<(std::ostream& os, const GpmsaComputerModelOptions& obj);

}  // End namespace QUESO

#endif // UQ_GCM_OPTIONS_H
