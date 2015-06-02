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

#include <boost/program_options.hpp>

#include <queso/GpmsaComputerModelOptions.h>
#include <queso/Miscellaneous.h>

namespace QUESO {

GcmOptionsValues::GcmOptionsValues()
  :
    m_prefix                                ("gcm_"),
    m_checkAgainstPreviousSample     (UQ_GCM_CHECK_AGAINST_PREVIOUS_SAMPLE_ODV        ),
    m_dataOutputFileName             (UQ_GCM_DATA_OUTPUT_FILE_NAME_ODV                ),
    m_dataOutputAllowAll             (UQ_GCM_DATA_OUTPUT_ALLOW_ALL_ODV                ),
    //m_dataOutputAllowedSet           (),
    m_priorSeqNumSamples             (UQ_GCM_PRIOR_SEQ_NUM_SAMPLES_ODV                ),
    m_priorSeqDataOutputFileName     (UQ_GCM_PRIOR_SEQ_DATA_OUTPUT_FILE_NAME_ODV      ),
    m_priorSeqDataOutputFileType     (UQ_GCM_PRIOR_SEQ_DATA_OUTPUT_FILE_TYPE_ODV      ),
    m_priorSeqDataOutputAllowAll     (UQ_GCM_PRIOR_SEQ_DATA_OUTPUT_ALLOW_ALL_ODV      ),
    //m_priorSeqDataOutputAllowedSet   (),
    m_nuggetValueForBtWyB            (UQ_GCM_NUGGET_VALUE_FOR_BT_WY_B_ODV             ),
    m_nuggetValueForBtWyBInv         (UQ_GCM_NUGGET_VALUE_FOR_BT_WY_B_INV_ODV         ),
    m_formCMatrix                    (UQ_GCM_FORM_C_MATRIX_ODV                        ),
    m_useTildeLogicForRankDefficientC(UQ_GCM_USE_TILDE_LOGIC_FOR_RANK_DEFFICIENT_C_ODV),
    m_predLag                        (UQ_GCM_PRED_LAG_ODV                             ),
    m_predVUsBySamplingRVs           (UQ_GCM_PRED_VUS_BY_SAMPLING_RVS_ODV             ),
    m_predVUsBySummingRVs            (UQ_GCM_PRED_VUS_BY_SUMMING_RVS_ODV              ),
    m_predVUsAtKeyPoints             (UQ_GCM_PRED_VUS_AT_KEY_POINTS_ODV               ),
    m_predWsBySamplingRVs            (UQ_GCM_PRED_WS_BY_SAMPLING_RVS_ODV              ),
    m_predWsBySummingRVs             (UQ_GCM_PRED_WS_BY_SUMMING_RVS_ODV               ),
    m_predWsAtKeyPoints              (UQ_GCM_PRED_WS_AT_KEY_POINTS_ODV                ),
    m_parser(NULL),
    m_option_help                           (m_prefix + "help"                           ),
    m_option_checkAgainstPreviousSample     (m_prefix + "checkAgainstPreviousSample"     ),
    m_option_dataOutputFileName             (m_prefix + "dataOutputFileName"             ),
    m_option_dataOutputAllowAll             (m_prefix + "dataOutputAllowAll"             ),
    m_option_dataOutputAllowedSet           (m_prefix + "dataOutputAllowedSet"           ),
    m_option_priorSeqNumSamples             (m_prefix + "priorSeqNumSamples"             ),
    m_option_priorSeqDataOutputFileName     (m_prefix + "priorSeqDataOutputFileName"     ),
    m_option_priorSeqDataOutputFileType     (m_prefix + "priorSeqDataOutputFileType"     ),
    m_option_priorSeqDataOutputAllowAll     (m_prefix + "priorSeqDataOutputAllowAll"     ),
    m_option_priorSeqDataOutputAllowedSet   (m_prefix + "priorSeqDataOutputAllowedSet"   ),
    m_option_nuggetValueForBtWyB            (m_prefix + "nuggetValueForBtWyB"            ),
    m_option_nuggetValueForBtWyBInv         (m_prefix + "nuggetValueForBtWyBInv"         ),
    m_option_formCMatrix                    (m_prefix + "formCMatrix"                    ),
    m_option_useTildeLogicForRankDefficientC(m_prefix + "useTildeLogicForRankDefficientC"),
    m_option_predLag                        (m_prefix + "predLag"                        ),
    m_option_predVUsBySamplingRVs           (m_prefix + "predVUsBySamplingRVs"           ),
    m_option_predVUsBySummingRVs            (m_prefix + "predVUsBySummingRVs"            ),
    m_option_predVUsAtKeyPoints             (m_prefix + "predVUsAtKeyPoints"             ),
    m_option_predWsBySamplingRVs            (m_prefix + "predWsBySamplingRVs"            ),
    m_option_predWsBySummingRVs             (m_prefix + "predWsBySummingRVs"             ),
    m_option_predWsAtKeyPoints              (m_prefix + "predWsAtKeyPoints"              )
{
}

GcmOptionsValues::GcmOptionsValues(const BaseEnvironment * env, const char *
    prefix)
  :
    m_prefix                                ((std::string)(prefix) + "gcm_"),
    m_checkAgainstPreviousSample     (UQ_GCM_CHECK_AGAINST_PREVIOUS_SAMPLE_ODV        ),
    m_dataOutputFileName             (UQ_GCM_DATA_OUTPUT_FILE_NAME_ODV                ),
    m_dataOutputAllowAll             (UQ_GCM_DATA_OUTPUT_ALLOW_ALL_ODV                ),
    //m_dataOutputAllowedSet           (),
    m_priorSeqNumSamples             (UQ_GCM_PRIOR_SEQ_NUM_SAMPLES_ODV                ),
    m_priorSeqDataOutputFileName     (UQ_GCM_PRIOR_SEQ_DATA_OUTPUT_FILE_NAME_ODV      ),
    m_priorSeqDataOutputFileType     (UQ_GCM_PRIOR_SEQ_DATA_OUTPUT_FILE_TYPE_ODV      ),
    m_priorSeqDataOutputAllowAll     (UQ_GCM_PRIOR_SEQ_DATA_OUTPUT_ALLOW_ALL_ODV      ),
    //m_priorSeqDataOutputAllowedSet   (),
    m_nuggetValueForBtWyB            (UQ_GCM_NUGGET_VALUE_FOR_BT_WY_B_ODV             ),
    m_nuggetValueForBtWyBInv         (UQ_GCM_NUGGET_VALUE_FOR_BT_WY_B_INV_ODV         ),
    m_formCMatrix                    (UQ_GCM_FORM_C_MATRIX_ODV                        ),
    m_useTildeLogicForRankDefficientC(UQ_GCM_USE_TILDE_LOGIC_FOR_RANK_DEFFICIENT_C_ODV),
    m_predLag                        (UQ_GCM_PRED_LAG_ODV                             ),
    m_predVUsBySamplingRVs           (UQ_GCM_PRED_VUS_BY_SAMPLING_RVS_ODV             ),
    m_predVUsBySummingRVs            (UQ_GCM_PRED_VUS_BY_SUMMING_RVS_ODV              ),
    m_predVUsAtKeyPoints             (UQ_GCM_PRED_VUS_AT_KEY_POINTS_ODV               ),
    m_predWsBySamplingRVs            (UQ_GCM_PRED_WS_BY_SAMPLING_RVS_ODV              ),
    m_predWsBySummingRVs             (UQ_GCM_PRED_WS_BY_SUMMING_RVS_ODV               ),
    m_predWsAtKeyPoints              (UQ_GCM_PRED_WS_AT_KEY_POINTS_ODV                ),
    m_parser(new BoostInputOptionsParser(env->optionsInputFileName())),
    m_option_help                           (m_prefix + "help"                           ),
    m_option_checkAgainstPreviousSample     (m_prefix + "checkAgainstPreviousSample"     ),
    m_option_dataOutputFileName             (m_prefix + "dataOutputFileName"             ),
    m_option_dataOutputAllowAll             (m_prefix + "dataOutputAllowAll"             ),
    m_option_dataOutputAllowedSet           (m_prefix + "dataOutputAllowedSet"           ),
    m_option_priorSeqNumSamples             (m_prefix + "priorSeqNumSamples"             ),
    m_option_priorSeqDataOutputFileName     (m_prefix + "priorSeqDataOutputFileName"     ),
    m_option_priorSeqDataOutputFileType     (m_prefix + "priorSeqDataOutputFileType"     ),
    m_option_priorSeqDataOutputAllowAll     (m_prefix + "priorSeqDataOutputAllowAll"     ),
    m_option_priorSeqDataOutputAllowedSet   (m_prefix + "priorSeqDataOutputAllowedSet"   ),
    m_option_nuggetValueForBtWyB            (m_prefix + "nuggetValueForBtWyB"            ),
    m_option_nuggetValueForBtWyBInv         (m_prefix + "nuggetValueForBtWyBInv"         ),
    m_option_formCMatrix                    (m_prefix + "formCMatrix"                    ),
    m_option_useTildeLogicForRankDefficientC(m_prefix + "useTildeLogicForRankDefficientC"),
    m_option_predLag                        (m_prefix + "predLag"                        ),
    m_option_predVUsBySamplingRVs           (m_prefix + "predVUsBySamplingRVs"           ),
    m_option_predVUsBySummingRVs            (m_prefix + "predVUsBySummingRVs"            ),
    m_option_predVUsAtKeyPoints             (m_prefix + "predVUsAtKeyPoints"             ),
    m_option_predWsBySamplingRVs            (m_prefix + "predWsBySamplingRVs"            ),
    m_option_predWsBySummingRVs             (m_prefix + "predWsBySummingRVs"             ),
    m_option_predWsAtKeyPoints              (m_prefix + "predWsAtKeyPoints"              )
{
  m_parser->registerOption(m_option_help.c_str(),                                                                                                                        "produce help message for mixed inverse problem");
  m_parser->registerOption<bool        >(m_option_checkAgainstPreviousSample,      UQ_GCM_CHECK_AGAINST_PREVIOUS_SAMPLE_ODV        , "check against previous sample"                 );
  m_parser->registerOption<std::string >(m_option_dataOutputFileName,              UQ_GCM_DATA_OUTPUT_FILE_NAME_ODV                , "name of data output file"                      );
  m_parser->registerOption<bool        >(m_option_dataOutputAllowAll,              UQ_GCM_DATA_OUTPUT_ALLOW_ALL_ODV                , "allow all or not"                              );
  m_parser->registerOption<std::string >(m_option_dataOutputAllowedSet,            UQ_GCM_DATA_OUTPUT_ALLOWED_SET_ODV              , "subEnvs that will write to data output file"   );
  m_parser->registerOption<unsigned int>(m_option_priorSeqNumSamples,              UQ_GCM_PRIOR_SEQ_NUM_SAMPLES_ODV                , "prior sequence size"                           );
  m_parser->registerOption<std::string >(m_option_priorSeqDataOutputFileName,      UQ_GCM_PRIOR_SEQ_DATA_OUTPUT_FILE_NAME_ODV      , "prior sequence data output filename"           );
  m_parser->registerOption<std::string >(m_option_priorSeqDataOutputFileType,      UQ_GCM_PRIOR_SEQ_DATA_OUTPUT_FILE_TYPE_ODV      , "prior sequence data output filetype"           );
  m_parser->registerOption<bool        >(m_option_priorSeqDataOutputAllowAll,      UQ_GCM_PRIOR_SEQ_DATA_OUTPUT_ALLOW_ALL_ODV      , "allow all or not"                              );
  m_parser->registerOption<std::string >(m_option_priorSeqDataOutputAllowedSet,    UQ_GCM_PRIOR_SEQ_DATA_OUTPUT_ALLOWED_SET_ODV    , "subEnvs that will write to data output file"   );
  m_parser->registerOption<double      >(m_option_nuggetValueForBtWyB,             UQ_GCM_NUGGET_VALUE_FOR_BT_WY_B_ODV             , "nugget value for Bt_Wy_W matrix"               );
  m_parser->registerOption<double      >(m_option_nuggetValueForBtWyBInv,          UQ_GCM_NUGGET_VALUE_FOR_BT_WY_B_INV_ODV         , "nugget value for Bt_Wy_W inverse matrix"       );
  m_parser->registerOption<double      >(m_option_formCMatrix,                     UQ_GCM_FORM_C_MATRIX_ODV                        , "form C matrix"                                 );
  m_parser->registerOption<bool        >(m_option_useTildeLogicForRankDefficientC, UQ_GCM_USE_TILDE_LOGIC_FOR_RANK_DEFFICIENT_C_ODV, "use tilde logic for rank defficient C"         );
  m_parser->registerOption<unsigned int>(m_option_predLag,                         UQ_GCM_PRED_LAG_ODV                             , "predLag"                                       );
  m_parser->registerOption<bool        >(m_option_predVUsBySamplingRVs,            UQ_GCM_PRED_VUS_BY_SAMPLING_RVS_ODV             , "predVUsBySamplingRVs"                          );
  m_parser->registerOption<bool        >(m_option_predVUsBySummingRVs,             UQ_GCM_PRED_VUS_BY_SUMMING_RVS_ODV              , "predVUsBySummingRVs"                           );
  m_parser->registerOption<bool        >(m_option_predVUsAtKeyPoints,              UQ_GCM_PRED_VUS_AT_KEY_POINTS_ODV               , "predVUsAtKeyPoints"                            );
  m_parser->registerOption<bool        >(m_option_predWsBySamplingRVs,             UQ_GCM_PRED_WS_BY_SAMPLING_RVS_ODV              , "predWsBySamplingRVs"                           );
  m_parser->registerOption<bool        >(m_option_predWsBySummingRVs,              UQ_GCM_PRED_WS_BY_SUMMING_RVS_ODV               , "predWsBySummingRVs"                            );
  m_parser->registerOption<bool        >(m_option_predWsAtKeyPoints,               UQ_GCM_PRED_WS_AT_KEY_POINTS_ODV                , "predWsAtKeyPoints"                             );

  m_parser->scanInputFile();

  m_parser->getOption<bool        >(m_option_checkAgainstPreviousSample,      m_checkAgainstPreviousSample);
  m_parser->getOption<std::string >(m_option_dataOutputFileName,              m_dataOutputFileName);
  m_parser->getOption<bool        >(m_option_dataOutputAllowAll,              m_dataOutputAllowAll);
  m_parser->getOption<std::set<unsigned int> >(m_option_dataOutputAllowedSet,            m_dataOutputAllowedSet);
  m_parser->getOption<unsigned int>(m_option_priorSeqNumSamples,              m_priorSeqNumSamples);
  m_parser->getOption<std::string >(m_option_priorSeqDataOutputFileName,      m_priorSeqDataOutputFileName);
  m_parser->getOption<std::string >(m_option_priorSeqDataOutputFileType,      m_priorSeqDataOutputFileType);
  m_parser->getOption<bool        >(m_option_priorSeqDataOutputAllowAll,      m_priorSeqDataOutputAllowAll);
  m_parser->getOption<std::set<unsigned int> >(m_option_priorSeqDataOutputAllowedSet,    m_priorSeqDataOutputAllowedSet);
  m_parser->getOption<double      >(m_option_nuggetValueForBtWyB,             m_nuggetValueForBtWyB);
  m_parser->getOption<double      >(m_option_nuggetValueForBtWyBInv,          m_nuggetValueForBtWyBInv);
  m_parser->getOption<double      >(m_option_formCMatrix,                     m_formCMatrix);
  m_parser->getOption<bool        >(m_option_useTildeLogicForRankDefficientC, m_useTildeLogicForRankDefficientC);
  m_parser->getOption<unsigned int>(m_option_predLag,                         m_predLag);
  m_parser->getOption<bool        >(m_option_predVUsBySamplingRVs,            m_predVUsBySamplingRVs);
  m_parser->getOption<bool        >(m_option_predVUsBySummingRVs,             m_predVUsBySummingRVs);
  m_parser->getOption<bool        >(m_option_predVUsAtKeyPoints,              m_predVUsAtKeyPoints);
  m_parser->getOption<bool        >(m_option_predWsBySamplingRVs,             m_predWsBySamplingRVs);
  m_parser->getOption<bool        >(m_option_predWsBySummingRVs,              m_predWsBySummingRVs);
  m_parser->getOption<bool        >(m_option_predWsAtKeyPoints,               m_predWsAtKeyPoints);
}

GcmOptionsValues::~GcmOptionsValues()
{
}

GcmOptionsValues::GcmOptionsValues(const GcmOptionsValues& src)
{
  this->copy(src);
}

GcmOptionsValues&
GcmOptionsValues::operator=(const GcmOptionsValues& rhs)
{
  this->copy(rhs);
  return *this;
}

// void
// GcmOptionsValues::defineOptions()
// {
//   (*m_optionsDescription).add_options()
//     (m_option_help.c_str(),                                                                                                                        "produce help message for mixed inverse problem")
//     (m_option_checkAgainstPreviousSample.c_str(),      boost::program_options::value<bool        >()->default_value(UQ_GCM_CHECK_AGAINST_PREVIOUS_SAMPLE_ODV        ), "check against previous sample"                 )
//     (m_option_dataOutputFileName.c_str(),              boost::program_options::value<std::string >()->default_value(UQ_GCM_DATA_OUTPUT_FILE_NAME_ODV                ), "name of data output file"                      )
//     (m_option_dataOutputAllowAll.c_str(),              boost::program_options::value<bool        >()->default_value(UQ_GCM_DATA_OUTPUT_ALLOW_ALL_ODV                ), "allow all or not"                              )
//     (m_option_dataOutputAllowedSet.c_str(),            boost::program_options::value<std::string >()->default_value(UQ_GCM_DATA_OUTPUT_ALLOWED_SET_ODV              ), "subEnvs that will write to data output file"   )
//     (m_option_priorSeqNumSamples.c_str(),              boost::program_options::value<unsigned int>()->default_value(UQ_GCM_PRIOR_SEQ_NUM_SAMPLES_ODV                ), "prior sequence size"                           )
//     (m_option_priorSeqDataOutputFileName.c_str(),      boost::program_options::value<std::string >()->default_value(UQ_GCM_PRIOR_SEQ_DATA_OUTPUT_FILE_NAME_ODV      ), "prior sequence data output filename"           )
//     (m_option_priorSeqDataOutputFileType.c_str(),      boost::program_options::value<std::string >()->default_value(UQ_GCM_PRIOR_SEQ_DATA_OUTPUT_FILE_TYPE_ODV      ), "prior sequence data output filetype"           )
//     (m_option_priorSeqDataOutputAllowAll.c_str(),      boost::program_options::value<bool        >()->default_value(UQ_GCM_PRIOR_SEQ_DATA_OUTPUT_ALLOW_ALL_ODV      ), "allow all or not"                              )
//     (m_option_priorSeqDataOutputAllowedSet.c_str(),    boost::program_options::value<std::string >()->default_value(UQ_GCM_PRIOR_SEQ_DATA_OUTPUT_ALLOWED_SET_ODV    ), "subEnvs that will write to data output file"   )
//     (m_option_nuggetValueForBtWyB.c_str(),             boost::program_options::value<double      >()->default_value(UQ_GCM_NUGGET_VALUE_FOR_BT_WY_B_ODV             ), "nugget value for Bt_Wy_W matrix"               )
//     (m_option_nuggetValueForBtWyBInv.c_str(),          boost::program_options::value<double      >()->default_value(UQ_GCM_NUGGET_VALUE_FOR_BT_WY_B_INV_ODV         ), "nugget value for Bt_Wy_W inverse matrix"       )
//     (m_option_formCMatrix.c_str(),                     boost::program_options::value<double      >()->default_value(UQ_GCM_FORM_C_MATRIX_ODV                        ), "form C matrix"                                 )
//     (m_option_useTildeLogicForRankDefficientC.c_str(), boost::program_options::value<bool        >()->default_value(UQ_GCM_USE_TILDE_LOGIC_FOR_RANK_DEFFICIENT_C_ODV), "use tilde logic for rank defficient C"         )
//     (m_option_predLag.c_str(),                         boost::program_options::value<unsigned int>()->default_value(UQ_GCM_PRED_LAG_ODV                             ), "predLag"                                       )
//     (m_option_predVUsBySamplingRVs.c_str(),            boost::program_options::value<bool        >()->default_value(UQ_GCM_PRED_VUS_BY_SAMPLING_RVS_ODV             ), "predVUsBySamplingRVs"                          )
//     (m_option_predVUsBySummingRVs.c_str(),             boost::program_options::value<bool        >()->default_value(UQ_GCM_PRED_VUS_BY_SUMMING_RVS_ODV              ), "predVUsBySummingRVs"                           )
//     (m_option_predVUsAtKeyPoints.c_str(),              boost::program_options::value<bool        >()->default_value(UQ_GCM_PRED_VUS_AT_KEY_POINTS_ODV               ), "predVUsAtKeyPoints"                            )
//     (m_option_predWsBySamplingRVs.c_str(),             boost::program_options::value<bool        >()->default_value(UQ_GCM_PRED_WS_BY_SAMPLING_RVS_ODV              ), "predWsBySamplingRVs"                           )
//     (m_option_predWsBySummingRVs.c_str(),              boost::program_options::value<bool        >()->default_value(UQ_GCM_PRED_WS_BY_SUMMING_RVS_ODV               ), "predWsBySummingRVs"                            )
//     (m_option_predWsAtKeyPoints.c_str(),               boost::program_options::value<bool        >()->default_value(UQ_GCM_PRED_WS_AT_KEY_POINTS_ODV                ), "predWsAtKeyPoints"                             )
//   ;
// }
//
// void
// GcmOptionsValues::getOptionValues()
// {
//   if ((*m_optionsMap).count(m_option_help)) {
//     if (m_env->subDisplayFile()) {
//       *m_env->subDisplayFile() << (*m_optionsDescription)
//                               << std::endl;
//     }
//   }
//
//   if ((*m_optionsMap).count(m_option_checkAgainstPreviousSample)) {
//     m_checkAgainstPreviousSample = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_checkAgainstPreviousSample]).as<bool>();
//   }
//
//   if ((*m_optionsMap).count(m_option_dataOutputFileName)) {
//     m_dataOutputFileName = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_dataOutputFileName]).as<std::string>();
//   }
//
//   if ((*m_optionsMap).count(m_option_dataOutputAllowAll)) {
//     m_dataOutputAllowAll = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_dataOutputAllowAll]).as<bool>();
//   }
//
//   if (m_dataOutputAllowAll) {
//     m_dataOutputAllowedSet.insert(m_env->subId());
//   }
//   else if ((*m_optionsMap).count(m_option_dataOutputAllowedSet)) {
//     m_dataOutputAllowedSet.clear();
//     std::vector<double> tmpAllow(0,0.);
//     std::string inputString = (*m_optionsMap)[m_option_dataOutputAllowedSet].as<std::string>();
//     MiscReadDoublesFromString(inputString,tmpAllow);
//
//     if (tmpAllow.size() > 0) {
//       for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
//         m_dataOutputAllowedSet.insert((unsigned int) tmpAllow[i]);
//       }
//     }
//   }
//
//   //std::cout << "In GpmsaComputerModelOptions::getMyOptionValues(), m_option_priorSeqNumSamples = " << m_option_priorSeqNumSamples << "___" << std::endl;
//   if ((*m_optionsMap).count(m_option_priorSeqNumSamples)) {
//     //std::cout << "In GpmsaComputerModelOptions::getMyOptionValues(), going to read m_option_priorSeqNumSamples..." << std::endl;
//     m_priorSeqNumSamples = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_priorSeqNumSamples]).as<unsigned int>();
//     //std::cout << "In GpmsaComputerModelOptions::getMyOptionValues(), just read m_option_priorSeqNumSamples = " << m_priorSeqNumSamples << std::endl;
//   }
//
//   if ((*m_optionsMap).count(m_option_priorSeqDataOutputFileName)) {
//     m_priorSeqDataOutputFileName = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_priorSeqDataOutputFileName]).as<std::string>();
//   }
//
//   if ((*m_optionsMap).count(m_option_priorSeqDataOutputFileType)) {
//     m_priorSeqDataOutputFileType = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_priorSeqDataOutputFileType]).as<std::string>();
//   }
//
//   if ((*m_optionsMap).count(m_option_priorSeqDataOutputAllowAll)) {
//     m_priorSeqDataOutputAllowAll = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_priorSeqDataOutputAllowAll]).as<bool>();
//   }
//
//   if (m_priorSeqDataOutputAllowAll) {
//     m_priorSeqDataOutputAllowedSet.insert(m_env->subId());
//   }
//   else if ((*m_optionsMap).count(m_option_priorSeqDataOutputAllowedSet)) {
//     m_priorSeqDataOutputAllowedSet.clear();
//     std::vector<double> tmpAllow(0,0.);
//     std::string inputString = (*m_optionsMap)[m_option_priorSeqDataOutputAllowedSet].as<std::string>();
//     MiscReadDoublesFromString(inputString,tmpAllow);
//
//     if (tmpAllow.size() > 0) {
//       for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
//         m_priorSeqDataOutputAllowedSet.insert((unsigned int) tmpAllow[i]);
//       }
//     }
//   }
//
//   if ((*m_optionsMap).count(m_option_nuggetValueForBtWyB)) {
//     m_nuggetValueForBtWyB = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_nuggetValueForBtWyB]).as<double>();
//   }
//
//   if ((*m_optionsMap).count(m_option_nuggetValueForBtWyBInv)) {
//     m_nuggetValueForBtWyBInv = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_nuggetValueForBtWyBInv]).as<double>();
//   }
//
//   if ((*m_optionsMap).count(m_option_formCMatrix)) {
//     m_formCMatrix = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_formCMatrix]).as<double>();
//   }
//
//   if ((*m_optionsMap).count(m_option_useTildeLogicForRankDefficientC)) {
//     m_useTildeLogicForRankDefficientC = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_useTildeLogicForRankDefficientC]).as<bool>();
//   }
//
//   if ((*m_optionsMap).count(m_option_predLag)) {
//     m_predLag = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_predLag]).as<unsigned int>();
//   }
//
//   if ((*m_optionsMap).count(m_option_predVUsBySamplingRVs)) {
//     m_predVUsBySamplingRVs = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_predVUsBySamplingRVs]).as<bool>();
//   }
//
//   if ((*m_optionsMap).count(m_option_predVUsBySummingRVs)) {
//     m_predVUsBySummingRVs = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_predVUsBySummingRVs]).as<bool>();
//   }
//
//   if ((*m_optionsMap).count(m_option_predVUsAtKeyPoints)) {
//     m_predVUsAtKeyPoints = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_predVUsAtKeyPoints]).as<bool>();
//   }
//
//   if ((*m_optionsMap).count(m_option_predWsBySamplingRVs)) {
//     m_predWsBySamplingRVs = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_predWsBySamplingRVs]).as<bool>();
//   }
//
//   if ((*m_optionsMap).count(m_option_predWsBySummingRVs)) {
//     m_predWsBySummingRVs = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_predWsBySummingRVs]).as<bool>();
//   }
//
//   if ((*m_optionsMap).count(m_option_predWsAtKeyPoints)) {
//     m_predWsAtKeyPoints = ((const boost::program_options::variable_value&) (*m_optionsMap)[m_option_predWsAtKeyPoints]).as<bool>();
//   }
// }

void
GcmOptionsValues::copy(const GcmOptionsValues& src)
{
  m_checkAgainstPreviousSample      = src.m_checkAgainstPreviousSample;
  m_dataOutputFileName              = src.m_dataOutputFileName;
  m_dataOutputAllowAll              = src.m_dataOutputAllowAll;
  m_dataOutputAllowedSet            = src.m_dataOutputAllowedSet;
  m_priorSeqNumSamples              = src.m_priorSeqNumSamples;
  m_priorSeqDataOutputFileName      = src.m_priorSeqDataOutputFileName;
  m_priorSeqDataOutputFileType      = src.m_priorSeqDataOutputFileType;
  m_priorSeqDataOutputAllowAll      = src.m_priorSeqDataOutputAllowAll;
  m_priorSeqDataOutputAllowedSet    = src.m_priorSeqDataOutputAllowedSet;
  m_nuggetValueForBtWyB             = src.m_nuggetValueForBtWyB;
  m_nuggetValueForBtWyBInv          = src.m_nuggetValueForBtWyBInv;
  m_formCMatrix                     = src.m_formCMatrix;
  m_useTildeLogicForRankDefficientC = src.m_useTildeLogicForRankDefficientC;
  m_predLag                         = src.m_predLag;
  m_predVUsBySamplingRVs            = src.m_predVUsBySamplingRVs;
  m_predVUsBySummingRVs             = src.m_predVUsBySummingRVs;
  m_predVUsAtKeyPoints              = src.m_predVUsAtKeyPoints;
  m_predWsBySamplingRVs             = src.m_predWsBySamplingRVs;
  m_predWsBySummingRVs              = src.m_predWsBySummingRVs;
  m_predWsAtKeyPoints               = src.m_predWsAtKeyPoints;

  return;
}

GpmsaComputerModelOptions::GpmsaComputerModelOptions(
  const BaseEnvironment& env,
  const char*                   prefix)
  :
  m_ov                                    (),
  m_prefix                                ((std::string)(prefix) + "gcm_"),
  m_env                                   (env),
  m_optionsDesc                           (new boost::program_options::options_description("Mixed Inverse Problem options")),
  m_option_help                           (m_prefix + "help"                           ),
  m_option_checkAgainstPreviousSample     (m_prefix + "checkAgainstPreviousSample"     ),
  m_option_dataOutputFileName             (m_prefix + "dataOutputFileName"             ),
  m_option_dataOutputAllowAll             (m_prefix + "dataOutputAllowAll"             ),
  m_option_dataOutputAllowedSet           (m_prefix + "dataOutputAllowedSet"           ),
  m_option_priorSeqNumSamples             (m_prefix + "priorSeqNumSamples"             ),
  m_option_priorSeqDataOutputFileName     (m_prefix + "priorSeqDataOutputFileName"     ),
  m_option_priorSeqDataOutputFileType     (m_prefix + "priorSeqDataOutputFileType"     ),
  m_option_priorSeqDataOutputAllowAll     (m_prefix + "priorSeqDataOutputAllowAll"     ),
  m_option_priorSeqDataOutputAllowedSet   (m_prefix + "priorSeqDataOutputAllowedSet"   ),
  m_option_nuggetValueForBtWyB            (m_prefix + "nuggetValueForBtWyB"            ),
  m_option_nuggetValueForBtWyBInv         (m_prefix + "nuggetValueForBtWyBInv"         ),
  m_option_formCMatrix                    (m_prefix + "formCMatrix"                    ),
  m_option_useTildeLogicForRankDefficientC(m_prefix + "useTildeLogicForRankDefficientC"),
  m_option_predLag                        (m_prefix + "predLag"                        ),
  m_option_predVUsBySamplingRVs           (m_prefix + "predVUsBySamplingRVs"           ),
  m_option_predVUsBySummingRVs            (m_prefix + "predVUsBySummingRVs"            ),
  m_option_predVUsAtKeyPoints             (m_prefix + "predVUsAtKeyPoints"             ),
  m_option_predWsBySamplingRVs            (m_prefix + "predWsBySamplingRVs"            ),
  m_option_predWsBySummingRVs             (m_prefix + "predWsBySummingRVs"             ),
  m_option_predWsAtKeyPoints              (m_prefix + "predWsAtKeyPoints"              )
{
  queso_deprecated();
  queso_require_not_equal_to_msg(m_env.optionsInputFileName(), "", "this constructor is incompatible with the abscense of an options input file");
}

GpmsaComputerModelOptions::GpmsaComputerModelOptions(
  const BaseEnvironment&  env,
  const char*                    prefix,
  const GcmOptionsValues& alternativeOptionsValues)
  :
  m_ov                                    (alternativeOptionsValues),
  m_prefix                                ((std::string)(prefix) + "gcm_"),
  m_env                                   (env),
  m_optionsDesc                           (NULL),
  m_option_help                           (m_prefix + "help"                           ),
  m_option_checkAgainstPreviousSample     (m_prefix + "checkAgainstPreviousSample"     ),
  m_option_dataOutputFileName             (m_prefix + "dataOutputFileName"             ),
  m_option_dataOutputAllowAll             (m_prefix + "dataOutputAllowAll"             ),
  m_option_dataOutputAllowedSet           (m_prefix + "dataOutputAllowedSet"           ),
  m_option_priorSeqNumSamples             (m_prefix + "priorSeqNumSamples"             ),
  m_option_priorSeqDataOutputFileName     (m_prefix + "priorSeqDataOutputFileName"     ),
  m_option_priorSeqDataOutputFileType     (m_prefix + "priorSeqDataOutputFileType"     ),
  m_option_priorSeqDataOutputAllowAll     (m_prefix + "priorSeqDataOutputAllowAll"     ),
  m_option_priorSeqDataOutputAllowedSet   (m_prefix + "priorSeqDataOutputAllowedSet"   ),
  m_option_nuggetValueForBtWyB            (m_prefix + "nuggetValueForBtWyB"            ),
  m_option_nuggetValueForBtWyBInv         (m_prefix + "nuggetValueForBtWyBInv"         ),
  m_option_formCMatrix                    (m_prefix + "formCMatrix"                    ),
  m_option_useTildeLogicForRankDefficientC(m_prefix + "useTildeLogicForRankDefficientC"),
  m_option_predLag                        (m_prefix + "predLag"                        ),
  m_option_predVUsBySamplingRVs           (m_prefix + "predVUsBySamplingRVs"           ),
  m_option_predVUsBySummingRVs            (m_prefix + "predVUsBySummingRVs"            ),
  m_option_predVUsAtKeyPoints             (m_prefix + "predVUsAtKeyPoints"             ),
  m_option_predWsBySamplingRVs            (m_prefix + "predWsBySamplingRVs"            ),
  m_option_predWsBySummingRVs             (m_prefix + "predWsBySummingRVs"             ),
  m_option_predWsAtKeyPoints              (m_prefix + "predWsAtKeyPoints"              )
{
  queso_deprecated();
  queso_require_equal_to_msg(m_env.optionsInputFileName(), "", "this constructor is incompatible with the existence of an options input file");

  if (m_env.subDisplayFile() != NULL) {
    *m_env.subDisplayFile() << "In GpmsaComputerModelOptions::constructor(2)"
                            << ": after setting values of options with prefix '" << m_prefix
                            << "', state of object is:"
                            << "\n" << *this
                            << std::endl;
  }
}

GpmsaComputerModelOptions::~GpmsaComputerModelOptions()
{
  queso_deprecated();

  if (m_optionsDesc) delete m_optionsDesc;
}

void
GpmsaComputerModelOptions::scanOptionsValues()
{
  queso_deprecated();
  queso_require_msg(m_optionsDesc, "m_optionsDesc variable is NULL");

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  //std::cout << "scan 000\n"
  //          << std::endl;
  getMyOptionValues              (*m_optionsDesc);
  //std::cout << "scan 001\n"
  //          << std::endl;

  if (m_env.subDisplayFile() != NULL) {
    *m_env.subDisplayFile() << "In GpmsaComputerModelOptions::scanOptionsValues()"
                            << ": after reading values of options with prefix '" << m_prefix
                            << "', state of  object is:"
                            << "\n" << *this
                            << std::endl;
  }

  return;
}

void
GpmsaComputerModelOptions::defineMyOptions(boost::program_options::options_description& optionsDesc) const
{
  queso_deprecated();

  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                                                                        "produce help message for mixed inverse problem")
    (m_option_checkAgainstPreviousSample.c_str(),      boost::program_options::value<bool        >()->default_value(UQ_GCM_CHECK_AGAINST_PREVIOUS_SAMPLE_ODV        ), "check against previous sample"                 )
    (m_option_dataOutputFileName.c_str(),              boost::program_options::value<std::string >()->default_value(UQ_GCM_DATA_OUTPUT_FILE_NAME_ODV                ), "name of data output file"                      )
    (m_option_dataOutputAllowAll.c_str(),              boost::program_options::value<bool        >()->default_value(UQ_GCM_DATA_OUTPUT_ALLOW_ALL_ODV                ), "allow all or not"                              )
    (m_option_dataOutputAllowedSet.c_str(),            boost::program_options::value<std::string >()->default_value(UQ_GCM_DATA_OUTPUT_ALLOWED_SET_ODV              ), "subEnvs that will write to data output file"   )
    (m_option_priorSeqNumSamples.c_str(),              boost::program_options::value<unsigned int>()->default_value(UQ_GCM_PRIOR_SEQ_NUM_SAMPLES_ODV                ), "prior sequence size"                           )
    (m_option_priorSeqDataOutputFileName.c_str(),      boost::program_options::value<std::string >()->default_value(UQ_GCM_PRIOR_SEQ_DATA_OUTPUT_FILE_NAME_ODV      ), "prior sequence data output filename"           )
    (m_option_priorSeqDataOutputFileType.c_str(),      boost::program_options::value<std::string >()->default_value(UQ_GCM_PRIOR_SEQ_DATA_OUTPUT_FILE_TYPE_ODV      ), "prior sequence data output filetype"           )
    (m_option_priorSeqDataOutputAllowAll.c_str(),      boost::program_options::value<bool        >()->default_value(UQ_GCM_PRIOR_SEQ_DATA_OUTPUT_ALLOW_ALL_ODV      ), "allow all or not"                              )
    (m_option_priorSeqDataOutputAllowedSet.c_str(),    boost::program_options::value<std::string >()->default_value(UQ_GCM_PRIOR_SEQ_DATA_OUTPUT_ALLOWED_SET_ODV    ), "subEnvs that will write to data output file"   )
    (m_option_nuggetValueForBtWyB.c_str(),             boost::program_options::value<double      >()->default_value(UQ_GCM_NUGGET_VALUE_FOR_BT_WY_B_ODV             ), "nugget value for Bt_Wy_W matrix"               )
    (m_option_nuggetValueForBtWyBInv.c_str(),          boost::program_options::value<double      >()->default_value(UQ_GCM_NUGGET_VALUE_FOR_BT_WY_B_INV_ODV         ), "nugget value for Bt_Wy_W inverse matrix"       )
    (m_option_formCMatrix.c_str(),                     boost::program_options::value<double      >()->default_value(UQ_GCM_FORM_C_MATRIX_ODV                        ), "form C matrix"                                 )
    (m_option_useTildeLogicForRankDefficientC.c_str(), boost::program_options::value<bool        >()->default_value(UQ_GCM_USE_TILDE_LOGIC_FOR_RANK_DEFFICIENT_C_ODV), "use tilde logic for rank defficient C"         )
    (m_option_predLag.c_str(),                         boost::program_options::value<unsigned int>()->default_value(UQ_GCM_PRED_LAG_ODV                             ), "predLag"                                       )
    (m_option_predVUsBySamplingRVs.c_str(),            boost::program_options::value<bool        >()->default_value(UQ_GCM_PRED_VUS_BY_SAMPLING_RVS_ODV             ), "predVUsBySamplingRVs"                          )
    (m_option_predVUsBySummingRVs.c_str(),             boost::program_options::value<bool        >()->default_value(UQ_GCM_PRED_VUS_BY_SUMMING_RVS_ODV              ), "predVUsBySummingRVs"                           )
    (m_option_predVUsAtKeyPoints.c_str(),              boost::program_options::value<bool        >()->default_value(UQ_GCM_PRED_VUS_AT_KEY_POINTS_ODV               ), "predVUsAtKeyPoints"                            )
    (m_option_predWsBySamplingRVs.c_str(),             boost::program_options::value<bool        >()->default_value(UQ_GCM_PRED_WS_BY_SAMPLING_RVS_ODV              ), "predWsBySamplingRVs"                           )
    (m_option_predWsBySummingRVs.c_str(),              boost::program_options::value<bool        >()->default_value(UQ_GCM_PRED_WS_BY_SUMMING_RVS_ODV               ), "predWsBySummingRVs"                            )
    (m_option_predWsAtKeyPoints.c_str(),               boost::program_options::value<bool        >()->default_value(UQ_GCM_PRED_WS_AT_KEY_POINTS_ODV                ), "predWsAtKeyPoints"                             )
  ;

  return;
}

void
GpmsaComputerModelOptions::getMyOptionValues(boost::program_options::options_description& optionsDesc)
{
  queso_deprecated();

  if (m_env.allOptionsMap().count(m_option_help)) {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << optionsDesc
                              << std::endl;
    }
  }

  if (m_env.allOptionsMap().count(m_option_checkAgainstPreviousSample)) {
    m_ov.m_checkAgainstPreviousSample = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_checkAgainstPreviousSample]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_dataOutputFileName)) {
    m_ov.m_dataOutputFileName = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_dataOutputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_dataOutputAllowAll)) {
    m_ov.m_dataOutputAllowAll = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_dataOutputAllowAll]).as<bool>();
  }

  if (m_ov.m_dataOutputAllowAll) {
    m_ov.m_dataOutputAllowedSet.insert(m_env.subId());
  }
  else if (m_env.allOptionsMap().count(m_option_dataOutputAllowedSet)) {
    m_ov.m_dataOutputAllowedSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_dataOutputAllowedSet].as<std::string>();
    MiscReadDoublesFromString(inputString,tmpAllow);

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_ov.m_dataOutputAllowedSet.insert((unsigned int) tmpAllow[i]);
      }
    }
  }

  //std::cout << "In GpmsaComputerModelOptions::getMyOptionValues(), m_option_priorSeqNumSamples = " << m_option_priorSeqNumSamples << "___" << std::endl;
  if (m_env.allOptionsMap().count(m_option_priorSeqNumSamples)) {
    //std::cout << "In GpmsaComputerModelOptions::getMyOptionValues(), going to read m_option_priorSeqNumSamples..." << std::endl;
    m_ov.m_priorSeqNumSamples = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_priorSeqNumSamples]).as<unsigned int>();
    //std::cout << "In GpmsaComputerModelOptions::getMyOptionValues(), just read m_option_priorSeqNumSamples = " << m_ov.m_priorSeqNumSamples << std::endl;
  }

  if (m_env.allOptionsMap().count(m_option_priorSeqDataOutputFileName)) {
    m_ov.m_priorSeqDataOutputFileName = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_priorSeqDataOutputFileName]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_priorSeqDataOutputFileType)) {
    m_ov.m_priorSeqDataOutputFileType = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_priorSeqDataOutputFileType]).as<std::string>();
  }

  if (m_env.allOptionsMap().count(m_option_priorSeqDataOutputAllowAll)) {
    m_ov.m_priorSeqDataOutputAllowAll = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_priorSeqDataOutputAllowAll]).as<bool>();
  }

  if (m_ov.m_priorSeqDataOutputAllowAll) {
    m_ov.m_priorSeqDataOutputAllowedSet.insert(m_env.subId());
  }
  else if (m_env.allOptionsMap().count(m_option_priorSeqDataOutputAllowedSet)) {
    m_ov.m_priorSeqDataOutputAllowedSet.clear();
    std::vector<double> tmpAllow(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_priorSeqDataOutputAllowedSet].as<std::string>();
    MiscReadDoublesFromString(inputString,tmpAllow);

    if (tmpAllow.size() > 0) {
      for (unsigned int i = 0; i < tmpAllow.size(); ++i) {
        m_ov.m_priorSeqDataOutputAllowedSet.insert((unsigned int) tmpAllow[i]);
      }
    }
  }

  if (m_env.allOptionsMap().count(m_option_nuggetValueForBtWyB)) {
    m_ov.m_nuggetValueForBtWyB = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_nuggetValueForBtWyB]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_nuggetValueForBtWyBInv)) {
    m_ov.m_nuggetValueForBtWyBInv = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_nuggetValueForBtWyBInv]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_formCMatrix)) {
    m_ov.m_formCMatrix = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_formCMatrix]).as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_useTildeLogicForRankDefficientC)) {
    m_ov.m_useTildeLogicForRankDefficientC = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_useTildeLogicForRankDefficientC]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_predLag)) {
    m_ov.m_predLag = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_predLag]).as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_predVUsBySamplingRVs)) {
    m_ov.m_predVUsBySamplingRVs = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_predVUsBySamplingRVs]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_predVUsBySummingRVs)) {
    m_ov.m_predVUsBySummingRVs = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_predVUsBySummingRVs]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_predVUsAtKeyPoints)) {
    m_ov.m_predVUsAtKeyPoints = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_predVUsAtKeyPoints]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_predWsBySamplingRVs)) {
    m_ov.m_predWsBySamplingRVs = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_predWsBySamplingRVs]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_predWsBySummingRVs)) {
    m_ov.m_predWsBySummingRVs = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_predWsBySummingRVs]).as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_predWsAtKeyPoints)) {
    m_ov.m_predWsAtKeyPoints = ((const boost::program_options::variable_value&) m_env.allOptionsMap()[m_option_predWsAtKeyPoints]).as<bool>();
  }

  return;
}

void
GpmsaComputerModelOptions::print(std::ostream& os) const
{
  queso_deprecated();

  os << "\n" << m_option_checkAgainstPreviousSample << " = " << m_ov.m_checkAgainstPreviousSample
     << "\n" << m_option_dataOutputFileName         << " = " << m_ov.m_dataOutputFileName
     << "\n" << m_option_dataOutputAllowAll         << " = " << m_ov.m_dataOutputAllowAll
     << "\n" << m_option_dataOutputAllowedSet       << " = ";
  for (std::set<unsigned int>::iterator setIt = m_ov.m_dataOutputAllowedSet.begin(); setIt != m_ov.m_dataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_priorSeqNumSamples              << " = " << m_ov.m_priorSeqNumSamples
     << "\n" << m_option_priorSeqDataOutputFileName      << " = " << m_ov.m_priorSeqDataOutputFileName
     << "\n" << m_option_priorSeqDataOutputFileType      << " = " << m_ov.m_priorSeqDataOutputFileType
     << "\n" << m_option_priorSeqDataOutputAllowAll      << " = " << m_ov.m_priorSeqDataOutputAllowAll
     << "\n" << m_option_priorSeqDataOutputAllowedSet    << " = ";
  for (std::set<unsigned int>::iterator setIt = m_ov.m_priorSeqDataOutputAllowedSet.begin(); setIt != m_ov.m_priorSeqDataOutputAllowedSet.end(); ++setIt) {
    os << *setIt << " ";
  }
  os << "\n" << m_option_nuggetValueForBtWyB             << " = " << m_ov.m_nuggetValueForBtWyB
     << "\n" << m_option_nuggetValueForBtWyBInv          << " = " << m_ov.m_nuggetValueForBtWyBInv
     << "\n" << m_option_formCMatrix                     << " = " << m_ov.m_formCMatrix
     << "\n" << m_option_useTildeLogicForRankDefficientC << " = " << m_ov.m_useTildeLogicForRankDefficientC
     << "\n" << m_option_predLag                         << " = " << m_ov.m_predLag
     << "\n" << m_option_predVUsBySamplingRVs            << " = " << m_ov.m_predVUsBySamplingRVs
     << "\n" << m_option_predVUsBySummingRVs             << " = " << m_ov.m_predVUsBySummingRVs
     << "\n" << m_option_predVUsAtKeyPoints              << " = " << m_ov.m_predVUsAtKeyPoints
     << "\n" << m_option_predWsBySamplingRVs             << " = " << m_ov.m_predWsBySamplingRVs
     << "\n" << m_option_predWsBySummingRVs              << " = " << m_ov.m_predWsBySummingRVs
     << "\n" << m_option_predWsAtKeyPoints               << " = " << m_ov.m_predWsAtKeyPoints
     << std::endl;

  return;
}

std::ostream& operator<<(std::ostream& os, const GpmsaComputerModelOptions& obj)
{
  queso_deprecated();

  obj.print(os);

  return os;
}

}  // End namespace QUESO
