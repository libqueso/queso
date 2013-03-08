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
// $Id$
//
//--------------------------------------------------------------------------

#include <uqDefines.h>
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS

#include <uqSequenceStatisticalOptions.h>
#include <uqMiscellaneous.h>

uqSsOptionsValuesClass::uqSsOptionsValuesClass()
  :
  m_initialDiscardedPortions(0                                           ),//,0.),
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
  m_meanMonitorPeriod       (UQ_SEQUENCE_MEAN_MONITOR_PERIOD_ODV         ),
  m_bmmRun                  (UQ_SEQUENCE_BMM_RUN_ODV                     ),
  m_bmmLengths              (0                                           ),//,0),
  m_fftCompute              (UQ_SEQUENCE_FFT_COMPUTE_ODV                 ),
  m_fftParamId              (UQ_SEQUENCE_FFT_PARAM_ID_ODV                ),
  m_fftSize                 (UQ_SEQUENCE_FFT_SIZE_ODV                    ),
  m_fftTestInversion        (UQ_SEQUENCE_FFT_TEST_INVERSION_ODV          ),
  m_fftWrite                (UQ_SEQUENCE_FFT_WRITE_ODV                   ),
  m_psdCompute              (UQ_SEQUENCE_PSD_COMPUTE_ODV                 ),
  m_psdNumBlocks            (UQ_SEQUENCE_PSD_NUM_BLOCKS_ODV              ),
  m_psdHopSizeRatio         (UQ_SEQUENCE_PSD_HOP_SIZE_RATIO_ODV          ),
  m_psdParamId              (UQ_SEQUENCE_PSD_PARAM_ID_ODV                ),
  m_psdWrite                (UQ_SEQUENCE_PSD_WRITE_ODV                   ),
  m_psdAtZeroCompute        (UQ_SEQUENCE_PSD_AT_ZERO_COMPUTE_ODV         ),
  m_psdAtZeroNumBlocks      (0                                           ),//,0),
  m_psdAtZeroHopSizeRatio   (UQ_SEQUENCE_PSD_AT_ZERO_HOP_SIZE_RATIO_ODV  ),
  m_psdAtZeroDisplay        (UQ_SEQUENCE_PSD_AT_ZERO_DISPLAY_ODV         ),
  m_psdAtZeroWrite          (UQ_SEQUENCE_PSD_AT_ZERO_WRITE_ODV           ),
  m_gewekeCompute           (UQ_SEQUENCE_GEWEKE_COMPUTE_ODV              ),
  m_gewekeNaRatio           (UQ_SEQUENCE_GEWEKE_NA_RATIO_ODV             ),
  m_gewekeNbRatio           (UQ_SEQUENCE_GEWEKE_NB_RATIO_ODV             ),
  m_gewekeDisplay           (UQ_SEQUENCE_GEWEKE_DISPLAY_ODV              ),
  m_gewekeWrite             (UQ_SEQUENCE_GEWEKE_WRITE_ODV                ),
  m_meanStaccCompute        (UQ_SEQUENCE_MEAN_STACC_COMPUTE_ODV          ),
  m_histCompute             (UQ_SEQUENCE_HIST_COMPUTE_ODV                ),
  m_histNumInternalBins     (UQ_SEQUENCE_HIST_NUM_INTERNAL_BINS_ODV      ),
  m_cdfStaccCompute         (UQ_SEQUENCE_CDF_STACC_COMPUTE_ODV           ),
  m_cdfStaccNumEvalPositions(UQ_SEQUENCE_CDF_STACC_NUM_EVAL_POSITIONS_ODV),
#endif
  m_autoCorrComputeViaDef   (UQ_SEQUENCE_AUTO_CORR_COMPUTE_VIA_DEF_ODV   ),
  m_autoCorrComputeViaFft   (UQ_SEQUENCE_AUTO_CORR_COMPUTE_VIA_FFT_ODV   ),
  m_autoCorrSecondLag       (UQ_SEQUENCE_AUTO_CORR_SECOND_LAG_ODV        ),
  m_autoCorrLagSpacing      (UQ_SEQUENCE_AUTO_CORR_LAG_SPACING_ODV       ),
  m_autoCorrNumLags         (UQ_SEQUENCE_AUTO_CORR_NUM_LAGS_ODV          ),
  m_autoCorrDisplay         (UQ_SEQUENCE_AUTO_CORR_DISPLAY_ODV           ),
  m_autoCorrWrite           (UQ_SEQUENCE_AUTO_CORR_WRITE_ODV             ),
  m_kdeCompute              (UQ_SEQUENCE_KDE_COMPUTE_ODV                 ),
  m_kdeNumEvalPositions     (UQ_SEQUENCE_KDE_NUM_EVAL_POSITIONS_ODV      ),
  m_covMatrixCompute        (UQ_SEQUENCE_COV_MATRIX_COMPUTE_ODV          ),
  m_corrMatrixCompute       (UQ_SEQUENCE_CORR_MATRIX_COMPUTE_ODV         )
{
}

uqSsOptionsValuesClass::~uqSsOptionsValuesClass()
{
}

uqSsOptionsValuesClass::uqSsOptionsValuesClass(const uqSsOptionsValuesClass& src)
{
  this->copy(src);
}

uqSsOptionsValuesClass&
uqSsOptionsValuesClass::operator=(const uqSsOptionsValuesClass& rhs)
{
  this->copy(rhs);
  return *this;
}

void
uqSsOptionsValuesClass::copy(const uqSsOptionsValuesClass& src)
{
  m_initialDiscardedPortions = src.m_initialDiscardedPortions;
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
  m_meanMonitorPeriod        = src.m_meanMonitorPeriod;
  m_bmmRun                   = src.m_bmmRun;
  m_bmmLengths               = src.m_bmmLengths;
  m_fftCompute               = src.m_fftCompute;
  m_fftParamId               = src.m_fftParamId;
  m_fftSize                  = src.m_fftSize;
  m_fftTestInversion         = src.m_fftTestInversion;
  m_fftWrite                 = src.m_fftWrite;
  m_psdCompute               = src.m_psdCompute;
  m_psdNumBlocks             = src.m_psdNumBlocks;
  m_psdHopSizeRatio          = src.m_psdHopSizeRatio;
  m_psdParamId               = src.m_psdParamId;
  m_psdWrite                 = src.m_psdWrite;
  m_psdAtZeroCompute         = src.m_psdAtZeroCompute;
  m_psdAtZeroNumBlocks       = src.m_psdAtZeroNumBlocks;
  m_psdAtZeroHopSizeRatio    = src.m_psdAtZeroHopSizeRatio;
  m_psdAtZeroDisplay         = src.m_psdAtZeroDisplay;
  m_psdAtZeroWrite           = src.m_psdAtZeroWrite;
  m_gewekeCompute            = src.m_gewekeCompute;
  m_gewekeNaRatio            = src.m_gewekeNaRatio;
  m_gewekeNbRatio            = src.m_gewekeNbRatio;
  m_gewekeDisplay            = src.m_gewekeDisplay;
  m_gewekeWrite              = src.m_gewekeWrite;
  m_meanStaccCompute         = src.m_meanStaccCompute;
  m_histCompute              = src.m_histCompute;
  m_histNumInternalBins      = src.m_histNumInternalBins;
  m_cdfStaccCompute          = src.m_cdfStaccCompute;
  m_cdfStaccNumEvalPositions = src.m_cdfStaccNumEvalPositions;
#endif
  m_autoCorrComputeViaDef    = src.m_autoCorrComputeViaDef;
  m_autoCorrComputeViaFft    = src.m_autoCorrComputeViaFft;
  m_autoCorrSecondLag        = src.m_autoCorrSecondLag;
  m_autoCorrLagSpacing       = src.m_autoCorrLagSpacing;
  m_autoCorrNumLags          = src.m_autoCorrNumLags;
  m_autoCorrDisplay          = src.m_autoCorrDisplay;
  m_autoCorrWrite            = src.m_autoCorrWrite;
  m_kdeCompute               = src.m_kdeCompute;
  m_kdeNumEvalPositions      = src.m_kdeNumEvalPositions;
  m_covMatrixCompute         = src.m_covMatrixCompute;
  m_corrMatrixCompute        = src.m_corrMatrixCompute;

  return;
}

uqSequenceStatisticalOptionsClass::uqSequenceStatisticalOptionsClass(
  /*! The QUESO environment                      */ const uqBaseEnvironmentClass& env,
  /*! Prefix for reading options from input file */ const std::string&            prefix)
  :
  m_ov                              (),
  m_prefix                          ((std::string)(prefix) + "stats_"),
  m_env                             (env),
  m_optionsDesc                     (new po::options_description("Chain statistical options")),
  m_option_help                     (m_prefix + "help"                     ),
  m_option_initialDiscardedPortions (m_prefix + "initialDiscardedPortions" ),
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
  m_option_mean_monitorPeriod       (m_prefix + "mean_monitorPeriod"       ),
  m_option_bmm_run                  (m_prefix + "bmm_run"                  ),
  m_option_bmm_lengths              (m_prefix + "bmm_lengths"              ),
  m_option_bmm_display              (m_prefix + "bmm_display"              ),
  m_option_bmm_write                (m_prefix + "bmm_write"                ),
  m_option_fft_compute              (m_prefix + "fft_compute"              ),
  m_option_fft_paramId              (m_prefix + "fft_paramId"              ),
  m_option_fft_size                 (m_prefix + "fft_size"                 ),
  m_option_fft_testInversion        (m_prefix + "fft_testInversion"        ),
  m_option_fft_write                (m_prefix + "fft_write"                ),
  m_option_psd_compute              (m_prefix + "psd_compute"              ),
  m_option_psd_numBlocks            (m_prefix + "psd_numBlocks"            ),
  m_option_psd_hopSizeRatio         (m_prefix + "psd_hopSizeRatio"         ),
  m_option_psd_paramId              (m_prefix + "psd_paramId"              ),
  m_option_psd_write                (m_prefix + "psd_write"                ),
  m_option_psdAtZero_compute        (m_prefix + "psdAtZero_compute"        ),
  m_option_psdAtZero_numBlocks      (m_prefix + "psdAtZero_numBlocks"      ),
  m_option_psdAtZero_hopSizeRatio   (m_prefix + "psdAtZero_hopSizeRatio"   ),
  m_option_psdAtZero_display        (m_prefix + "psdAtZero_display"        ),
  m_option_psdAtZero_write          (m_prefix + "psdAtZero_write"          ),
  m_option_geweke_compute           (m_prefix + "geweke_compute"           ),
  m_option_geweke_naRatio           (m_prefix + "geweke_naRatio"           ),
  m_option_geweke_nbRatio           (m_prefix + "geweke_nbRatio"           ),
  m_option_geweke_display           (m_prefix + "geweke_display"           ),
  m_option_geweke_write             (m_prefix + "geweke_write"             ),
  m_option_meanStacc_compute        (m_prefix + "meanStacc_compute"        ),
  m_option_hist_compute             (m_prefix + "hist_compute"             ),
  m_option_hist_numInternalBins     (m_prefix + "hist_numInternalBins"     ),
  m_option_cdfStacc_compute         (m_prefix + "cdfStacc_compute"         ),
  m_option_cdfStacc_numEvalPositions(m_prefix + "cdfStacc_numEvalPositions"),
#endif
  m_option_autoCorr_computeViaDef   (m_prefix + "autoCorr_computeViaDef"   ),
  m_option_autoCorr_computeViaFft   (m_prefix + "autoCorr_computeViaFft"   ),
  m_option_autoCorr_secondLag       (m_prefix + "autoCorr_secondLag"       ),
  m_option_autoCorr_lagSpacing      (m_prefix + "autoCorr_lagSpacing"      ),
  m_option_autoCorr_numLags         (m_prefix + "autoCorr_numLags"         ),
  m_option_autoCorr_display         (m_prefix + "autoCorr_display"         ),
  m_option_autoCorr_write           (m_prefix + "autoCorr_write"           ),
  m_option_kde_compute              (m_prefix + "kde_compute"              ),
  m_option_kde_numEvalPositions     (m_prefix + "kde_numEvalPositions"     ),
  m_option_covMatrix_compute        (m_prefix + "covMatrix_compute"        ),
  m_option_corrMatrix_compute       (m_prefix + "corrMatrix_compute"       )
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering uqSequenceStatisticalOptions::constructor(1)"
                            << ", prefix = " << m_prefix
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(m_env.optionsInputFileName() == "",
                      m_env.worldRank(),
                      "uqSequenceStatisticalOptionsClass::constructor(1)",
                      "this constructor is incompatible with the abscense of an options input file");

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "After reading values of options with prefix '" << m_prefix
                            << "', state of uqSequenceStatisticalOptionsClass object is:"
                            << "\n" << *this
                            << std::endl;
  }

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving uqSequenceStatisticalOptions::constructor(1)"
                            << ", prefix = " << m_prefix
                            << std::endl;
  }
}

uqSequenceStatisticalOptionsClass::uqSequenceStatisticalOptionsClass(
  /*! The QUESO environment                      */ const uqBaseEnvironmentClass& env,
  /*! Prefix for reading options from input file */ const std::string&            prefix,
                                                    const uqSsOptionsValuesClass& alternativeOptionsValues)
  :
  m_ov                              (alternativeOptionsValues),
  m_prefix                          ((std::string)(prefix) + "stats_"),
  m_env                             (env),
  m_optionsDesc                     (NULL),
  m_option_help                     (m_prefix + "help"                     ),
  m_option_initialDiscardedPortions (m_prefix + "initialDiscardedPortions" ),
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
  m_option_mean_monitorPeriod       (m_prefix + "mean_monitorPeriod"       ),
  m_option_bmm_run                  (m_prefix + "bmm_run"                  ),
  m_option_bmm_lengths              (m_prefix + "bmm_lengths"              ),
  m_option_bmm_display              (m_prefix + "bmm_display"              ),
  m_option_bmm_write                (m_prefix + "bmm_write"                ),
  m_option_fft_compute              (m_prefix + "fft_compute"              ),
  m_option_fft_paramId              (m_prefix + "fft_paramId"              ),
  m_option_fft_size                 (m_prefix + "fft_size"                 ),
  m_option_fft_testInversion        (m_prefix + "fft_testInversion"        ),
  m_option_fft_write                (m_prefix + "fft_write"                ),
  m_option_psd_compute              (m_prefix + "psd_compute"              ),
  m_option_psd_numBlocks            (m_prefix + "psd_numBlocks"            ),
  m_option_psd_hopSizeRatio         (m_prefix + "psd_hopSizeRatio"         ),
  m_option_psd_paramId              (m_prefix + "psd_paramId"              ),
  m_option_psd_write                (m_prefix + "psd_write"                ),
  m_option_psdAtZero_compute        (m_prefix + "psdAtZero_compute"        ),
  m_option_psdAtZero_numBlocks      (m_prefix + "psdAtZero_numBlocks"      ),
  m_option_psdAtZero_hopSizeRatio   (m_prefix + "psdAtZero_hopSizeRatio"   ),
  m_option_psdAtZero_display        (m_prefix + "psdAtZero_display"        ),
  m_option_psdAtZero_write          (m_prefix + "psdAtZero_write"          ),
  m_option_geweke_compute           (m_prefix + "geweke_compute"           ),
  m_option_geweke_naRatio           (m_prefix + "geweke_naRatio"           ),
  m_option_geweke_nbRatio           (m_prefix + "geweke_nbRatio"           ),
  m_option_geweke_display           (m_prefix + "geweke_display"           ),
  m_option_geweke_write             (m_prefix + "geweke_write"             ),
  m_option_meanStacc_compute        (m_prefix + "meanStacc_compute"        ),
  m_option_hist_compute             (m_prefix + "hist_compute"             ),
  m_option_hist_numInternalBins     (m_prefix + "hist_numInternalBins"     ),
  m_option_cdfStacc_compute         (m_prefix + "cdfStacc_compute"         ),
  m_option_cdfStacc_numEvalPositions(m_prefix + "cdfStacc_numEvalPositions"),
#endif
  m_option_autoCorr_computeViaDef   (m_prefix + "autoCorr_computeViaDef"   ),
  m_option_autoCorr_computeViaFft   (m_prefix + "autoCorr_computeViaFft"   ),
  m_option_autoCorr_secondLag       (m_prefix + "autoCorr_secondLag"       ),
  m_option_autoCorr_lagSpacing      (m_prefix + "autoCorr_lagSpacing"      ),
  m_option_autoCorr_numLags         (m_prefix + "autoCorr_numLags"         ),
  m_option_autoCorr_display         (m_prefix + "autoCorr_display"         ),
  m_option_autoCorr_write           (m_prefix + "autoCorr_write"           ),
  m_option_kde_compute              (m_prefix + "kde_compute"              ),
  m_option_kde_numEvalPositions     (m_prefix + "kde_numEvalPositions"     ),
  m_option_covMatrix_compute        (m_prefix + "covMatrix_compute"        ),
  m_option_corrMatrix_compute       (m_prefix + "corrMatrix_compute"       )
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering uqSequenceStatisticalOptions::constructor(2)"
                            << ", prefix = " << m_prefix
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(m_env.optionsInputFileName() != "",
                      m_env.worldRank(),
                      "uqSequenceStatisticalOptionsClass::constructor(2)",
                      "this constructor is incompatible with the existence of an options input file");

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "After setting values of options with prefix '" << m_prefix
                            << "', state of uqSequenceStatisticalOptionsClass object is:"
                            << "\n" << *this
                            << std::endl;
  }

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving uqSequenceStatisticalOptions::constructor(2)"
                            << ", prefix = " << m_prefix
                            << std::endl;
  }
}

uqSequenceStatisticalOptionsClass::~uqSequenceStatisticalOptionsClass()
{
  if (m_optionsDesc) delete m_optionsDesc;
}

void
uqSequenceStatisticalOptionsClass::defineMyOptions(
  po::options_description& optionsDesc) const
{
  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                                                                    "produce help message for chain statistical options"             )
    (m_option_initialDiscardedPortions.c_str(),       po::value<std::string >()->default_value(UQ_SEQUENCE_INITIAL_DISCARDED_PORTIONS_ODV      ), "list of initial discarded portions for chain statistics"        )
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
    (m_option_mean_monitorPeriod.c_str(),             po::value<unsigned int>()->default_value(UQ_SEQUENCE_MEAN_MONITOR_PERIOD_ODV             ), "period for monitoring mean"                                     )
    (m_option_bmm_run.c_str(),                        po::value<bool        >()->default_value(UQ_SEQUENCE_BMM_RUN_ODV                         ), "compute variance of sample mean with batch means method"        )
    (m_option_bmm_lengths.c_str(),                    po::value<std::string >()->default_value(UQ_SEQUENCE_BMM_LENGTHS_ODV                     ), "list of batch lenghts for BMM"                                  )
    (m_option_fft_compute.c_str(),                    po::value<bool        >()->default_value(UQ_SEQUENCE_FFT_COMPUTE_ODV                     ), "compute fft"                                                    )
    (m_option_fft_paramId.c_str(),                    po::value<unsigned int>()->default_value(UQ_SEQUENCE_FFT_PARAM_ID_ODV                    ), "parameter id for fft computations"                              )
    (m_option_fft_size.c_str(),                       po::value<unsigned int>()->default_value(UQ_SEQUENCE_FFT_SIZE_ODV                        ), "fft size"                                                       )
    (m_option_fft_testInversion.c_str(),              po::value<bool        >()->default_value(UQ_SEQUENCE_FFT_TEST_INVERSION_ODV              ), "test fft inversion"                                             )
    (m_option_fft_write.c_str(),                      po::value<bool        >()->default_value(UQ_SEQUENCE_FFT_WRITE_ODV                       ), "write fft"                                                      )
    (m_option_psd_compute.c_str(),                    po::value<bool        >()->default_value(UQ_SEQUENCE_PSD_COMPUTE_ODV                     ), "compute psd"                                                    )
    (m_option_psd_numBlocks.c_str(),                  po::value<unsigned int>()->default_value(UQ_SEQUENCE_PSD_NUM_BLOCKS_ODV                  ), "number of blocks for computation of psd"                        )
    (m_option_psd_hopSizeRatio.c_str(),               po::value<double      >()->default_value(UQ_SEQUENCE_PSD_HOP_SIZE_RATIO_ODV              ), "hop size ratio for psd"                                         )
    (m_option_psd_paramId.c_str(),                    po::value<unsigned int>()->default_value(UQ_SEQUENCE_PSD_PARAM_ID_ODV                    ), "parameter id for psd computations"                              )
    (m_option_psd_write.c_str(),                      po::value<bool        >()->default_value(UQ_SEQUENCE_PSD_WRITE_ODV                       ), "write psd"                                                      )
    (m_option_psdAtZero_compute.c_str(),              po::value<bool        >()->default_value(UQ_SEQUENCE_PSD_AT_ZERO_COMPUTE_ODV             ), "compute power spectral densities"                               )
    (m_option_psdAtZero_numBlocks.c_str(),            po::value<std::string >()->default_value(UQ_SEQUENCE_PSD_AT_ZERO_NUM_BLOCKS_ODV          ), "list of numbers of blocks for computation of psd at zero"       )
    (m_option_psdAtZero_hopSizeRatio.c_str(),         po::value<double      >()->default_value(UQ_SEQUENCE_PSD_AT_ZERO_HOP_SIZE_RATIO_ODV      ), "hop size ratio for psd at zero"                                 )
    (m_option_psdAtZero_display.c_str(),              po::value<bool        >()->default_value(UQ_SEQUENCE_PSD_AT_ZERO_DISPLAY_ODV             ), "display computed psd at frequency zero on screen"               )
    (m_option_psdAtZero_write.c_str(),                po::value<bool        >()->default_value(UQ_SEQUENCE_PSD_AT_ZERO_WRITE_ODV               ), "write computed psd at frequency zero to the output file"        )
    (m_option_geweke_compute.c_str(),                 po::value<bool        >()->default_value(UQ_SEQUENCE_GEWEKE_COMPUTE_ODV                  ), "compute Geweke coefficients"                                    )
    (m_option_geweke_naRatio.c_str(),                 po::value<double      >()->default_value(UQ_SEQUENCE_GEWEKE_NA_RATIO_ODV                 ), "ratio NA for Geweke"                                            ) 
    (m_option_geweke_nbRatio.c_str(),                 po::value<double      >()->default_value(UQ_SEQUENCE_GEWEKE_NB_RATIO_ODV                 ), "ratio NB for Geweke"                                            )
    (m_option_geweke_display.c_str(),                 po::value<bool        >()->default_value(UQ_SEQUENCE_GEWEKE_DISPLAY_ODV                  ), "display computed Geweke on screen"                              )
    (m_option_geweke_write.c_str(),                   po::value<bool        >()->default_value(UQ_SEQUENCE_GEWEKE_WRITE_ODV                    ), "write computed Geweke to the output file"                       )
    (m_option_meanStacc_compute.c_str(),              po::value<bool        >()->default_value(UQ_SEQUENCE_MEAN_STACC_COMPUTE_ODV              ), "compute statistical accuracy of mean"                           )
    (m_option_hist_compute.c_str(),                   po::value<bool        >()->default_value(UQ_SEQUENCE_HIST_COMPUTE_ODV                    ), "compute histograms"                                             )
    (m_option_hist_numInternalBins.c_str(),           po::value<unsigned int>()->default_value(UQ_SEQUENCE_HIST_NUM_INTERNAL_BINS_ODV          ), "number of internal bins"                                        )
    (m_option_cdfStacc_compute.c_str(),               po::value<bool        >()->default_value(UQ_SEQUENCE_CDF_STACC_COMPUTE_ODV               ), "compute statisical accuracy of cdf"                             )
    (m_option_cdfStacc_numEvalPositions.c_str(),      po::value<unsigned int>()->default_value(UQ_SEQUENCE_CDF_STACC_NUM_EVAL_POSITIONS_ODV    ), "number of evaluations points for statistical accuracy of cdf"   )
#endif
    (m_option_autoCorr_computeViaDef.c_str(),         po::value<bool        >()->default_value(UQ_SEQUENCE_AUTO_CORR_COMPUTE_VIA_DEF_ODV       ), "compute correlations via definition"                            )
    (m_option_autoCorr_computeViaFft.c_str(),         po::value<bool        >()->default_value(UQ_SEQUENCE_AUTO_CORR_COMPUTE_VIA_FFT_ODV       ), "compute correlations via fft"                                   )
    (m_option_autoCorr_secondLag.c_str(),             po::value<unsigned int>()->default_value(UQ_SEQUENCE_AUTO_CORR_SECOND_LAG_ODV            ), "second lag for computation of autocorrelations"                 )
    (m_option_autoCorr_lagSpacing.c_str(),            po::value<unsigned int>()->default_value(UQ_SEQUENCE_AUTO_CORR_LAG_SPACING_ODV           ), "lag spacing for computation of autocorrelations"                )
    (m_option_autoCorr_numLags.c_str(),               po::value<unsigned int>()->default_value(UQ_SEQUENCE_AUTO_CORR_NUM_LAGS_ODV              ), "number of lags for computation of autocorrelations"             )
    (m_option_autoCorr_display.c_str(),               po::value<bool        >()->default_value(UQ_SEQUENCE_AUTO_CORR_DISPLAY_ODV               ), "display computed autocorrelations on the screen"                )
    (m_option_autoCorr_write.c_str(),                 po::value<bool        >()->default_value(UQ_SEQUENCE_AUTO_CORR_WRITE_ODV                 ), "write computed autocorrelations to the output file"             )
    (m_option_kde_compute.c_str(),                    po::value<bool        >()->default_value(UQ_SEQUENCE_KDE_COMPUTE_ODV                     ), "compute kernel density estimators"                              )
    (m_option_kde_numEvalPositions.c_str(),           po::value<unsigned int>()->default_value(UQ_SEQUENCE_KDE_NUM_EVAL_POSITIONS_ODV          ), "number of evaluation positions"                                 )
    (m_option_covMatrix_compute.c_str(),              po::value<bool        >()->default_value(UQ_SEQUENCE_COV_MATRIX_COMPUTE_ODV              ), "compute covariance matrix"                                      )
    (m_option_corrMatrix_compute.c_str(),             po::value<bool        >()->default_value(UQ_SEQUENCE_CORR_MATRIX_COMPUTE_ODV             ), "compute correlation matrix"                                     )
  ;

  return;
}

void
uqSequenceStatisticalOptionsClass::getMyOptionValues(
  po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help)) {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << optionsDesc
                             << std::endl;
    }
  }

  if (m_env.allOptionsMap().count(m_option_initialDiscardedPortions)) {
    m_ov.m_initialDiscardedPortions.clear();
    std::vector<double> tmpPortions(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_initialDiscardedPortions].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpPortions);
    //if (m_env.subDisplayFile()) {
    //  *m_env.subDisplayFile() << "In uqSequenceStatisticalOptionsClass::getMyOptionValues(): percents = ";
    //  for (unsigned int i = 0; i < tmpPortions.size(); ++i) {
    //    *m_env.subDisplayFile() << " " << tmpPortions[i];
    //  }
    //  *m_env.subDisplayFile() << std::endl;
    //}

    if (tmpPortions.size() > 0) {
      m_ov.m_initialDiscardedPortions.resize(tmpPortions.size(),0.);
      for (unsigned int i = 0; i < m_ov.m_initialDiscardedPortions.size(); ++i) {
        m_ov.m_initialDiscardedPortions[i] = tmpPortions[i];
      }
    }
  }

#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
  if (m_env.allOptionsMap().count(m_option_mean_monitorPeriod)) {
    m_ov.m_meanMonitorPeriod = m_env.allOptionsMap()[m_option_mean_monitorPeriod].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_bmm_run)) {
    m_ov.m_bmmRun = m_env.allOptionsMap()[m_option_bmm_run].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_bmm_lengths)) {
    m_ov.m_bmmLengths.clear();
    std::vector<double> tmpLengths(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_bmm_lengths].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpLengths);
    //if (m_env.subDisplayFile()) {
    //  *m_env.subDisplayFile() << "In uqSequenceStatisticalOptionsClass::getMyOptionValues(): lengths for BMM = ";
    //  for (unsigned int i = 0; i < tmpLengths.size(); ++i) {
    //    *m_env.subDisplayFile() << " " << tmpLengths[i];
    //  }
    //  *m_env.subDisplayFile() << std::endl;
    //}

    if (tmpLengths.size() > 0) {
      m_ov.m_bmmLengths.resize(tmpLengths.size(),0);
      for (unsigned int i = 0; i < m_ov.m_bmmLengths.size(); ++i) {
        m_ov.m_bmmLengths[i] = (unsigned int) tmpLengths[i];
      }
    }
  }

  if (m_env.allOptionsMap().count(m_option_fft_compute)) {
    m_ov.m_fftCompute = m_env.allOptionsMap()[m_option_fft_compute].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_fft_paramId)) {
    m_ov.m_fftParamId = m_env.allOptionsMap()[m_option_fft_paramId].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_fft_size)) {
    m_ov.m_fftSize = m_env.allOptionsMap()[m_option_fft_size].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_fft_testInversion)) {
    m_ov.m_fftTestInversion = m_env.allOptionsMap()[m_option_fft_testInversion].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_fft_write)) {
    m_ov.m_fftWrite = m_env.allOptionsMap()[m_option_fft_write].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_psd_compute)) {
    m_ov.m_psdCompute = m_env.allOptionsMap()[m_option_psd_compute].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_psd_numBlocks)) {
    m_ov.m_psdNumBlocks = m_env.allOptionsMap()[m_option_psd_numBlocks].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_psd_hopSizeRatio)) {
    m_ov.m_psdHopSizeRatio = m_env.allOptionsMap()[m_option_psd_hopSizeRatio].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_psd_paramId)) {
    m_ov.m_psdParamId = m_env.allOptionsMap()[m_option_psd_paramId].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_psd_write)) {
    m_ov.m_psdWrite = m_env.allOptionsMap()[m_option_psd_write].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_psdAtZero_compute)) {
    m_ov.m_psdAtZeroCompute = m_env.allOptionsMap()[m_option_psdAtZero_compute].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_psdAtZero_numBlocks)) {
    m_ov.m_psdAtZeroNumBlocks.clear();
    std::vector<double> tmpNumBlocks(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_psdAtZero_numBlocks].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpNumBlocks);
    //if (m_env.subDisplayFile()) {
    //  *m_env.subDisplayFile() << "In uqSequenceStatisticalOptionsClass::getMyOptionValues(): numBlocks for psdAtZero = ";
    //  for (unsigned int i = 0; i < tmpNumBlocks.size(); ++i) {
    //    *m_env.subDisplayFile() << " " << numBlocks[i];
    //  }
    //  *m_env.subDisplayFile() << std::endl;
    //}

    if (tmpNumBlocks.size() > 0) {
      m_ov.m_psdAtZeroNumBlocks.resize(tmpNumBlocks.size(),0);
      for (unsigned int i = 0; i < m_ov.m_psdAtZeroNumBlocks.size(); ++i) {
        m_ov.m_psdAtZeroNumBlocks[i] = (unsigned int) tmpNumBlocks[i];
      }
    }
  }

  if (m_env.allOptionsMap().count(m_option_psdAtZero_hopSizeRatio)) {
    m_ov.m_psdAtZeroHopSizeRatio = m_env.allOptionsMap()[m_option_psdAtZero_hopSizeRatio].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_psdAtZero_display)) {
    m_ov.m_psdAtZeroDisplay = m_env.allOptionsMap()[m_option_psdAtZero_display].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_psdAtZero_write)) {
    m_ov.m_psdAtZeroWrite = m_env.allOptionsMap()[m_option_psdAtZero_write].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_geweke_compute)) {
    m_ov.m_gewekeCompute = m_env.allOptionsMap()[m_option_geweke_compute].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_geweke_naRatio)) {
    m_ov.m_gewekeNaRatio = m_env.allOptionsMap()[m_option_geweke_naRatio].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_geweke_nbRatio)) {
    m_ov.m_gewekeNbRatio = m_env.allOptionsMap()[m_option_geweke_nbRatio].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_geweke_display)) {
    m_ov.m_gewekeDisplay = m_env.allOptionsMap()[m_option_geweke_display].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_geweke_write)) {
    m_ov.m_gewekeWrite = m_env.allOptionsMap()[m_option_geweke_write].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_meanStacc_compute)) {
    m_ov.m_meanStaccCompute = m_env.allOptionsMap()[m_option_meanStacc_compute].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_hist_compute)) {
    m_ov.m_histCompute = m_env.allOptionsMap()[m_option_hist_compute].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_hist_numInternalBins)) {
    m_ov.m_histNumInternalBins = m_env.allOptionsMap()[m_option_hist_numInternalBins].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_cdfStacc_compute)) {
    m_ov.m_cdfStaccCompute = m_env.allOptionsMap()[m_option_cdfStacc_compute].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_cdfStacc_numEvalPositions)) {
    m_ov.m_cdfStaccNumEvalPositions = m_env.allOptionsMap()[m_option_cdfStacc_numEvalPositions].as<unsigned int>();
  }
#endif
  if (m_env.allOptionsMap().count(m_option_autoCorr_computeViaDef)) {
    m_ov.m_autoCorrComputeViaDef = m_env.allOptionsMap()[m_option_autoCorr_computeViaDef].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_autoCorr_computeViaFft)) {
    m_ov.m_autoCorrComputeViaFft = m_env.allOptionsMap()[m_option_autoCorr_computeViaFft].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_autoCorr_secondLag)) {
    m_ov.m_autoCorrSecondLag = m_env.allOptionsMap()[m_option_autoCorr_secondLag].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_autoCorr_lagSpacing)) {
    m_ov.m_autoCorrLagSpacing = m_env.allOptionsMap()[m_option_autoCorr_lagSpacing].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_autoCorr_numLags)) {
    m_ov.m_autoCorrNumLags = m_env.allOptionsMap()[m_option_autoCorr_numLags].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_autoCorr_display)) {
    m_ov.m_autoCorrDisplay = m_env.allOptionsMap()[m_option_autoCorr_display].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_autoCorr_write)) {
    m_ov.m_autoCorrWrite = m_env.allOptionsMap()[m_option_autoCorr_write].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_kde_compute)) {
    m_ov.m_kdeCompute = m_env.allOptionsMap()[m_option_kde_compute].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_kde_numEvalPositions)) {
    m_ov.m_kdeNumEvalPositions = m_env.allOptionsMap()[m_option_kde_numEvalPositions].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_covMatrix_compute)) {
    m_ov.m_covMatrixCompute = m_env.allOptionsMap()[m_option_covMatrix_compute].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_corrMatrix_compute)) {
    m_ov.m_corrMatrixCompute = m_env.allOptionsMap()[m_option_corrMatrix_compute].as<bool>();
  }

  return;
}

const std::vector<double>&
uqSequenceStatisticalOptionsClass::initialDiscardedPortions() const
{
  return m_ov.m_initialDiscardedPortions;
}

#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
unsigned int
uqSequenceStatisticalOptionsClass::meanMonitorPeriod() const
{
  return m_ov.m_meanMonitorPeriod;
}

bool
uqSequenceStatisticalOptionsClass::bmmRun() const
{
  return m_ov.m_bmmRun;
}

const std::vector<unsigned int>&
uqSequenceStatisticalOptionsClass::bmmLengths() const
{
  return m_ov.m_bmmLengths;
}

bool
uqSequenceStatisticalOptionsClass::bmmDisplay() const
{
  return m_ov.m_bmmDisplay;
}

bool
uqSequenceStatisticalOptionsClass::bmmWrite() const
{
  return m_ov.m_bmmWrite;
}

bool
uqSequenceStatisticalOptionsClass::fftCompute() const
{
  return m_ov.m_fftCompute;
}

unsigned int
uqSequenceStatisticalOptionsClass::fftParamId() const
{
  return m_ov.m_fftParamId;
}

unsigned int
uqSequenceStatisticalOptionsClass::fftSize() const
{
  return m_ov.m_fftSize;
}

bool
uqSequenceStatisticalOptionsClass::fftTestInversion() const
{
  return m_ov.m_fftTestInversion;
}

bool
uqSequenceStatisticalOptionsClass::fftWrite() const
{
  return m_ov.m_fftWrite;
}

bool
uqSequenceStatisticalOptionsClass::psdCompute() const
{
  return m_ov.m_psdCompute;
}

unsigned int
uqSequenceStatisticalOptionsClass::psdNumBlocks() const
{
  return m_ov.m_psdNumBlocks;
}

double
uqSequenceStatisticalOptionsClass::psdHopSizeRatio() const
{
  return m_ov.m_psdHopSizeRatio;
}

unsigned int
uqSequenceStatisticalOptionsClass::psdParamId() const
{
  return m_ov.m_psdParamId;
}

bool
uqSequenceStatisticalOptionsClass::psdWrite() const
{
  return m_ov.m_psdWrite;
}

bool
uqSequenceStatisticalOptionsClass::psdAtZeroCompute() const
{
  return m_ov.m_psdAtZeroCompute;
}

const std::vector<unsigned int>&
uqSequenceStatisticalOptionsClass::psdAtZeroNumBlocks() const
{
  return m_ov.m_psdAtZeroNumBlocks;
}

double
uqSequenceStatisticalOptionsClass::psdAtZeroHopSizeRatio() const
{
  return m_ov.m_psdAtZeroHopSizeRatio;
}

bool
uqSequenceStatisticalOptionsClass::psdAtZeroDisplay() const
{
  return m_ov.m_psdAtZeroDisplay;
}

bool
uqSequenceStatisticalOptionsClass::psdAtZeroWrite() const
{
  return m_ov.m_psdAtZeroWrite;
}

bool
uqSequenceStatisticalOptionsClass::gewekeCompute() const
{
  return m_ov.m_gewekeCompute;
}

double
uqSequenceStatisticalOptionsClass::gewekeNaRatio() const
{
  return m_ov.m_gewekeNaRatio;
}

double
uqSequenceStatisticalOptionsClass::gewekeNbRatio() const
{
  return m_ov.m_gewekeNbRatio;
}

bool
uqSequenceStatisticalOptionsClass::gewekeDisplay() const
{
  return m_ov.m_gewekeDisplay;
}

bool
uqSequenceStatisticalOptionsClass::gewekeWrite() const
{
  return m_ov.m_gewekeWrite;
}

bool
uqSequenceStatisticalOptionsClass::meanStaccCompute() const
{
  return m_ov.m_meanStaccCompute;
}

bool
uqSequenceStatisticalOptionsClass::histCompute() const
{
  return m_ov.m_histCompute;
}

unsigned int
uqSequenceStatisticalOptionsClass::histNumInternalBins() const
{
  return m_ov.m_histNumInternalBins;
}

bool
uqSequenceStatisticalOptionsClass::cdfStaccCompute() const
{
  return m_ov.m_cdfStaccCompute;
}

unsigned int
uqSequenceStatisticalOptionsClass::cdfStaccNumEvalPositions() const
{
  return m_ov.m_cdfStaccNumEvalPositions;
}
#endif
bool
uqSequenceStatisticalOptionsClass::autoCorrComputeViaDef() const
{
  return m_ov.m_autoCorrComputeViaDef;
}

bool
uqSequenceStatisticalOptionsClass::autoCorrComputeViaFft() const
{
  return m_ov.m_autoCorrComputeViaFft;
}

unsigned int
uqSequenceStatisticalOptionsClass::autoCorrSecondLag() const
{
  return m_ov.m_autoCorrSecondLag;
}

unsigned int
uqSequenceStatisticalOptionsClass::autoCorrLagSpacing() const
{
  return m_ov.m_autoCorrLagSpacing;
}

unsigned int
uqSequenceStatisticalOptionsClass::autoCorrNumLags() const
{
  return m_ov.m_autoCorrNumLags;
}

bool
uqSequenceStatisticalOptionsClass::autoCorrDisplay() const
{
  return m_ov.m_autoCorrDisplay;
}

bool
uqSequenceStatisticalOptionsClass::autoCorrWrite() const
{
  return m_ov.m_autoCorrWrite;
}

bool
uqSequenceStatisticalOptionsClass::kdeCompute() const
{
  return m_ov.m_kdeCompute;
}

unsigned int
uqSequenceStatisticalOptionsClass::kdeNumEvalPositions() const
{
  return m_ov.m_kdeNumEvalPositions;
}

bool
uqSequenceStatisticalOptionsClass::covMatrixCompute() const
{
  return m_ov.m_covMatrixCompute;
}

bool
uqSequenceStatisticalOptionsClass::corrMatrixCompute() const
{
  return m_ov.m_corrMatrixCompute;
}

void
uqSequenceStatisticalOptionsClass::print(std::ostream& os) const
{
  os << "\n" << m_option_initialDiscardedPortions << " = ";
  for (unsigned int i = 0; i < m_ov.m_initialDiscardedPortions.size(); ++i) {
    os << m_ov.m_initialDiscardedPortions[i] << " ";
  }
  os
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
     << "\n" << m_option_mean_monitorPeriod << " = " << m_ov.m_meanMonitorPeriod
     << "\n" << m_option_bmm_run            << " = " << m_ov.m_bmmRun
     << "\n" << m_option_bmm_lengths        << " = ";
  for (unsigned int i = 0; i < m_ov.m_bmmLengths.size(); ++i) {
    os << m_ov.m_bmmLengths[i] << " ";
  }
  os << "\n" << m_option_fft_compute         << " = " << m_ov.m_fftCompute
     << "\n" << m_option_fft_paramId         << " = " << m_ov.m_fftParamId
     << "\n" << m_option_fft_size            << " = " << m_ov.m_fftSize
     << "\n" << m_option_fft_testInversion   << " = " << m_ov.m_fftTestInversion
     << "\n" << m_option_fft_write           << " = " << m_ov.m_fftWrite
     << "\n" << m_option_psd_compute         << " = " << m_ov.m_psdCompute
     << "\n" << m_option_psd_paramId         << " = " << m_ov.m_psdParamId
     << "\n" << m_option_psd_numBlocks       << " = " << m_ov.m_psdNumBlocks
     << "\n" << m_option_psd_hopSizeRatio    << " = " << m_ov.m_psdHopSizeRatio
     << "\n" << m_option_psd_write           << " = " << m_ov.m_psdWrite
     << "\n" << m_option_psdAtZero_compute   << " = " << m_ov.m_psdAtZeroCompute
     << "\n" << m_option_psdAtZero_numBlocks << " = ";
  for (unsigned int i = 0; i < m_ov.m_psdAtZeroNumBlocks.size(); ++i) {
    os << m_ov.m_psdAtZeroNumBlocks[i] << " ";
  }
  os << "\n" << m_option_psdAtZero_hopSizeRatio    << " = " << m_ov.m_psdAtZeroHopSizeRatio
     << "\n" << m_option_psdAtZero_display         << " = " << m_ov.m_psdAtZeroDisplay
     << "\n" << m_option_psdAtZero_write           << " = " << m_ov.m_psdAtZeroWrite
     << "\n" << m_option_geweke_compute            << " = " << m_ov.m_gewekeCompute
     << "\n" << m_option_geweke_naRatio            << " = " << m_ov.m_gewekeNaRatio
     << "\n" << m_option_geweke_nbRatio            << " = " << m_ov.m_gewekeNbRatio
     << "\n" << m_option_geweke_display            << " = " << m_ov.m_gewekeDisplay
     << "\n" << m_option_geweke_write              << " = " << m_ov.m_gewekeWrite
     << "\n" << m_option_meanStacc_compute         << " = " << m_ov.m_meanStaccCompute
     << "\n" << m_option_hist_compute              << " = " << m_ov.m_histCompute
     << "\n" << m_option_hist_numInternalBins      << " = " << m_ov.m_histNumInternalBins
     << "\n" << m_option_cdfStacc_compute          << " = " << m_ov.m_cdfStaccCompute
     << "\n" << m_option_cdfStacc_numEvalPositions << " = " << m_ov.m_cdfStaccNumEvalPositions
#endif
     << "\n" << m_option_autoCorr_computeViaDef    << " = " << m_ov.m_autoCorrComputeViaDef
     << "\n" << m_option_autoCorr_computeViaFft    << " = " << m_ov.m_autoCorrComputeViaFft
     << "\n" << m_option_autoCorr_secondLag        << " = " << m_ov.m_autoCorrSecondLag
     << "\n" << m_option_autoCorr_lagSpacing       << " = " << m_ov.m_autoCorrLagSpacing
     << "\n" << m_option_autoCorr_numLags          << " = " << m_ov.m_autoCorrNumLags
     << "\n" << m_option_autoCorr_display          << " = " << m_ov.m_autoCorrDisplay
     << "\n" << m_option_autoCorr_write            << " = " << m_ov.m_autoCorrWrite
     << "\n" << m_option_kde_compute               << " = " << m_ov.m_kdeCompute
     << "\n" << m_option_kde_numEvalPositions      << " = " << m_ov.m_kdeNumEvalPositions
     << "\n" << m_option_covMatrix_compute         << " = " << m_ov.m_covMatrixCompute
     << "\n" << m_option_corrMatrix_compute        << " = " << m_ov.m_corrMatrixCompute
     << std::endl;

  return;
}

std::ostream&
operator<<(std::ostream& os, const uqSequenceStatisticalOptionsClass& obj)
{
  obj.print(os);

  return os;
}

#endif // ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
