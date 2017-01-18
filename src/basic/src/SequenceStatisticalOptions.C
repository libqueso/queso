//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2017 The PECOS Development Team
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

#include <queso/Defines.h>
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

#include <queso/SequenceStatisticalOptions.h>
#include <queso/Miscellaneous.h>

namespace QUESO {

SsOptionsValues::SsOptionsValues()
  :
    m_prefix                          ("stats_"),
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
    m_corrMatrixCompute       (UQ_SEQUENCE_CORR_MATRIX_COMPUTE_ODV         ),
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
}

SsOptionsValues::SsOptionsValues(const BaseEnvironment * env, const char *
    prefix)
  :
    m_prefix                          ((std::string)(prefix) + "stats_"),
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
    m_corrMatrixCompute       (UQ_SEQUENCE_CORR_MATRIX_COMPUTE_ODV         ),
    m_parser(new BoostInputOptionsParser(env->optionsInputFileName())),
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
  (m_option_help,                                                                                                                    "produce help message for chain statistical options"             )
  (m_option_initialDiscardedPortions,       boost::program_options::value<std::string >()->default_value(UQ_SEQUENCE_INITIAL_DISCARDED_PORTIONS_ODV      ), "list of initial discarded portions for chain statistics"        )
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
  m_parser->registerOption(m_option_mean_monitorPeriod,             UQ_SEQUENCE_MEAN_MONITOR_PERIOD_ODV, "period for monitoring mean"                                     );
  m_parser->registerOption(m_option_bmm_run,                        UQ_SEQUENCE_BMM_RUN_ODV            , "compute variance of sample mean with batch means method"        );
  m_parser->registerOption<std::string>(m_option_bmm_lengths,                    UQ_SEQUENCE_BMM_LENGTHS_ODV        , "list of batch lenghts for BMM"                                  );
  m_parser->registerOption(m_option_fft_compute,                    UQ_SEQUENCE_FFT_COMPUTE_ODV        , "compute fft"                                                    );
  m_parser->registerOption(m_option_fft_paramId,                    UQ_SEQUENCE_FFT_PARAM_ID_ODV       , "parameter id for fft computations"                              );
  m_parser->registerOption(m_option_fft_size,                       UQ_SEQUENCE_FFT_SIZE_ODV           , "fft size"                                                       );
  m_parser->registerOption(m_option_fft_testInversion,              UQ_SEQUENCE_FFT_TEST_INVERSION_ODV , "test fft inversion"                                             );
  m_parser->registerOption(m_option_fft_write,                      UQ_SEQUENCE_FFT_WRITE_ODV          , "write fft"                                                      );
  m_parser->registerOption(m_option_psd_compute,                    UQ_SEQUENCE_PSD_COMPUTE_ODV        , "compute psd"                                                    );
  m_parser->registerOption(m_option_psd_numBlocks,                  UQ_SEQUENCE_PSD_NUM_BLOCKS_ODV     , "number of blocks for computation of psd"                        );
  m_parser->registerOption(m_option_psd_hopSizeRatio,               UQ_SEQUENCE_PSD_HOP_SIZE_RATIO_ODV , "hop size ratio for psd"                                         );
  m_parser->registerOption(m_option_psd_paramId,                    UQ_SEQUENCE_PSD_PARAM_ID_ODV       , "parameter id for psd computations"                              );
  m_parser->registerOption(m_option_psd_write,                      UQ_SEQUENCE_PSD_WRITE_ODV          , "write psd"                                                      );
  m_parser->registerOption(m_option_psdAtZero_compute,              UQ_SEQUENCE_PSD_AT_ZERO_COMPUTE_ODV, "compute power spectral densities"                               );
  m_parser->registerOption<std::string>(m_option_psdAtZero_numBlocks,            UQ_SEQUENCE_PSD_AT_ZERO_NUM_BLOCKS_ODV          , "list of numbers of blocks for computation of psd at zero"       );
  m_parser->registerOption(m_option_psdAtZero_hopSizeRatio,         UQ_SEQUENCE_PSD_AT_ZERO_HOP_SIZE_RATIO_ODV      , "hop size ratio for psd at zero"                                 );
  m_parser->registerOption(m_option_psdAtZero_display,              UQ_SEQUENCE_PSD_AT_ZERO_DISPLAY_ODV, "display computed psd at frequency zero on screen"               );
  m_parser->registerOption(m_option_psdAtZero_write,                UQ_SEQUENCE_PSD_AT_ZERO_WRITE_ODV  , "write computed psd at frequency zero to the output file"        );
  m_parser->registerOption(m_option_geweke_compute,                 UQ_SEQUENCE_GEWEKE_COMPUTE_ODV     , "compute Geweke coefficients"                                    );
  m_parser->registerOption(m_option_geweke_naRatio,                 UQ_SEQUENCE_GEWEKE_NA_RATIO_ODV    , "ratio NA for Geweke"                                            );
  m_parser->registerOption(m_option_geweke_nbRatio,                 UQ_SEQUENCE_GEWEKE_NB_RATIO_ODV    , "ratio NB for Geweke"                                            );
  m_parser->registerOption(m_option_geweke_display,                 UQ_SEQUENCE_GEWEKE_DISPLAY_ODV     , "display computed Geweke on screen"                              );
  m_parser->registerOption(m_option_geweke_write,                   UQ_SEQUENCE_GEWEKE_WRITE_ODV       , "write computed Geweke to the output file"                       );
  m_parser->registerOption(m_option_meanStacc_compute,              UQ_SEQUENCE_MEAN_STACC_COMPUTE_ODV , "compute statistical accuracy of mean"                           );
  m_parser->registerOption(m_option_hist_compute,                   UQ_SEQUENCE_HIST_COMPUTE_ODV       , "compute histograms"                                             );
  m_parser->registerOption(m_option_hist_numInternalBins,           UQ_SEQUENCE_HIST_NUM_INTERNAL_BINS_ODV          , "number of internal bins"                                        );
  m_parser->registerOption(m_option_cdfStacc_compute,               UQ_SEQUENCE_CDF_STACC_COMPUTE_ODV  , "compute statisical accuracy of cdf"                             );
  m_parser->registerOption(m_option_cdfStacc_numEvalPositions,      UQ_SEQUENCE_CDF_STACC_NUM_EVAL_POSITIONS_ODV    , "number of evaluations points for statistical accuracy of cdf"   );
#endif
  m_parser->registerOption(m_option_autoCorr_computeViaDef,         boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_AUTO_CORR_COMPUTE_VIA_DEF_ODV       ), "compute correlations via definition"                            )
  m_parser->registerOption(m_option_autoCorr_computeViaFft,         boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_AUTO_CORR_COMPUTE_VIA_FFT_ODV       ), "compute correlations via fft"                                   )
  m_parser->registerOption(m_option_autoCorr_secondLag,             boost::program_options::value<unsigned int>()->default_value(UQ_SEQUENCE_AUTO_CORR_SECOND_LAG_ODV            ), "second lag for computation of autocorrelations"                 )
  m_parser->registerOption(m_option_autoCorr_lagSpacing,            boost::program_options::value<unsigned int>()->default_value(UQ_SEQUENCE_AUTO_CORR_LAG_SPACING_ODV           ), "lag spacing for computation of autocorrelations"                )
  m_parser->registerOption(m_option_autoCorr_numLags,               boost::program_options::value<unsigned int>()->default_value(UQ_SEQUENCE_AUTO_CORR_NUM_LAGS_ODV              ), "number of lags for computation of autocorrelations"             )
  m_parser->registerOption(m_option_autoCorr_display,               boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_AUTO_CORR_DISPLAY_ODV               ), "display computed autocorrelations on the screen"                )
  m_parser->registerOption(m_option_autoCorr_write,                 boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_AUTO_CORR_WRITE_ODV                 ), "write computed autocorrelations to the output file"             )
  m_parser->registerOption(m_option_kde_compute,                    boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_KDE_COMPUTE_ODV                     ), "compute kernel density estimators"                              )
  m_parser->registerOption(m_option_kde_numEvalPositions,           boost::program_options::value<unsigned int>()->default_value(UQ_SEQUENCE_KDE_NUM_EVAL_POSITIONS_ODV          ), "number of evaluation positions"                                 )
  m_parser->registerOption(m_option_covMatrix_compute,              boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_COV_MATRIX_COMPUTE_ODV              ), "compute covariance matrix"                                      )
  m_parser->registerOption(m_option_corrMatrix_compute,             boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_CORR_MATRIX_COMPUTE_ODV             ), "compute correlation matrix"                                     )
}

SsOptionsValues::~SsOptionsValues()
{
}

SsOptionsValues::SsOptionsValues(const SsOptionsValues& src)
{
  this->copy(src);
}

SsOptionsValues&
SsOptionsValues::operator=(const SsOptionsValues& rhs)
{
  this->copy(rhs);
  return *this;
}

void
SsOptionsValues::defineOptions()
{
  (*m_optionsDescription).add_options()
    (m_option_help.c_str(),                                                                                                                    "produce help message for chain statistical options"             )
    (m_option_initialDiscardedPortions.c_str(),       boost::program_options::value<std::string >()->default_value(UQ_SEQUENCE_INITIAL_DISCARDED_PORTIONS_ODV      ), "list of initial discarded portions for chain statistics"        )
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
    (m_option_mean_monitorPeriod.c_str(),             boost::program_options::value<unsigned int>()->default_value(UQ_SEQUENCE_MEAN_MONITOR_PERIOD_ODV             ), "period for monitoring mean"                                     )
    (m_option_bmm_run.c_str(),                        boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_BMM_RUN_ODV                         ), "compute variance of sample mean with batch means method"        )
    (m_option_bmm_lengths.c_str(),                    boost::program_options::value<std::string >()->default_value(UQ_SEQUENCE_BMM_LENGTHS_ODV                     ), "list of batch lenghts for BMM"                                  )
    (m_option_fft_compute.c_str(),                    boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_FFT_COMPUTE_ODV                     ), "compute fft"                                                    )
    (m_option_fft_paramId.c_str(),                    boost::program_options::value<unsigned int>()->default_value(UQ_SEQUENCE_FFT_PARAM_ID_ODV                    ), "parameter id for fft computations"                              )
    (m_option_fft_size.c_str(),                       boost::program_options::value<unsigned int>()->default_value(UQ_SEQUENCE_FFT_SIZE_ODV                        ), "fft size"                                                       )
    (m_option_fft_testInversion.c_str(),              boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_FFT_TEST_INVERSION_ODV              ), "test fft inversion"                                             )
    (m_option_fft_write.c_str(),                      boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_FFT_WRITE_ODV                       ), "write fft"                                                      )
    (m_option_psd_compute.c_str(),                    boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_PSD_COMPUTE_ODV                     ), "compute psd"                                                    )
    (m_option_psd_numBlocks.c_str(),                  boost::program_options::value<unsigned int>()->default_value(UQ_SEQUENCE_PSD_NUM_BLOCKS_ODV                  ), "number of blocks for computation of psd"                        )
    (m_option_psd_hopSizeRatio.c_str(),               boost::program_options::value<double      >()->default_value(UQ_SEQUENCE_PSD_HOP_SIZE_RATIO_ODV              ), "hop size ratio for psd"                                         )
    (m_option_psd_paramId.c_str(),                    boost::program_options::value<unsigned int>()->default_value(UQ_SEQUENCE_PSD_PARAM_ID_ODV                    ), "parameter id for psd computations"                              )
    (m_option_psd_write.c_str(),                      boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_PSD_WRITE_ODV                       ), "write psd"                                                      )
    (m_option_psdAtZero_compute.c_str(),              boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_PSD_AT_ZERO_COMPUTE_ODV             ), "compute power spectral densities"                               )
    (m_option_psdAtZero_numBlocks.c_str(),            boost::program_options::value<std::string >()->default_value(UQ_SEQUENCE_PSD_AT_ZERO_NUM_BLOCKS_ODV          ), "list of numbers of blocks for computation of psd at zero"       )
    (m_option_psdAtZero_hopSizeRatio.c_str(),         boost::program_options::value<double      >()->default_value(UQ_SEQUENCE_PSD_AT_ZERO_HOP_SIZE_RATIO_ODV      ), "hop size ratio for psd at zero"                                 )
    (m_option_psdAtZero_display.c_str(),              boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_PSD_AT_ZERO_DISPLAY_ODV             ), "display computed psd at frequency zero on screen"               )
    (m_option_psdAtZero_write.c_str(),                boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_PSD_AT_ZERO_WRITE_ODV               ), "write computed psd at frequency zero to the output file"        )
    (m_option_geweke_compute.c_str(),                 boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_GEWEKE_COMPUTE_ODV                  ), "compute Geweke coefficients"                                    )
    (m_option_geweke_naRatio.c_str(),                 boost::program_options::value<double      >()->default_value(UQ_SEQUENCE_GEWEKE_NA_RATIO_ODV                 ), "ratio NA for Geweke"                                            )
    (m_option_geweke_nbRatio.c_str(),                 boost::program_options::value<double      >()->default_value(UQ_SEQUENCE_GEWEKE_NB_RATIO_ODV                 ), "ratio NB for Geweke"                                            )
    (m_option_geweke_display.c_str(),                 boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_GEWEKE_DISPLAY_ODV                  ), "display computed Geweke on screen"                              )
    (m_option_geweke_write.c_str(),                   boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_GEWEKE_WRITE_ODV                    ), "write computed Geweke to the output file"                       )
    (m_option_meanStacc_compute.c_str(),              boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_MEAN_STACC_COMPUTE_ODV              ), "compute statistical accuracy of mean"                           )
    (m_option_hist_compute.c_str(),                   boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_HIST_COMPUTE_ODV                    ), "compute histograms"                                             )
    (m_option_hist_numInternalBins.c_str(),           boost::program_options::value<unsigned int>()->default_value(UQ_SEQUENCE_HIST_NUM_INTERNAL_BINS_ODV          ), "number of internal bins"                                        )
    (m_option_cdfStacc_compute.c_str(),               boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_CDF_STACC_COMPUTE_ODV               ), "compute statisical accuracy of cdf"                             )
    (m_option_cdfStacc_numEvalPositions.c_str(),      boost::program_options::value<unsigned int>()->default_value(UQ_SEQUENCE_CDF_STACC_NUM_EVAL_POSITIONS_ODV    ), "number of evaluations points for statistical accuracy of cdf"   )
#endif
    (m_option_autoCorr_computeViaDef.c_str(),         boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_AUTO_CORR_COMPUTE_VIA_DEF_ODV       ), "compute correlations via definition"                            )
    (m_option_autoCorr_computeViaFft.c_str(),         boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_AUTO_CORR_COMPUTE_VIA_FFT_ODV       ), "compute correlations via fft"                                   )
    (m_option_autoCorr_secondLag.c_str(),             boost::program_options::value<unsigned int>()->default_value(UQ_SEQUENCE_AUTO_CORR_SECOND_LAG_ODV            ), "second lag for computation of autocorrelations"                 )
    (m_option_autoCorr_lagSpacing.c_str(),            boost::program_options::value<unsigned int>()->default_value(UQ_SEQUENCE_AUTO_CORR_LAG_SPACING_ODV           ), "lag spacing for computation of autocorrelations"                )
    (m_option_autoCorr_numLags.c_str(),               boost::program_options::value<unsigned int>()->default_value(UQ_SEQUENCE_AUTO_CORR_NUM_LAGS_ODV              ), "number of lags for computation of autocorrelations"             )
    (m_option_autoCorr_display.c_str(),               boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_AUTO_CORR_DISPLAY_ODV               ), "display computed autocorrelations on the screen"                )
    (m_option_autoCorr_write.c_str(),                 boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_AUTO_CORR_WRITE_ODV                 ), "write computed autocorrelations to the output file"             )
    (m_option_kde_compute.c_str(),                    boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_KDE_COMPUTE_ODV                     ), "compute kernel density estimators"                              )
    (m_option_kde_numEvalPositions.c_str(),           boost::program_options::value<unsigned int>()->default_value(UQ_SEQUENCE_KDE_NUM_EVAL_POSITIONS_ODV          ), "number of evaluation positions"                                 )
    (m_option_covMatrix_compute.c_str(),              boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_COV_MATRIX_COMPUTE_ODV              ), "compute covariance matrix"                                      )
    (m_option_corrMatrix_compute.c_str(),             boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_CORR_MATRIX_COMPUTE_ODV             ), "compute correlation matrix"                                     )
  ;
}

void
SsOptionsValues::getOptionValues()
{
  if ((*m_optionsMap).count(m_option_help)) {
    if (m_env->subDisplayFile()) {
      *m_env->subDisplayFile() << (*m_optionsDescription)
                             << std::endl;
    }
  }

  if ((*m_optionsMap).count(m_option_initialDiscardedPortions)) {
    m_initialDiscardedPortions.clear();
    std::vector<double> tmpPortions(0,0.);
    std::string inputString = (*m_optionsMap)[m_option_initialDiscardedPortions].as<std::string>();
    MiscReadDoublesFromString(inputString,tmpPortions);
    //if (m_env->subDisplayFile()) {
    //  *m_env->subDisplayFile() << "In SequenceStatisticalOptions::getMyOptionValues(): percents = ";
    //  for (unsigned int i = 0; i < tmpPortions.size(); ++i) {
    //    *m_env->subDisplayFile() << " " << tmpPortions[i];
    //  }
    //  *m_env->subDisplayFile() << std::endl;
    //}

    if (tmpPortions.size() > 0) {
      m_initialDiscardedPortions.resize(tmpPortions.size(),0.);
      for (unsigned int i = 0; i < m_initialDiscardedPortions.size(); ++i) {
        m_initialDiscardedPortions[i] = tmpPortions[i];
      }
    }
  }

#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
  if ((*m_optionsMap).count(m_option_mean_monitorPeriod)) {
    m_meanMonitorPeriod = (*m_optionsMap)[m_option_mean_monitorPeriod].as<unsigned int>();
  }

  if ((*m_optionsMap).count(m_option_bmm_run)) {
    m_bmmRun = (*m_optionsMap)[m_option_bmm_run].as<bool>();
  }

  if ((*m_optionsMap).count(m_option_bmm_lengths)) {
    m_bmmLengths.clear();
    std::vector<double> tmpLengths(0,0.);
    std::string inputString = (*m_optionsMap)[m_option_bmm_lengths].as<std::string>();
    MiscReadDoublesFromString(inputString,tmpLengths);
    //if (m_env->subDisplayFile()) {
    //  *m_env->subDisplayFile() << "In SequenceStatisticalOptions::getMyOptionValues(): lengths for BMM = ";
    //  for (unsigned int i = 0; i < tmpLengths.size(); ++i) {
    //    *m_env->subDisplayFile() << " " << tmpLengths[i];
    //  }
    //  *m_env->subDisplayFile() << std::endl;
    //}

    if (tmpLengths.size() > 0) {
      m_bmmLengths.resize(tmpLengths.size(),0);
      for (unsigned int i = 0; i < m_bmmLengths.size(); ++i) {
        m_bmmLengths[i] = (unsigned int) tmpLengths[i];
      }
    }
  }

  if ((*m_optionsMap).count(m_option_fft_compute)) {
    m_fftCompute = (*m_optionsMap)[m_option_fft_compute].as<bool>();
  }

  if ((*m_optionsMap).count(m_option_fft_paramId)) {
    m_fftParamId = (*m_optionsMap)[m_option_fft_paramId].as<unsigned int>();
  }

  if ((*m_optionsMap).count(m_option_fft_size)) {
    m_fftSize = (*m_optionsMap)[m_option_fft_size].as<unsigned int>();
  }

  if ((*m_optionsMap).count(m_option_fft_testInversion)) {
    m_fftTestInversion = (*m_optionsMap)[m_option_fft_testInversion].as<bool>();
  }

  if ((*m_optionsMap).count(m_option_fft_write)) {
    m_fftWrite = (*m_optionsMap)[m_option_fft_write].as<bool>();
  }

  if ((*m_optionsMap).count(m_option_psd_compute)) {
    m_psdCompute = (*m_optionsMap)[m_option_psd_compute].as<bool>();
  }

  if ((*m_optionsMap).count(m_option_psd_numBlocks)) {
    m_psdNumBlocks = (*m_optionsMap)[m_option_psd_numBlocks].as<unsigned int>();
  }

  if ((*m_optionsMap).count(m_option_psd_hopSizeRatio)) {
    m_psdHopSizeRatio = (*m_optionsMap)[m_option_psd_hopSizeRatio].as<double>();
  }

  if ((*m_optionsMap).count(m_option_psd_paramId)) {
    m_psdParamId = (*m_optionsMap)[m_option_psd_paramId].as<unsigned int>();
  }

  if ((*m_optionsMap).count(m_option_psd_write)) {
    m_psdWrite = (*m_optionsMap)[m_option_psd_write].as<bool>();
  }

  if ((*m_optionsMap).count(m_option_psdAtZero_compute)) {
    m_psdAtZeroCompute = (*m_optionsMap)[m_option_psdAtZero_compute].as<bool>();
  }

  if ((*m_optionsMap).count(m_option_psdAtZero_numBlocks)) {
    m_psdAtZeroNumBlocks.clear();
    std::vector<double> tmpNumBlocks(0,0.);
    std::string inputString = (*m_optionsMap)[m_option_psdAtZero_numBlocks].as<std::string>();
    MiscReadDoublesFromString(inputString,tmpNumBlocks);
    //if (m_env->subDisplayFile()) {
    //  *m_env->subDisplayFile() << "In SequenceStatisticalOptions::getMyOptionValues(): numBlocks for psdAtZero = ";
    //  for (unsigned int i = 0; i < tmpNumBlocks.size(); ++i) {
    //    *m_env->subDisplayFile() << " " << numBlocks[i];
    //  }
    //  *m_env->subDisplayFile() << std::endl;
    //}

    if (tmpNumBlocks.size() > 0) {
      m_psdAtZeroNumBlocks.resize(tmpNumBlocks.size(),0);
      for (unsigned int i = 0; i < m_psdAtZeroNumBlocks.size(); ++i) {
        m_psdAtZeroNumBlocks[i] = (unsigned int) tmpNumBlocks[i];
      }
    }
  }

  if ((*m_optionsMap).count(m_option_psdAtZero_hopSizeRatio)) {
    m_psdAtZeroHopSizeRatio = (*m_optionsMap)[m_option_psdAtZero_hopSizeRatio].as<double>();
  }

  if ((*m_optionsMap).count(m_option_psdAtZero_display)) {
    m_psdAtZeroDisplay = (*m_optionsMap)[m_option_psdAtZero_display].as<bool>();
  }

  if ((*m_optionsMap).count(m_option_psdAtZero_write)) {
    m_psdAtZeroWrite = (*m_optionsMap)[m_option_psdAtZero_write].as<bool>();
  }

  if ((*m_optionsMap).count(m_option_geweke_compute)) {
    m_gewekeCompute = (*m_optionsMap)[m_option_geweke_compute].as<bool>();
  }

  if ((*m_optionsMap).count(m_option_geweke_naRatio)) {
    m_gewekeNaRatio = (*m_optionsMap)[m_option_geweke_naRatio].as<double>();
  }

  if ((*m_optionsMap).count(m_option_geweke_nbRatio)) {
    m_gewekeNbRatio = (*m_optionsMap)[m_option_geweke_nbRatio].as<double>();
  }

  if ((*m_optionsMap).count(m_option_geweke_display)) {
    m_gewekeDisplay = (*m_optionsMap)[m_option_geweke_display].as<bool>();
  }

  if ((*m_optionsMap).count(m_option_geweke_write)) {
    m_gewekeWrite = (*m_optionsMap)[m_option_geweke_write].as<bool>();
  }

  if ((*m_optionsMap).count(m_option_meanStacc_compute)) {
    m_meanStaccCompute = (*m_optionsMap)[m_option_meanStacc_compute].as<bool>();
  }

  if ((*m_optionsMap).count(m_option_hist_compute)) {
    m_histCompute = (*m_optionsMap)[m_option_hist_compute].as<bool>();
  }

  if ((*m_optionsMap).count(m_option_hist_numInternalBins)) {
    m_histNumInternalBins = (*m_optionsMap)[m_option_hist_numInternalBins].as<unsigned int>();
  }

  if ((*m_optionsMap).count(m_option_cdfStacc_compute)) {
    m_cdfStaccCompute = (*m_optionsMap)[m_option_cdfStacc_compute].as<bool>();
  }

  if ((*m_optionsMap).count(m_option_cdfStacc_numEvalPositions)) {
    m_cdfStaccNumEvalPositions = (*m_optionsMap)[m_option_cdfStacc_numEvalPositions].as<unsigned int>();
  }
#endif
  if ((*m_optionsMap).count(m_option_autoCorr_computeViaDef)) {
    m_autoCorrComputeViaDef = (*m_optionsMap)[m_option_autoCorr_computeViaDef].as<bool>();
  }

  if ((*m_optionsMap).count(m_option_autoCorr_computeViaFft)) {
    m_autoCorrComputeViaFft = (*m_optionsMap)[m_option_autoCorr_computeViaFft].as<bool>();
  }

  if ((*m_optionsMap).count(m_option_autoCorr_secondLag)) {
    m_autoCorrSecondLag = (*m_optionsMap)[m_option_autoCorr_secondLag].as<unsigned int>();
  }

  if ((*m_optionsMap).count(m_option_autoCorr_lagSpacing)) {
    m_autoCorrLagSpacing = (*m_optionsMap)[m_option_autoCorr_lagSpacing].as<unsigned int>();
  }

  if ((*m_optionsMap).count(m_option_autoCorr_numLags)) {
    m_autoCorrNumLags = (*m_optionsMap)[m_option_autoCorr_numLags].as<unsigned int>();
  }

  if ((*m_optionsMap).count(m_option_autoCorr_display)) {
    m_autoCorrDisplay = (*m_optionsMap)[m_option_autoCorr_display].as<bool>();
  }

  if ((*m_optionsMap).count(m_option_autoCorr_write)) {
    m_autoCorrWrite = (*m_optionsMap)[m_option_autoCorr_write].as<bool>();
  }

  if ((*m_optionsMap).count(m_option_kde_compute)) {
    m_kdeCompute = (*m_optionsMap)[m_option_kde_compute].as<bool>();
  }

  if ((*m_optionsMap).count(m_option_kde_numEvalPositions)) {
    m_kdeNumEvalPositions = (*m_optionsMap)[m_option_kde_numEvalPositions].as<unsigned int>();
  }

  if ((*m_optionsMap).count(m_option_covMatrix_compute)) {
    m_covMatrixCompute = (*m_optionsMap)[m_option_covMatrix_compute].as<bool>();
  }

  if ((*m_optionsMap).count(m_option_corrMatrix_compute)) {
    m_corrMatrixCompute = (*m_optionsMap)[m_option_corrMatrix_compute].as<bool>();
  }
}

void
SsOptionsValues::copy(const SsOptionsValues& src)
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

SequenceStatisticalOptions::SequenceStatisticalOptions(
  /*! The QUESO environment                      */ const BaseEnvironment& env,
  /*! Prefix for reading options from input file */ const std::string&            prefix)
  :
  m_ov                              (),
  m_prefix                          ((std::string)(prefix) + "stats_"),
  m_env                             (env),
  m_optionsDesc                     (new boost::program_options::options_description("Chain statistical options")),
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
  queso_deprecated();

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering SequenceStatisticalOptions::constructor(1)"
                            << ", prefix = " << m_prefix
                            << std::endl;
  }

  queso_require_not_equal_to_msg(m_env.optionsInputFileName(), "", "this constructor is incompatible with the abscense of an options input file");

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "After reading values of options with prefix '" << m_prefix
                            << "', state of SequenceStatisticalOptions object is:"
                            << "\n" << *this
                            << std::endl;
  }

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving SequenceStatisticalOptions::constructor(1)"
                            << ", prefix = " << m_prefix
                            << std::endl;
  }
}

SequenceStatisticalOptions::SequenceStatisticalOptions(
  /*! The QUESO environment                      */ const BaseEnvironment& env,
  /*! Prefix for reading options from input file */ const std::string&            prefix,
                                                    const SsOptionsValues& alternativeOptionsValues)
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
  queso_deprecated();

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering SequenceStatisticalOptions::constructor(2)"
                            << ", prefix = " << m_prefix
                            << std::endl;
  }

  queso_require_equal_to_msg(m_env.optionsInputFileName(), "", "this constructor is incompatible with the existence of an options input file");

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "After setting values of options with prefix '" << m_prefix
                            << "', state of SequenceStatisticalOptions object is:"
                            << "\n" << *this
                            << std::endl;
  }

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving SequenceStatisticalOptions::constructor(2)"
                            << ", prefix = " << m_prefix
                            << std::endl;
  }
}

SequenceStatisticalOptions::~SequenceStatisticalOptions()
{
  queso_deprecated();

  if (m_optionsDesc) delete m_optionsDesc;
}

void
SequenceStatisticalOptions::defineMyOptions(
  boost::program_options::options_description& optionsDesc) const
{
  queso_deprecated();

  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                                                                    "produce help message for chain statistical options"             )
    (m_option_initialDiscardedPortions.c_str(),       boost::program_options::value<std::string >()->default_value(UQ_SEQUENCE_INITIAL_DISCARDED_PORTIONS_ODV      ), "list of initial discarded portions for chain statistics"        )
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
    (m_option_mean_monitorPeriod.c_str(),             boost::program_options::value<unsigned int>()->default_value(UQ_SEQUENCE_MEAN_MONITOR_PERIOD_ODV             ), "period for monitoring mean"                                     )
    (m_option_bmm_run.c_str(),                        boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_BMM_RUN_ODV                         ), "compute variance of sample mean with batch means method"        )
    (m_option_bmm_lengths.c_str(),                    boost::program_options::value<std::string >()->default_value(UQ_SEQUENCE_BMM_LENGTHS_ODV                     ), "list of batch lenghts for BMM"                                  )
    (m_option_fft_compute.c_str(),                    boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_FFT_COMPUTE_ODV                     ), "compute fft"                                                    )
    (m_option_fft_paramId.c_str(),                    boost::program_options::value<unsigned int>()->default_value(UQ_SEQUENCE_FFT_PARAM_ID_ODV                    ), "parameter id for fft computations"                              )
    (m_option_fft_size.c_str(),                       boost::program_options::value<unsigned int>()->default_value(UQ_SEQUENCE_FFT_SIZE_ODV                        ), "fft size"                                                       )
    (m_option_fft_testInversion.c_str(),              boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_FFT_TEST_INVERSION_ODV              ), "test fft inversion"                                             )
    (m_option_fft_write.c_str(),                      boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_FFT_WRITE_ODV                       ), "write fft"                                                      )
    (m_option_psd_compute.c_str(),                    boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_PSD_COMPUTE_ODV                     ), "compute psd"                                                    )
    (m_option_psd_numBlocks.c_str(),                  boost::program_options::value<unsigned int>()->default_value(UQ_SEQUENCE_PSD_NUM_BLOCKS_ODV                  ), "number of blocks for computation of psd"                        )
    (m_option_psd_hopSizeRatio.c_str(),               boost::program_options::value<double      >()->default_value(UQ_SEQUENCE_PSD_HOP_SIZE_RATIO_ODV              ), "hop size ratio for psd"                                         )
    (m_option_psd_paramId.c_str(),                    boost::program_options::value<unsigned int>()->default_value(UQ_SEQUENCE_PSD_PARAM_ID_ODV                    ), "parameter id for psd computations"                              )
    (m_option_psd_write.c_str(),                      boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_PSD_WRITE_ODV                       ), "write psd"                                                      )
    (m_option_psdAtZero_compute.c_str(),              boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_PSD_AT_ZERO_COMPUTE_ODV             ), "compute power spectral densities"                               )
    (m_option_psdAtZero_numBlocks.c_str(),            boost::program_options::value<std::string >()->default_value(UQ_SEQUENCE_PSD_AT_ZERO_NUM_BLOCKS_ODV          ), "list of numbers of blocks for computation of psd at zero"       )
    (m_option_psdAtZero_hopSizeRatio.c_str(),         boost::program_options::value<double      >()->default_value(UQ_SEQUENCE_PSD_AT_ZERO_HOP_SIZE_RATIO_ODV      ), "hop size ratio for psd at zero"                                 )
    (m_option_psdAtZero_display.c_str(),              boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_PSD_AT_ZERO_DISPLAY_ODV             ), "display computed psd at frequency zero on screen"               )
    (m_option_psdAtZero_write.c_str(),                boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_PSD_AT_ZERO_WRITE_ODV               ), "write computed psd at frequency zero to the output file"        )
    (m_option_geweke_compute.c_str(),                 boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_GEWEKE_COMPUTE_ODV                  ), "compute Geweke coefficients"                                    )
    (m_option_geweke_naRatio.c_str(),                 boost::program_options::value<double      >()->default_value(UQ_SEQUENCE_GEWEKE_NA_RATIO_ODV                 ), "ratio NA for Geweke"                                            )
    (m_option_geweke_nbRatio.c_str(),                 boost::program_options::value<double      >()->default_value(UQ_SEQUENCE_GEWEKE_NB_RATIO_ODV                 ), "ratio NB for Geweke"                                            )
    (m_option_geweke_display.c_str(),                 boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_GEWEKE_DISPLAY_ODV                  ), "display computed Geweke on screen"                              )
    (m_option_geweke_write.c_str(),                   boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_GEWEKE_WRITE_ODV                    ), "write computed Geweke to the output file"                       )
    (m_option_meanStacc_compute.c_str(),              boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_MEAN_STACC_COMPUTE_ODV              ), "compute statistical accuracy of mean"                           )
    (m_option_hist_compute.c_str(),                   boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_HIST_COMPUTE_ODV                    ), "compute histograms"                                             )
    (m_option_hist_numInternalBins.c_str(),           boost::program_options::value<unsigned int>()->default_value(UQ_SEQUENCE_HIST_NUM_INTERNAL_BINS_ODV          ), "number of internal bins"                                        )
    (m_option_cdfStacc_compute.c_str(),               boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_CDF_STACC_COMPUTE_ODV               ), "compute statisical accuracy of cdf"                             )
    (m_option_cdfStacc_numEvalPositions.c_str(),      boost::program_options::value<unsigned int>()->default_value(UQ_SEQUENCE_CDF_STACC_NUM_EVAL_POSITIONS_ODV    ), "number of evaluations points for statistical accuracy of cdf"   )
#endif
    (m_option_autoCorr_computeViaDef.c_str(),         boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_AUTO_CORR_COMPUTE_VIA_DEF_ODV       ), "compute correlations via definition"                            )
    (m_option_autoCorr_computeViaFft.c_str(),         boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_AUTO_CORR_COMPUTE_VIA_FFT_ODV       ), "compute correlations via fft"                                   )
    (m_option_autoCorr_secondLag.c_str(),             boost::program_options::value<unsigned int>()->default_value(UQ_SEQUENCE_AUTO_CORR_SECOND_LAG_ODV            ), "second lag for computation of autocorrelations"                 )
    (m_option_autoCorr_lagSpacing.c_str(),            boost::program_options::value<unsigned int>()->default_value(UQ_SEQUENCE_AUTO_CORR_LAG_SPACING_ODV           ), "lag spacing for computation of autocorrelations"                )
    (m_option_autoCorr_numLags.c_str(),               boost::program_options::value<unsigned int>()->default_value(UQ_SEQUENCE_AUTO_CORR_NUM_LAGS_ODV              ), "number of lags for computation of autocorrelations"             )
    (m_option_autoCorr_display.c_str(),               boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_AUTO_CORR_DISPLAY_ODV               ), "display computed autocorrelations on the screen"                )
    (m_option_autoCorr_write.c_str(),                 boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_AUTO_CORR_WRITE_ODV                 ), "write computed autocorrelations to the output file"             )
    (m_option_kde_compute.c_str(),                    boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_KDE_COMPUTE_ODV                     ), "compute kernel density estimators"                              )
    (m_option_kde_numEvalPositions.c_str(),           boost::program_options::value<unsigned int>()->default_value(UQ_SEQUENCE_KDE_NUM_EVAL_POSITIONS_ODV          ), "number of evaluation positions"                                 )
    (m_option_covMatrix_compute.c_str(),              boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_COV_MATRIX_COMPUTE_ODV              ), "compute covariance matrix"                                      )
    (m_option_corrMatrix_compute.c_str(),             boost::program_options::value<bool        >()->default_value(UQ_SEQUENCE_CORR_MATRIX_COMPUTE_ODV             ), "compute correlation matrix"                                     )
  ;

  return;
}

void
SequenceStatisticalOptions::getMyOptionValues(
  boost::program_options::options_description& optionsDesc)
{
  queso_deprecated();

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
    MiscReadDoublesFromString(inputString,tmpPortions);
    //if (m_env.subDisplayFile()) {
    //  *m_env.subDisplayFile() << "In SequenceStatisticalOptions::getMyOptionValues(): percents = ";
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
    MiscReadDoublesFromString(inputString,tmpLengths);
    //if (m_env.subDisplayFile()) {
    //  *m_env.subDisplayFile() << "In SequenceStatisticalOptions::getMyOptionValues(): lengths for BMM = ";
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
    MiscReadDoublesFromString(inputString,tmpNumBlocks);
    //if (m_env.subDisplayFile()) {
    //  *m_env.subDisplayFile() << "In SequenceStatisticalOptions::getMyOptionValues(): numBlocks for psdAtZero = ";
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
SequenceStatisticalOptions::initialDiscardedPortions() const
{
  queso_deprecated();

  return m_ov.m_initialDiscardedPortions;
}

#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
unsigned int
SequenceStatisticalOptions::meanMonitorPeriod() const
{
  queso_deprecated();

  return m_ov.m_meanMonitorPeriod;
}

bool
SequenceStatisticalOptions::bmmRun() const
{
  queso_deprecated();

  return m_ov.m_bmmRun;
}

const std::vector<unsigned int>&
SequenceStatisticalOptions::bmmLengths() const
{
  queso_deprecated();

  return m_ov.m_bmmLengths;
}

bool
SequenceStatisticalOptions::bmmDisplay() const
{
  queso_deprecated();

  return m_ov.m_bmmDisplay;
}

bool
SequenceStatisticalOptions::bmmWrite() const
{
  queso_deprecated();

  return m_ov.m_bmmWrite;
}

bool
SequenceStatisticalOptions::fftCompute() const
{
  queso_deprecated();

  return m_ov.m_fftCompute;
}

unsigned int
SequenceStatisticalOptions::fftParamId() const
{
  queso_deprecated();

  return m_ov.m_fftParamId;
}

unsigned int
SequenceStatisticalOptions::fftSize() const
{
  queso_deprecated();

  return m_ov.m_fftSize;
}

bool
SequenceStatisticalOptions::fftTestInversion() const
{
  queso_deprecated();

  return m_ov.m_fftTestInversion;
}

bool
SequenceStatisticalOptions::fftWrite() const
{
  queso_deprecated();

  return m_ov.m_fftWrite;
}

bool
SequenceStatisticalOptions::psdCompute() const
{
  queso_deprecated();

  return m_ov.m_psdCompute;
}

unsigned int
SequenceStatisticalOptions::psdNumBlocks() const
{
  queso_deprecated();

  return m_ov.m_psdNumBlocks;
}

double
SequenceStatisticalOptions::psdHopSizeRatio() const
{
  queso_deprecated();

  return m_ov.m_psdHopSizeRatio;
}

unsigned int
SequenceStatisticalOptions::psdParamId() const
{
  queso_deprecated();

  return m_ov.m_psdParamId;
}

bool
SequenceStatisticalOptions::psdWrite() const
{
  queso_deprecated();

  return m_ov.m_psdWrite;
}

bool
SequenceStatisticalOptions::psdAtZeroCompute() const
{
  queso_deprecated();

  return m_ov.m_psdAtZeroCompute;
}

const std::vector<unsigned int>&
SequenceStatisticalOptions::psdAtZeroNumBlocks() const
{
  queso_deprecated();

  return m_ov.m_psdAtZeroNumBlocks;
}

double
SequenceStatisticalOptions::psdAtZeroHopSizeRatio() const
{
  queso_deprecated();

  return m_ov.m_psdAtZeroHopSizeRatio;
}

bool
SequenceStatisticalOptions::psdAtZeroDisplay() const
{
  queso_deprecated();

  return m_ov.m_psdAtZeroDisplay;
}

bool
SequenceStatisticalOptions::psdAtZeroWrite() const
{
  queso_deprecated();

  return m_ov.m_psdAtZeroWrite;
}

bool
SequenceStatisticalOptions::gewekeCompute() const
{
  queso_deprecated();

  return m_ov.m_gewekeCompute;
}

double
SequenceStatisticalOptions::gewekeNaRatio() const
{
  queso_deprecated();

  return m_ov.m_gewekeNaRatio;
}

double
SequenceStatisticalOptions::gewekeNbRatio() const
{
  queso_deprecated();

  return m_ov.m_gewekeNbRatio;
}

bool
SequenceStatisticalOptions::gewekeDisplay() const
{
  queso_deprecated();

  return m_ov.m_gewekeDisplay;
}

bool
SequenceStatisticalOptions::gewekeWrite() const
{
  queso_deprecated();

  return m_ov.m_gewekeWrite;
}

bool
SequenceStatisticalOptions::meanStaccCompute() const
{
  queso_deprecated();

  return m_ov.m_meanStaccCompute;
}

bool
SequenceStatisticalOptions::histCompute() const
{
  queso_deprecated();

  return m_ov.m_histCompute;
}

unsigned int
SequenceStatisticalOptions::histNumInternalBins() const
{
  queso_deprecated();

  return m_ov.m_histNumInternalBins;
}

bool
SequenceStatisticalOptions::cdfStaccCompute() const
{
  queso_deprecated();

  return m_ov.m_cdfStaccCompute;
}

unsigned int
SequenceStatisticalOptions::cdfStaccNumEvalPositions() const
{
  queso_deprecated();

  return m_ov.m_cdfStaccNumEvalPositions;
}
#endif
bool
SequenceStatisticalOptions::autoCorrComputeViaDef() const
{
  queso_deprecated();

  return m_ov.m_autoCorrComputeViaDef;
}

bool
SequenceStatisticalOptions::autoCorrComputeViaFft() const
{
  queso_deprecated();

  return m_ov.m_autoCorrComputeViaFft;
}

unsigned int
SequenceStatisticalOptions::autoCorrSecondLag() const
{
  queso_deprecated();

  return m_ov.m_autoCorrSecondLag;
}

unsigned int
SequenceStatisticalOptions::autoCorrLagSpacing() const
{
  queso_deprecated();

  return m_ov.m_autoCorrLagSpacing;
}

unsigned int
SequenceStatisticalOptions::autoCorrNumLags() const
{
  queso_deprecated();

  return m_ov.m_autoCorrNumLags;
}

bool
SequenceStatisticalOptions::autoCorrDisplay() const
{
  queso_deprecated();

  return m_ov.m_autoCorrDisplay;
}

bool
SequenceStatisticalOptions::autoCorrWrite() const
{
  queso_deprecated();

  return m_ov.m_autoCorrWrite;
}

bool
SequenceStatisticalOptions::kdeCompute() const
{
  queso_deprecated();

  return m_ov.m_kdeCompute;
}

unsigned int
SequenceStatisticalOptions::kdeNumEvalPositions() const
{
  queso_deprecated();

  return m_ov.m_kdeNumEvalPositions;
}

bool
SequenceStatisticalOptions::covMatrixCompute() const
{
  queso_deprecated();

  return m_ov.m_covMatrixCompute;
}

bool
SequenceStatisticalOptions::corrMatrixCompute() const
{
  queso_deprecated();

  return m_ov.m_corrMatrixCompute;
}

void
SequenceStatisticalOptions::print(std::ostream& os) const
{
  queso_deprecated();

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
operator<<(std::ostream& os, const SequenceStatisticalOptions& obj)
{
  queso_deprecated();

  obj.print(os);

  return os;
}

}  // End namespace QUESO

#endif // ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
