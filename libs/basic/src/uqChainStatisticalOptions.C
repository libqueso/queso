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

#include <uqChainStatisticalOptions.h>
#include <uqMiscellaneous.h>

uqChainStatisticalOptionsClass::uqChainStatisticalOptionsClass(
  const uqBaseEnvironmentClass& env,
  const std::string&        prefix)
  :
  m_env                     (env),
  m_prefix                  ((std::string)(prefix) + "stats_"),
  m_optionsDesc             (new po::options_description("Chain statistical options")),
  m_initialDiscardedPortions(0),//,0.),
  m_bmmRun                  (UQ_CHAIN_BMM_RUN_ODV),
  m_bmmLengths              (0),//,0),
  m_fftCompute              (UQ_CHAIN_FFT_COMPUTE_ODV),
  m_fftParamId              (UQ_CHAIN_FFT_PARAM_ID_ODV),
  m_fftSize                 (UQ_CHAIN_FFT_SIZE_ODV),
  m_fftTestInversion        (UQ_CHAIN_FFT_TEST_INVERSION_ODV),
  m_fftWrite                (UQ_CHAIN_FFT_WRITE_ODV),
  m_psdCompute              (UQ_CHAIN_PSD_COMPUTE_ODV),
  m_psdNumBlocks            (UQ_CHAIN_PSD_NUM_BLOCKS_ODV),
  m_psdHopSizeRatio         (UQ_CHAIN_PSD_HOP_SIZE_RATIO_ODV),
  m_psdParamId              (UQ_CHAIN_PSD_PARAM_ID_ODV),
  m_psdWrite                (UQ_CHAIN_PSD_WRITE_ODV),
  m_psdAtZeroCompute        (UQ_CHAIN_PSD_AT_ZERO_COMPUTE_ODV),
  m_psdAtZeroNumBlocks      (0),//,0),
  m_psdAtZeroHopSizeRatio   (UQ_CHAIN_PSD_AT_ZERO_HOP_SIZE_RATIO_ODV),
  m_psdAtZeroDisplay        (UQ_CHAIN_PSD_AT_ZERO_DISPLAY_ODV),
  m_psdAtZeroWrite          (UQ_CHAIN_PSD_AT_ZERO_WRITE_ODV),
  m_gewekeCompute           (UQ_CHAIN_GEWEKE_COMPUTE_ODV),
  m_gewekeNaRatio           (UQ_CHAIN_GEWEKE_NA_RATIO_ODV),
  m_gewekeNbRatio           (UQ_CHAIN_GEWEKE_NB_RATIO_ODV),
  m_gewekeDisplay           (UQ_CHAIN_GEWEKE_DISPLAY_ODV),
  m_gewekeWrite             (UQ_CHAIN_GEWEKE_WRITE_ODV),
  m_autoCorrComputeViaDef   (UQ_CHAIN_AUTO_CORR_COMPUTE_VIA_DEF_ODV),
  m_autoCorrComputeViaFft   (UQ_CHAIN_AUTO_CORR_COMPUTE_VIA_FFT_ODV),
  m_autoCorrSecondLag       (UQ_CHAIN_AUTO_CORR_SECOND_LAG_ODV),
  m_autoCorrLagSpacing      (UQ_CHAIN_AUTO_CORR_LAG_SPACING_ODV),
  m_autoCorrNumLags         (UQ_CHAIN_AUTO_CORR_NUM_LAGS_ODV),
  m_autoCorrDisplay         (UQ_CHAIN_AUTO_CORR_DISPLAY_ODV),
  m_autoCorrWrite           (UQ_CHAIN_AUTO_CORR_WRITE_ODV),
  m_histCompute             (UQ_CHAIN_HIST_COMPUTE_ODV),
  m_histNumInternalBins     (UQ_CHAIN_HIST_NUM_INTERNAL_BINS_ODV),
  m_cdfStaccCompute         (UQ_CHAIN_CDF_STACC_COMPUTE_ODV),
  m_cdfStaccNumEvalPositions(UQ_CHAIN_CDF_STACC_NUM_EVAL_POSITIONS_ODV),
  m_kdeCompute              (UQ_CHAIN_KDE_COMPUTE_ODV),
  m_kdeNumEvalPositions     (UQ_CHAIN_KDE_NUM_EVAL_POSITIONS_ODV),
  m_covMatrixCompute        (UQ_CHAIN_COV_MATRIX_COMPUTE_ODV),
  m_corrMatrixCompute       (UQ_CHAIN_CORR_MATRIX_COMPUTE_ODV)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering uqChainStatisticalOptions::constructor()"
                           << ", prefix = " << m_prefix
                           << std::endl;
  }

  m_option_help                      = m_prefix + "help";

  m_option_initialDiscardedPortions  = m_prefix + "initialDiscardedPortions";

  m_option_bmm_run                   = m_prefix + "bmm_run";
  m_option_bmm_lengths               = m_prefix + "bmm_lengths";
  m_option_bmm_display               = m_prefix + "bmm_display";
  m_option_bmm_write                 = m_prefix + "bmm_write";

  m_option_fft_compute               = m_prefix + "fft_compute";
  m_option_fft_paramId               = m_prefix + "fft_paramId";
  m_option_fft_size                  = m_prefix + "fft_size";
  m_option_fft_testInversion         = m_prefix + "fft_testInversion";
  m_option_fft_write                 = m_prefix + "fft_write";

  m_option_psd_compute               = m_prefix + "psd_compute";
  m_option_psd_numBlocks             = m_prefix + "psd_numBlocks";
  m_option_psd_hopSizeRatio          = m_prefix + "psd_hopSizeRatio";
  m_option_psd_paramId               = m_prefix + "psd_paramId";
  m_option_psd_write                 = m_prefix + "psd_write";

  m_option_psdAtZero_compute         = m_prefix + "psdAtZero_compute";
  m_option_psdAtZero_numBlocks       = m_prefix + "psdAtZero_numBlocks";
  m_option_psdAtZero_hopSizeRatio    = m_prefix + "psdAtZero_hopSizeRatio";
  m_option_psdAtZero_display         = m_prefix + "psdAtZero_display";
  m_option_psdAtZero_write           = m_prefix + "psdAtZero_write";

  m_option_geweke_compute            = m_prefix + "geweke_compute";
  m_option_geweke_naRatio            = m_prefix + "geweke_naRatio";
  m_option_geweke_nbRatio            = m_prefix + "geweke_nbRatio";
  m_option_geweke_display            = m_prefix + "geweke_display";
  m_option_geweke_write              = m_prefix + "geweke_write";

  m_option_autoCorr_computeViaDef    = m_prefix + "autoCorr_computeViaDef";
  m_option_autoCorr_computeViaFft    = m_prefix + "autoCorr_computeViaFft";
  m_option_autoCorr_secondLag        = m_prefix + "autoCorr_secondLag";
  m_option_autoCorr_lagSpacing       = m_prefix + "autoCorr_lagSpacing";
  m_option_autoCorr_numLags          = m_prefix + "autoCorr_numLags";
  m_option_autoCorr_display          = m_prefix + "autoCorr_display";
  m_option_autoCorr_write            = m_prefix + "autoCorr_write";

  m_option_hist_compute              = m_prefix + "hist_compute";
  m_option_hist_numInternalBins      = m_prefix + "hist_numInternalBins";

  m_option_cdfStacc_compute          = m_prefix + "cdfStacc_compute";
  m_option_cdfStacc_numEvalPositions = m_prefix + "cdfStacc_numEvalPositions";

  m_option_kde_compute               = m_prefix + "kde_compute";
  m_option_kde_numEvalPositions      = m_prefix + "kde_numEvalPositions";

  m_option_covMatrix_compute         = m_prefix + "covMatrix_compute";
  m_option_corrMatrix_compute        = m_prefix + "corrMatrix_compute";

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "After getting values of options with prefix '" << m_prefix
                           << "', state of uqChainStatisticalOptionsClass object is:"
                           << "\n" << *this
                           << std::endl;
  }

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving uqChainStatisticalOptions::constructor()"
                           << ", prefix = " << m_prefix
                           << std::endl;
  }
}

uqChainStatisticalOptionsClass::~uqChainStatisticalOptionsClass()
{
  if (m_optionsDesc) delete m_optionsDesc;
}

void
uqChainStatisticalOptionsClass::defineMyOptions(
  po::options_description& optionsDesc) const
{
  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                                                                    "produce help message for chain statistical options"             )
    (m_option_initialDiscardedPortions.c_str(),       po::value<std::string >()->default_value(UQ_CHAIN_INITIAL_DISCARDED_PORTIONS_ODV      ), "list of initial discarded portions for chain statistics"        )
    (m_option_bmm_run.c_str(),                        po::value<bool        >()->default_value(UQ_CHAIN_BMM_RUN_ODV                         ), "compute variance of sample mean with batch means method"        )
    (m_option_bmm_lengths.c_str(),                    po::value<std::string >()->default_value(UQ_CHAIN_BMM_LENGTHS_ODV                     ), "list of batch lenghts for BMM"                                  )
    (m_option_fft_compute.c_str(),                    po::value<bool        >()->default_value(UQ_CHAIN_FFT_COMPUTE_ODV                     ), "compute fft"                                                    )
    (m_option_fft_paramId.c_str(),                    po::value<unsigned int>()->default_value(UQ_CHAIN_FFT_PARAM_ID_ODV                    ), "parameter id for fft computations"                              )
    (m_option_fft_size.c_str(),                       po::value<unsigned int>()->default_value(UQ_CHAIN_FFT_SIZE_ODV                        ), "fft size"                                                       )
    (m_option_fft_testInversion.c_str(),              po::value<bool        >()->default_value(UQ_CHAIN_FFT_TEST_INVERSION_ODV              ), "test fft inversion"                                             )
    (m_option_fft_write.c_str(),                      po::value<bool        >()->default_value(UQ_CHAIN_FFT_WRITE_ODV                       ), "write fft"                                                      )
    (m_option_psd_compute.c_str(),                    po::value<bool        >()->default_value(UQ_CHAIN_PSD_COMPUTE_ODV                     ), "compute psd"                                                    )
    (m_option_psd_numBlocks.c_str(),                  po::value<unsigned int>()->default_value(UQ_CHAIN_PSD_NUM_BLOCKS_ODV                  ), "number of blocks for psd"                                       )
    (m_option_psd_hopSizeRatio.c_str(),               po::value<double      >()->default_value(UQ_CHAIN_PSD_HOP_SIZE_RATIO_ODV              ), "hop size ratio for psd"                                         )
    (m_option_psd_paramId.c_str(),                    po::value<unsigned int>()->default_value(UQ_CHAIN_PSD_PARAM_ID_ODV                    ), "parameter id for psd computations"                              )
    (m_option_psd_write.c_str(),                      po::value<bool        >()->default_value(UQ_CHAIN_PSD_WRITE_ODV                       ), "write psd"                                                      )
    (m_option_psdAtZero_compute.c_str(),              po::value<bool        >()->default_value(UQ_CHAIN_PSD_AT_ZERO_COMPUTE_ODV             ), "compute power spectral densities"                               )
    (m_option_psdAtZero_numBlocks.c_str(),            po::value<std::string >()->default_value(UQ_CHAIN_PSD_AT_ZERO_NUM_BLOCKS_ODV          ), "list of numbers of blocks for computation of psd at zero"       )
    (m_option_psdAtZero_hopSizeRatio.c_str(),         po::value<double      >()->default_value(UQ_CHAIN_PSD_AT_ZERO_HOP_SIZE_RATIO_ODV      ), "hop size ratio for psd at zero"                                 )
    (m_option_psdAtZero_display.c_str(),              po::value<bool        >()->default_value(UQ_CHAIN_PSD_AT_ZERO_DISPLAY_ODV             ), "display computed psd at frequency zero on screen"               )
    (m_option_psdAtZero_write.c_str(),                po::value<bool        >()->default_value(UQ_CHAIN_PSD_AT_ZERO_WRITE_ODV               ), "write computed psd at frequency zero to the output file"        )
    (m_option_geweke_compute.c_str(),                 po::value<bool        >()->default_value(UQ_CHAIN_GEWEKE_COMPUTE_ODV                  ), "compute Geweke coefficients"                                    )
    (m_option_geweke_naRatio.c_str(),                 po::value<double      >()->default_value(UQ_CHAIN_GEWEKE_NA_RATIO_ODV                 ), "ratio NA for Geweke"                                            ) 
    (m_option_geweke_nbRatio.c_str(),                 po::value<double      >()->default_value(UQ_CHAIN_GEWEKE_NB_RATIO_ODV                 ), "ratio NB for Geweke"                                            )
    (m_option_geweke_display.c_str(),                 po::value<bool        >()->default_value(UQ_CHAIN_GEWEKE_DISPLAY_ODV                  ), "display computed Geweke on screen"                              )
    (m_option_geweke_write.c_str(),                   po::value<bool        >()->default_value(UQ_CHAIN_GEWEKE_WRITE_ODV                    ), "write computed Geweke to the output file"                       )
    (m_option_autoCorr_computeViaDef.c_str(),         po::value<bool        >()->default_value(UQ_CHAIN_AUTO_CORR_COMPUTE_VIA_DEF_ODV       ), "compute correlations via definition"                            )
    (m_option_autoCorr_computeViaFft.c_str(),         po::value<bool        >()->default_value(UQ_CHAIN_AUTO_CORR_COMPUTE_VIA_FFT_ODV       ), "compute correlations via fft"                                   )
    (m_option_autoCorr_secondLag.c_str(),             po::value<unsigned int>()->default_value(UQ_CHAIN_AUTO_CORR_SECOND_LAG_ODV            ), "second lag for computation of autocorrelations"                 )
    (m_option_autoCorr_lagSpacing.c_str(),            po::value<unsigned int>()->default_value(UQ_CHAIN_AUTO_CORR_LAG_SPACING_ODV           ), "lag spacing for computation of autocorrelations"                )
    (m_option_autoCorr_numLags.c_str(),               po::value<unsigned int>()->default_value(UQ_CHAIN_AUTO_CORR_NUM_LAGS_ODV              ), "number of lags for computation of autocorrelations"             )
    (m_option_autoCorr_display.c_str(),               po::value<bool        >()->default_value(UQ_CHAIN_AUTO_CORR_DISPLAY_ODV               ), "display computed autocorrelations on the screen"                )
    (m_option_autoCorr_write.c_str(),                 po::value<bool        >()->default_value(UQ_CHAIN_AUTO_CORR_WRITE_ODV                 ), "write computed autocorrelations to the output file"             )
    (m_option_hist_compute.c_str(),                   po::value<bool        >()->default_value(UQ_CHAIN_HIST_COMPUTE_ODV                    ), "compute histograms"                                             )
    (m_option_hist_numInternalBins.c_str(),           po::value<unsigned int>()->default_value(UQ_CHAIN_HIST_NUM_INTERNAL_BINS_ODV          ), "number of internal bins"                                        )
    (m_option_cdfStacc_compute.c_str(),               po::value<bool        >()->default_value(UQ_CHAIN_CDF_STACC_COMPUTE_ODV               ), "compute histograms"                                             )
    (m_option_cdfStacc_numEvalPositions.c_str(),      po::value<unsigned int>()->default_value(UQ_CHAIN_CDF_STACC_NUM_EVAL_POSITIONS_ODV    ), "number of internal bins"                                        )
    (m_option_kde_compute.c_str(),                    po::value<bool        >()->default_value(UQ_CHAIN_KDE_COMPUTE_ODV                     ), "compute kernel density estimators"                              )
    (m_option_kde_numEvalPositions.c_str(),           po::value<unsigned int>()->default_value(UQ_CHAIN_KDE_NUM_EVAL_POSITIONS_ODV          ), "number of evaluation positions"                                 )
    (m_option_covMatrix_compute.c_str(),              po::value<bool        >()->default_value(UQ_CHAIN_COV_MATRIX_COMPUTE_ODV              ), "compute covariance matrix"                                      )
    (m_option_corrMatrix_compute.c_str(),             po::value<bool        >()->default_value(UQ_CHAIN_CORR_MATRIX_COMPUTE_ODV             ), "compute correlation matrix"                                     )
  ;

  return;
}

void
uqChainStatisticalOptionsClass::getMyOptionValues(
  po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << optionsDesc
                             << std::endl;
    }
  }

  if (m_env.allOptionsMap().count(m_option_initialDiscardedPortions.c_str())) {
    m_initialDiscardedPortions.clear();
    std::vector<double> tmpPortions(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_initialDiscardedPortions.c_str()].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpPortions);
    //if (m_env.subDisplayFile()) {
    //  *m_env.subDisplayFile() << "In uqChainStatisticalOptionsClass::getMyOptionValues(): percents = ";
    //  for (unsigned int i = 0; i < tmpPortions.size(); ++i) {
    //    *m_env.subDisplayFile() << " " << tmpPortions[i];
    //  }
    //  *m_env.subDisplayFile() << std::endl;
    //}

    if (tmpPortions.size() > 0) {
      m_initialDiscardedPortions.resize(tmpPortions.size(),0.);
      for (unsigned int i = 0; i < m_initialDiscardedPortions.size(); ++i) {
        m_initialDiscardedPortions[i] = tmpPortions[i];
      }
    }
  }

  if (m_env.allOptionsMap().count(m_option_bmm_run.c_str())) {
    m_bmmRun = m_env.allOptionsMap()[m_option_bmm_run.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_bmm_lengths.c_str())) {
    m_bmmLengths.clear();
    std::vector<double> tmpLengths(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_bmm_lengths.c_str()].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpLengths);
    //if (m_env.subDisplayFile()) {
    //  *m_env.subDisplayFile() << "In uqChainStatisticalOptionsClass::getMyOptionValues(): lengths for BMM = ";
    //  for (unsigned int i = 0; i < tmpLengths.size(); ++i) {
    //    *m_env.subDisplayFile() << " " << tmpLengths[i];
    //  }
    //  *m_env.subDisplayFile() << std::endl;
    //}

    if (tmpLengths.size() > 0) {
      m_bmmLengths.resize(tmpLengths.size(),0);
      for (unsigned int i = 0; i < m_bmmLengths.size(); ++i) {
        m_bmmLengths[i] = (unsigned int) tmpLengths[i];
      }
    }
  }

  if (m_env.allOptionsMap().count(m_option_fft_compute.c_str())) {
    m_fftCompute = m_env.allOptionsMap()[m_option_fft_compute.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_fft_paramId.c_str())) {
    m_fftParamId = m_env.allOptionsMap()[m_option_fft_paramId.c_str()].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_fft_size.c_str())) {
    m_fftSize = m_env.allOptionsMap()[m_option_fft_size.c_str()].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_fft_testInversion.c_str())) {
    m_fftTestInversion = m_env.allOptionsMap()[m_option_fft_testInversion.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_fft_write.c_str())) {
    m_fftWrite = m_env.allOptionsMap()[m_option_fft_write.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_psd_compute.c_str())) {
    m_psdCompute = m_env.allOptionsMap()[m_option_psd_compute.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_psd_numBlocks.c_str())) {
    m_psdNumBlocks = m_env.allOptionsMap()[m_option_psd_numBlocks.c_str()].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_psd_hopSizeRatio.c_str())) {
    m_psdHopSizeRatio = m_env.allOptionsMap()[m_option_psd_hopSizeRatio.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_psd_paramId.c_str())) {
    m_psdParamId = m_env.allOptionsMap()[m_option_psd_paramId.c_str()].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_psd_write.c_str())) {
    m_psdWrite = m_env.allOptionsMap()[m_option_psd_write.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_psdAtZero_compute.c_str())) {
    m_psdAtZeroCompute = m_env.allOptionsMap()[m_option_psdAtZero_compute.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_psdAtZero_numBlocks.c_str())) {
    m_psdAtZeroNumBlocks.clear();
    std::vector<double> tmpNumBlocks(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_psdAtZero_numBlocks.c_str()].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpNumBlocks);
    //if (m_env.subDisplayFile()) {
    //  *m_env.subDisplayFile() << "In uqChainStatisticalOptionsClass::getMyOptionValues(): numBlocks for psdAtZero = ";
    //  for (unsigned int i = 0; i < tmpNumBlocks.size(); ++i) {
    //    *m_env.subDisplayFile() << " " << numBlocks[i];
    //  }
    //  *m_env.subDisplayFile() << std::endl;
    //}

    if (tmpNumBlocks.size() > 0) {
      m_psdAtZeroNumBlocks.resize(tmpNumBlocks.size(),0);
      for (unsigned int i = 0; i < m_psdAtZeroNumBlocks.size(); ++i) {
        m_psdAtZeroNumBlocks[i] = (unsigned int) tmpNumBlocks[i];
      }
    }
  }

  if (m_env.allOptionsMap().count(m_option_psdAtZero_hopSizeRatio.c_str())) {
    m_psdAtZeroHopSizeRatio = m_env.allOptionsMap()[m_option_psdAtZero_hopSizeRatio.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_psdAtZero_display.c_str())) {
    m_psdAtZeroDisplay = m_env.allOptionsMap()[m_option_psdAtZero_display.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_psdAtZero_write.c_str())) {
    m_psdAtZeroWrite = m_env.allOptionsMap()[m_option_psdAtZero_write.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_geweke_compute.c_str())) {
    m_gewekeCompute = m_env.allOptionsMap()[m_option_geweke_compute.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_geweke_naRatio.c_str())) {
    m_gewekeNaRatio = m_env.allOptionsMap()[m_option_geweke_naRatio.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_geweke_nbRatio.c_str())) {
    m_gewekeNbRatio = m_env.allOptionsMap()[m_option_geweke_nbRatio.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_geweke_display.c_str())) {
    m_gewekeDisplay = m_env.allOptionsMap()[m_option_geweke_display.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_geweke_write.c_str())) {
    m_gewekeWrite = m_env.allOptionsMap()[m_option_geweke_write.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_autoCorr_computeViaDef.c_str())) {
    m_autoCorrComputeViaDef = m_env.allOptionsMap()[m_option_autoCorr_computeViaDef.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_autoCorr_computeViaFft.c_str())) {
    m_autoCorrComputeViaFft = m_env.allOptionsMap()[m_option_autoCorr_computeViaFft.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_autoCorr_secondLag.c_str())) {
    m_autoCorrSecondLag = m_env.allOptionsMap()[m_option_autoCorr_secondLag.c_str()].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_autoCorr_lagSpacing.c_str())) {
    m_autoCorrLagSpacing = m_env.allOptionsMap()[m_option_autoCorr_lagSpacing.c_str()].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_autoCorr_numLags.c_str())) {
    m_autoCorrNumLags = m_env.allOptionsMap()[m_option_autoCorr_numLags.c_str()].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_autoCorr_display.c_str())) {
    m_autoCorrDisplay = m_env.allOptionsMap()[m_option_autoCorr_display.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_autoCorr_write.c_str())) {
    m_autoCorrWrite = m_env.allOptionsMap()[m_option_autoCorr_write.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_hist_compute.c_str())) {
    m_histCompute = m_env.allOptionsMap()[m_option_hist_compute.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_hist_numInternalBins.c_str())) {
    m_histNumInternalBins = m_env.allOptionsMap()[m_option_hist_numInternalBins.c_str()].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_cdfStacc_compute.c_str())) {
    m_cdfStaccCompute = m_env.allOptionsMap()[m_option_cdfStacc_compute.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_cdfStacc_numEvalPositions.c_str())) {
    m_cdfStaccNumEvalPositions = m_env.allOptionsMap()[m_option_cdfStacc_numEvalPositions.c_str()].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_kde_compute.c_str())) {
    m_kdeCompute = m_env.allOptionsMap()[m_option_kde_compute.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_kde_numEvalPositions.c_str())) {
    m_kdeNumEvalPositions = m_env.allOptionsMap()[m_option_kde_numEvalPositions.c_str()].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_covMatrix_compute.c_str())) {
    m_covMatrixCompute = m_env.allOptionsMap()[m_option_covMatrix_compute.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_corrMatrix_compute.c_str())) {
    m_corrMatrixCompute = m_env.allOptionsMap()[m_option_corrMatrix_compute.c_str()].as<bool>();
  }

  return;
}

const std::vector<double>&
uqChainStatisticalOptionsClass::initialDiscardedPortions() const
{
  return m_initialDiscardedPortions;
}

bool
uqChainStatisticalOptionsClass::bmmRun() const
{
  return m_bmmRun;
}

const std::vector<unsigned int>&
uqChainStatisticalOptionsClass::bmmLengths() const
{
  return m_bmmLengths;
}

bool
uqChainStatisticalOptionsClass::bmmDisplay() const
{
  return m_bmmDisplay;
}

bool
uqChainStatisticalOptionsClass::bmmWrite() const
{
  return m_bmmWrite;
}

bool
uqChainStatisticalOptionsClass::fftCompute() const
{
  return m_fftCompute;
}

unsigned int
uqChainStatisticalOptionsClass::fftParamId() const
{
  return m_fftParamId;
}

unsigned int
uqChainStatisticalOptionsClass::fftSize() const
{
  return m_fftSize;
}

bool
uqChainStatisticalOptionsClass::fftTestInversion() const
{
  return m_fftTestInversion;
}

bool
uqChainStatisticalOptionsClass::fftWrite() const
{
  return m_fftWrite;
}

bool
uqChainStatisticalOptionsClass::psdCompute() const
{
  return m_psdCompute;
}

unsigned int
uqChainStatisticalOptionsClass::psdNumBlocks() const
{
  return m_psdNumBlocks;
}

double
uqChainStatisticalOptionsClass::psdHopSizeRatio() const
{
  return m_psdHopSizeRatio;
}

unsigned int
uqChainStatisticalOptionsClass::psdParamId() const
{
  return m_psdParamId;
}

bool
uqChainStatisticalOptionsClass::psdWrite() const
{
  return m_psdWrite;
}

bool
uqChainStatisticalOptionsClass::psdAtZeroCompute() const
{
  return m_psdAtZeroCompute;
}

const std::vector<unsigned int>&
uqChainStatisticalOptionsClass::psdAtZeroNumBlocks() const
{
  return m_psdAtZeroNumBlocks;
}

double
uqChainStatisticalOptionsClass::psdAtZeroHopSizeRatio() const
{
  return m_psdAtZeroHopSizeRatio;
}

bool
uqChainStatisticalOptionsClass::psdAtZeroDisplay() const
{
  return m_psdAtZeroDisplay;
}

bool
uqChainStatisticalOptionsClass::psdAtZeroWrite() const
{
  return m_psdAtZeroWrite;
}

bool
uqChainStatisticalOptionsClass::gewekeCompute() const
{
  return m_gewekeCompute;
}

double
uqChainStatisticalOptionsClass::gewekeNaRatio() const
{
  return m_gewekeNaRatio;
}

double
uqChainStatisticalOptionsClass::gewekeNbRatio() const
{
  return m_gewekeNbRatio;
}

bool
uqChainStatisticalOptionsClass::gewekeDisplay() const
{
  return m_gewekeDisplay;
}

bool
uqChainStatisticalOptionsClass::gewekeWrite() const
{
  return m_gewekeWrite;
}

bool
uqChainStatisticalOptionsClass::autoCorrComputeViaDef() const
{
  return m_autoCorrComputeViaDef;
}

bool
uqChainStatisticalOptionsClass::autoCorrComputeViaFft() const
{
  return m_autoCorrComputeViaFft;
}

unsigned int
uqChainStatisticalOptionsClass::autoCorrSecondLag() const
{
  return m_autoCorrSecondLag;
}

unsigned int
uqChainStatisticalOptionsClass::autoCorrLagSpacing() const
{
  return m_autoCorrLagSpacing;
}

unsigned int
uqChainStatisticalOptionsClass::autoCorrNumLags() const
{
  return m_autoCorrNumLags;
}

bool
uqChainStatisticalOptionsClass::autoCorrDisplay() const
{
  return m_autoCorrDisplay;
}

bool
uqChainStatisticalOptionsClass::autoCorrWrite() const
{
  return m_autoCorrWrite;
}

bool
uqChainStatisticalOptionsClass::histCompute() const
{
  return m_histCompute;
}

unsigned int
uqChainStatisticalOptionsClass::histNumInternalBins() const
{
  return m_histNumInternalBins;
}

bool
uqChainStatisticalOptionsClass::cdfStaccCompute() const
{
  return m_cdfStaccCompute;
}

unsigned int
uqChainStatisticalOptionsClass::cdfStaccNumEvalPositions() const
{
  return m_cdfStaccNumEvalPositions;
}

bool
uqChainStatisticalOptionsClass::kdeCompute() const
{
  return m_kdeCompute;
}

unsigned int
uqChainStatisticalOptionsClass::kdeNumEvalPositions() const
{
  return m_kdeNumEvalPositions;
}

bool
uqChainStatisticalOptionsClass::covMatrixCompute() const
{
  return m_covMatrixCompute;
}

bool
uqChainStatisticalOptionsClass::corrMatrixCompute() const
{
  return m_corrMatrixCompute;
}

void
uqChainStatisticalOptionsClass::print(std::ostream& os) const
{
  os << "\n" << m_option_initialDiscardedPortions << " = ";
  for (unsigned int i = 0; i < m_initialDiscardedPortions.size(); ++i) {
    os << m_initialDiscardedPortions[i] << " ";
  }
  os << "\n" << m_option_bmm_run     << " = " << m_bmmRun
     << "\n" << m_option_bmm_lengths << " = ";
  for (unsigned int i = 0; i < m_bmmLengths.size(); ++i) {
    os << m_bmmLengths[i] << " ";
  }
  os << "\n" << m_option_fft_compute         << " = " << m_fftCompute
     << "\n" << m_option_fft_paramId         << " = " << m_fftParamId
     << "\n" << m_option_fft_size            << " = " << m_fftSize
     << "\n" << m_option_fft_testInversion   << " = " << m_fftTestInversion
     << "\n" << m_option_fft_write           << " = " << m_fftWrite
     << "\n" << m_option_psd_compute         << " = " << m_psdCompute
     << "\n" << m_option_psd_paramId         << " = " << m_psdParamId
     << "\n" << m_option_psd_numBlocks       << " = " << m_psdNumBlocks
     << "\n" << m_option_psd_hopSizeRatio    << " = " << m_psdHopSizeRatio
     << "\n" << m_option_psd_write           << " = " << m_psdWrite
     << "\n" << m_option_psdAtZero_compute   << " = " << m_psdAtZeroCompute
     << "\n" << m_option_psdAtZero_numBlocks << " = ";
  for (unsigned int i = 0; i < m_psdAtZeroNumBlocks.size(); ++i) {
    os << m_psdAtZeroNumBlocks[i] << " ";
  }
  os << "\n" << m_option_psdAtZero_hopSizeRatio    << " = " << m_psdAtZeroHopSizeRatio
     << "\n" << m_option_psdAtZero_display         << " = " << m_psdAtZeroDisplay
     << "\n" << m_option_psdAtZero_write           << " = " << m_psdAtZeroWrite
     << "\n" << m_option_geweke_compute            << " = " << m_gewekeCompute
     << "\n" << m_option_geweke_naRatio            << " = " << m_gewekeNaRatio
     << "\n" << m_option_geweke_nbRatio            << " = " << m_gewekeNbRatio
     << "\n" << m_option_geweke_display            << " = " << m_gewekeDisplay
     << "\n" << m_option_geweke_write              << " = " << m_gewekeWrite
     << "\n" << m_option_autoCorr_computeViaDef    << " = " << m_autoCorrComputeViaDef
     << "\n" << m_option_autoCorr_computeViaFft    << " = " << m_autoCorrComputeViaFft
     << "\n" << m_option_autoCorr_secondLag        << " = " << m_autoCorrSecondLag
     << "\n" << m_option_autoCorr_lagSpacing       << " = " << m_autoCorrLagSpacing
     << "\n" << m_option_autoCorr_numLags          << " = " << m_autoCorrNumLags
     << "\n" << m_option_autoCorr_display          << " = " << m_autoCorrDisplay
     << "\n" << m_option_autoCorr_write            << " = " << m_autoCorrWrite
     << "\n" << m_option_hist_compute              << " = " << m_histCompute
     << "\n" << m_option_hist_numInternalBins      << " = " << m_histNumInternalBins
     << "\n" << m_option_cdfStacc_compute          << " = " << m_cdfStaccCompute
     << "\n" << m_option_cdfStacc_numEvalPositions << " = " << m_cdfStaccNumEvalPositions
     << "\n" << m_option_kde_compute               << " = " << m_kdeCompute
     << "\n" << m_option_kde_numEvalPositions      << " = " << m_kdeNumEvalPositions
     << "\n" << m_option_covMatrix_compute         << " = " << m_covMatrixCompute
     << "\n" << m_option_corrMatrix_compute        << " = " << m_corrMatrixCompute
     << std::endl;

  return;
}

std::ostream&
operator<<(std::ostream& os, const uqChainStatisticalOptionsClass& obj)
{
  obj.print(os);

  return os;
}
