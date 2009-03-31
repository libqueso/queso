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

#ifndef __UQ_CHAIN_STATISTICAL_OPTIONS_H__
#define __UQ_CHAIN_STATISTICAL_OPTIONS_H__

#include <uqEnvironment.h>

#define UQ_MCMC_INITIAL_DISCARDED_PORTIONS_ODV "100."
#define UQ_MCMC_BMM_RUN_ODV                    0
#define UQ_MCMC_BMM_LENGTHS_ODV                "0"
#define UQ_MCMC_BMM_DISPLAY_ODV                0
#define UQ_MCMC_BMM_WRITE_ODV                  0
#define UQ_MCMC_FFT_COMPUTE_ODV                0
#define UQ_MCMC_FFT_PARAM_ID_ODV               0
#define UQ_MCMC_FFT_SIZE_ODV                   2048
#define UQ_MCMC_FFT_TEST_INVERSION_ODV         0
#define UQ_MCMC_FFT_WRITE_ODV                  0
#define UQ_MCMC_PSD_COMPUTE_ODV                0
#define UQ_MCMC_PSD_NUM_BLOCKS_ODV             0
#define UQ_MCMC_PSD_HOP_SIZE_RATIO_ODV         0.
#define UQ_MCMC_PSD_PARAM_ID_ODV               0
#define UQ_MCMC_PSD_WRITE_ODV                  0
#define UQ_MCMC_PSD_AT_ZERO_COMPUTE_ODV        0
#define UQ_MCMC_PSD_AT_ZERO_NUM_BLOCKS_ODV     "8"
#define UQ_MCMC_PSD_AT_ZERO_HOP_SIZE_RATIO_ODV .5
#define UQ_MCMC_PSD_AT_ZERO_DISPLAY_ODV        0
#define UQ_MCMC_PSD_AT_ZERO_WRITE_ODV          0
#define UQ_MCMC_GEWEKE_COMPUTE_ODV             0
#define UQ_MCMC_GEWEKE_NA_RATIO_ODV            .1
#define UQ_MCMC_GEWEKE_NB_RATIO_ODV            .5
#define UQ_MCMC_GEWEKE_DISPLAY_ODV             0
#define UQ_MCMC_GEWEKE_WRITE_ODV               0
#define UQ_MCMC_AUTO_CORR_COMPUTE_VIA_DEF_ODV  0
#define UQ_MCMC_AUTO_CORR_COMPUTE_VIA_FFT_ODV  0
#define UQ_MCMC_AUTO_CORR_SECOND_LAG_ODV       0
#define UQ_MCMC_AUTO_CORR_LAG_SPACING_ODV      0
#define UQ_MCMC_AUTO_CORR_NUM_LAGS_ODV         0
#define UQ_MCMC_AUTO_CORR_DISPLAY_ODV          0
#define UQ_MCMC_AUTO_CORR_WRITE_ODV            0
#define UQ_MCMC_HIST_COMPUTE_ODV               0
#define UQ_MCMC_HIST_NUM_INTERNAL_BINS_ODV     100
#define UQ_MCMC_KDE_COMPUTE_ODV                0
#define UQ_MCMC_KDE_NUM_EVAL_POSITIONS_ODV     100

/*! A templated class that stores statistical options for a chain.
 */
class uqChainStatisticalOptionsClass
{
public:
  uqChainStatisticalOptionsClass(const uqBaseEnvironmentClass& env,     /*! The QUESO toolkit environment.                 */
                                 const std::string&        prefix); /*! Prefix for reading parameters from input file. */
 ~uqChainStatisticalOptionsClass();

  const std::vector<double>&       initialDiscardedPortions() const;

        bool                       bmmRun    () const;
  const std::vector<unsigned int>& bmmLengths() const;
        bool                       bmmDisplay() const;
        bool                       bmmWrite  () const;

        bool                       fftCompute      () const;
        unsigned int               fftParamId      () const;
        unsigned int               fftSize         () const;
        bool                       fftTestInversion() const;
        bool                       fftWrite        () const;

        bool                       psdCompute     () const;
        unsigned int               psdNumBlocks   () const;
        double                     psdHopSizeRatio() const;
        unsigned int               psdParamId     () const;
        bool                       psdWrite       () const;

        bool                       psdAtZeroCompute     () const;
  const std::vector<unsigned int>& psdAtZeroNumBlocks   () const;
        double                     psdAtZeroHopSizeRatio() const;
        bool                       psdAtZeroDisplay     () const;
        bool                       psdAtZeroWrite       () const;

        bool                       gewekeCompute() const;
        double                     gewekeNaRatio() const;
        double                     gewekeNbRatio() const;
        bool                       gewekeDisplay() const;
        bool                       gewekeWrite  () const;

        bool                       autoCorrComputeViaDef() const;
        bool                       autoCorrComputeViaFft() const;
        unsigned int               autoCorrSecondLag    () const;
        unsigned int               autoCorrLagSpacing   () const;
        unsigned int               autoCorrNumLags      () const;
        bool                       autoCorrDisplay      () const;
        bool                       autoCorrWrite        () const;

        bool                       histCompute        () const;
        unsigned int               histNumInternalBins() const;

        bool                       kdeCompute         () const;
        unsigned int               kdeNumEvalPositions() const;

  void print                    (std::ostream& os) const;

private:
  void   defineMyOptions        (po::options_description& optionsDesc) const;
  void   getMyOptionValues      (po::options_description& optionsDesc);

  const uqBaseEnvironmentClass& m_env;
        std::string         m_prefix;

  std::string m_option_help;
  std::string m_option_initialDiscardedPortions;
  std::string m_option_bmm_run;
  std::string m_option_bmm_lengths;
  std::string m_option_bmm_display;
  std::string m_option_bmm_write;
  std::string m_option_fft_compute;
  std::string m_option_fft_paramId;
  std::string m_option_fft_size;
  std::string m_option_fft_testInversion;
  std::string m_option_fft_write;
  std::string m_option_psd_compute;
  std::string m_option_psd_numBlocks;
  std::string m_option_psd_hopSizeRatio;
  std::string m_option_psd_paramId;
  std::string m_option_psd_write;
  std::string m_option_psdAtZero_compute;
  std::string m_option_psdAtZero_numBlocks;
  std::string m_option_psdAtZero_hopSizeRatio;
  std::string m_option_psdAtZero_display;
  std::string m_option_psdAtZero_write;
  std::string m_option_geweke_compute;
  std::string m_option_geweke_naRatio;
  std::string m_option_geweke_nbRatio;
  std::string m_option_geweke_display;
  std::string m_option_geweke_write;
  std::string m_option_autoCorr_computeViaDef;
  std::string m_option_autoCorr_computeViaFft;
  std::string m_option_autoCorr_secondLag;
  std::string m_option_autoCorr_lagSpacing;
  std::string m_option_autoCorr_numLags;
  std::string m_option_autoCorr_display;
  std::string m_option_autoCorr_write;
  std::string m_option_hist_compute;
  std::string m_option_hist_numInternalBins;
  std::string m_option_kde_compute;
  std::string m_option_kde_numEvalPositions;

  po::options_description*  m_optionsDesc;

  std::vector<double>       m_initialDiscardedPortions;

  bool                      m_bmmRun;
  std::vector<unsigned int> m_bmmLengths;
  bool                      m_bmmDisplay;
  bool                      m_bmmWrite;

  bool                      m_fftCompute;
  unsigned int              m_fftParamId;
  unsigned int              m_fftSize;
  bool                      m_fftTestInversion;
  bool                      m_fftWrite;

  bool                      m_psdCompute;
  unsigned int              m_psdNumBlocks;
  double                    m_psdHopSizeRatio;
  unsigned int              m_psdParamId;
  bool                      m_psdWrite;

  bool                      m_psdAtZeroCompute;
  std::vector<unsigned int> m_psdAtZeroNumBlocks;
  double                    m_psdAtZeroHopSizeRatio;
  bool                      m_psdAtZeroDisplay;
  bool                      m_psdAtZeroWrite;

  bool                      m_gewekeCompute;
  double                    m_gewekeNaRatio;
  double                    m_gewekeNbRatio;
  bool                      m_gewekeDisplay;
  bool                      m_gewekeWrite;

  bool                      m_autoCorrComputeViaDef;
  bool                      m_autoCorrComputeViaFft;
  unsigned int              m_autoCorrSecondLag;
  unsigned int              m_autoCorrLagSpacing;
  unsigned int              m_autoCorrNumLags;
  bool                      m_autoCorrDisplay;
  bool                      m_autoCorrWrite;

  bool                      m_histCompute;
  unsigned int              m_histNumInternalBins;

  bool                      m_kdeCompute;
  unsigned int              m_kdeNumEvalPositions;
};

std::ostream& operator<<(std::ostream& os, const uqChainStatisticalOptionsClass& obj);

#endif // __UQ_CHAIN_STATISTICAL_OPTIONS_H__

