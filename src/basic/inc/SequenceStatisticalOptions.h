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

#ifndef UQ_SEQUENCE_STATISTICAL_OPTIONS_H
#define UQ_SEQUENCE_STATISTICAL_OPTIONS_H

#include <queso/Defines.h>

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS

#include <queso/Environment.h>
#include <queso/BoostInputOptionsParser.h>

#define UQ_SEQUENCE_INITIAL_DISCARDED_PORTIONS_ODV   "0."
#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
#define UQ_SEQUENCE_MEAN_MONITOR_PERIOD_ODV          0
#define UQ_SEQUENCE_BMM_RUN_ODV                      0
#define UQ_SEQUENCE_BMM_LENGTHS_ODV                  "0"
#define UQ_SEQUENCE_BMM_DISPLAY_ODV                  0
#define UQ_SEQUENCE_BMM_WRITE_ODV                    0
#define UQ_SEQUENCE_FFT_COMPUTE_ODV                  0
#define UQ_SEQUENCE_FFT_PARAM_ID_ODV                 0
#define UQ_SEQUENCE_FFT_SIZE_ODV                     2048
#define UQ_SEQUENCE_FFT_TEST_INVERSION_ODV           0
#define UQ_SEQUENCE_FFT_WRITE_ODV                    0
#define UQ_SEQUENCE_PSD_COMPUTE_ODV                  0
#define UQ_SEQUENCE_PSD_NUM_BLOCKS_ODV               8
#define UQ_SEQUENCE_PSD_HOP_SIZE_RATIO_ODV           0.
#define UQ_SEQUENCE_PSD_PARAM_ID_ODV                 0
#define UQ_SEQUENCE_PSD_WRITE_ODV                    0
#define UQ_SEQUENCE_PSD_AT_ZERO_COMPUTE_ODV          0
#define UQ_SEQUENCE_PSD_AT_ZERO_NUM_BLOCKS_ODV       "8"
#define UQ_SEQUENCE_PSD_AT_ZERO_HOP_SIZE_RATIO_ODV   .5
#define UQ_SEQUENCE_PSD_AT_ZERO_DISPLAY_ODV          0
#define UQ_SEQUENCE_PSD_AT_ZERO_WRITE_ODV            0
#define UQ_SEQUENCE_GEWEKE_COMPUTE_ODV               0
#define UQ_SEQUENCE_GEWEKE_NA_RATIO_ODV              .1
#define UQ_SEQUENCE_GEWEKE_NB_RATIO_ODV              .5
#define UQ_SEQUENCE_GEWEKE_DISPLAY_ODV               0
#define UQ_SEQUENCE_GEWEKE_WRITE_ODV                 0
#define UQ_SEQUENCE_MEAN_STACC_COMPUTE_ODV           0
#define UQ_SEQUENCE_HIST_COMPUTE_ODV                 0
#define UQ_SEQUENCE_HIST_NUM_INTERNAL_BINS_ODV       100
#define UQ_SEQUENCE_CDF_STACC_COMPUTE_ODV            0
#define UQ_SEQUENCE_CDF_STACC_NUM_EVAL_POSITIONS_ODV 50
#endif
#define UQ_SEQUENCE_AUTO_CORR_COMPUTE_VIA_DEF_ODV    0
#define UQ_SEQUENCE_AUTO_CORR_COMPUTE_VIA_FFT_ODV    0
#define UQ_SEQUENCE_AUTO_CORR_SECOND_LAG_ODV         0
#define UQ_SEQUENCE_AUTO_CORR_LAG_SPACING_ODV        0
#define UQ_SEQUENCE_AUTO_CORR_NUM_LAGS_ODV           0
#define UQ_SEQUENCE_AUTO_CORR_DISPLAY_ODV            0
#define UQ_SEQUENCE_AUTO_CORR_WRITE_ODV              0
#define UQ_SEQUENCE_KDE_COMPUTE_ODV                  0
#define UQ_SEQUENCE_KDE_NUM_EVAL_POSITIONS_ODV       100
#define UQ_SEQUENCE_COV_MATRIX_COMPUTE_ODV           0
#define UQ_SEQUENCE_CORR_MATRIX_COMPUTE_ODV          0

namespace boost {
  namespace program_options {
    class options_description;
  }
}

/*!\file SequenceStatisticalOptions.h
 * \brief A templated class that stores default statistical options
 *
 * \class SsOptionsValues
 * \brief A templated class that stores default statistical options for a sequence of vectors, e.g.
 *    a Markov chain, a Monte Carlo input sequence, or a Monte Carlo output sequence.
 */

class SsOptionsValues
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  /*! It assigns to the variables the pre-defined options for a sequence of data (scalars; vectors).*/
  SsOptionsValues            ();

  //! Prefix constructor.
  /*! Uses the prefix to read options from an input file. */
  SsOptionsValues(const BaseEnvironment * env, const char * prefix);

  //! Copy  constructor.
  /*! It assigns to \c this' variables, the same values of the variable of \c src.*/
  SsOptionsValues            (const SsOptionsValues& src);

  //! Destructor.
  virtual ~SsOptionsValues            ();
  //@}

  //! @name Set methods
  //@{
  //! Assignment operator; it copies \c rhs to \c this.
  SsOptionsValues& operator= (const SsOptionsValues& rhs);
  //@}

  //! @name Public attributes
  //@{

  std::string               m_prefix;

  //! Stores the initial  discarded portion of the chain.
  std::vector<double>       m_initialDiscardedPortions;

  //! Whether or not compute autocorrelation via definition.
  bool                      m_autoCorrComputeViaDef;

   //! Whether or not compute autocorrelation via FFT.
  bool                      m_autoCorrComputeViaFft;

  //! Second lag of the autocorrelation.
  unsigned int              m_autoCorrSecondLag;

  //! Lag spacing of the autocorrelation.
  unsigned int              m_autoCorrLagSpacing;

  //! Number of lags of the autocorrelation.
  unsigned int              m_autoCorrNumLags;

  //! Whether or not display autocorrelation.
  bool                      m_autoCorrDisplay;

  //! Whether or not write autocorrelation to file.
  bool                      m_autoCorrWrite;

  //! Whether or not compute kernel density estimate (kde).
  bool                      m_kdeCompute;

  //! Number of positions to evaluate kde.
  unsigned int              m_kdeNumEvalPositions;

  //! Whether or not compute covariance matrix.
  bool                      m_covMatrixCompute;

  //! Whether or not compute correlation matrix.
  bool                      m_corrMatrixCompute;

#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
  unsigned int              m_meanMonitorPeriod;

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

  bool                      m_meanStaccCompute;

  bool                      m_histCompute;
  unsigned int              m_histNumInternalBins;

  bool                      m_cdfStaccCompute;
  unsigned int              m_cdfStaccNumEvalPositions;
#endif
  //@}
  // end public attributes
private:
  BoostInputOptionsParser * m_parser;

  std::string                   m_option_help;
  std::string                   m_option_initialDiscardedPortions;

  std::string                   m_option_autoCorr_computeViaDef;
  std::string                   m_option_autoCorr_computeViaFft;
  std::string                   m_option_autoCorr_secondLag;
  std::string                   m_option_autoCorr_lagSpacing;
  std::string                   m_option_autoCorr_numLags;
  std::string                   m_option_autoCorr_display;
  std::string                   m_option_autoCorr_write;
  std::string                   m_option_kde_compute;
  std::string                   m_option_kde_numEvalPositions;
  std::string                   m_option_covMatrix_compute;
  std::string                   m_option_corrMatrix_compute;

#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
  std::string                   m_option_mean_monitorPeriod;
  std::string                   m_option_bmm_run;
  std::string                   m_option_bmm_lengths;
  std::string                   m_option_bmm_display;
  std::string                   m_option_bmm_write;
  std::string                   m_option_fft_compute;
  std::string                   m_option_fft_paramId;
  std::string                   m_option_fft_size;
  std::string                   m_option_fft_testInversion;
  std::string                   m_option_fft_write;
  std::string                   m_option_psd_compute;
  std::string                   m_option_psd_numBlocks;
  std::string                   m_option_psd_hopSizeRatio;
  std::string                   m_option_psd_paramId;
  std::string                   m_option_psd_write;
  std::string                   m_option_psdAtZero_compute;
  std::string                   m_option_psdAtZero_numBlocks;
  std::string                   m_option_psdAtZero_hopSizeRatio;
  std::string                   m_option_psdAtZero_display;
  std::string                   m_option_psdAtZero_write;
  std::string                   m_option_geweke_compute;
  std::string                   m_option_geweke_naRatio;
  std::string                   m_option_geweke_nbRatio;
  std::string                   m_option_geweke_display;
  std::string                   m_option_geweke_write;
  std::string                   m_option_meanStacc_compute;
  std::string                   m_option_hist_compute;
  std::string                   m_option_hist_numInternalBins;
  std::string                   m_option_cdfStacc_compute;
  std::string                   m_option_cdfStacc_numEvalPositions;
#endif

  //! Copies the option values from \c src to \c this.
  void copy(const SsOptionsValues& src);
};

//------------------------------------------------------------------------------------------------

/*!\class SequenceStatisticalOptions
 * \brief  A templated class that stores statistical options (optionally read from an input file)
 *
 *  A templated class that stores statistical options for a sequence of vectors, e.g.
 *    a Markov chain, a Monte Carlo input sequence, or a Monte Carlo output sequence.
 */


class SequenceStatisticalOptions
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor: reads options from the input file.
  /*! It assigns to the variables the identifying strings (prefixes) for a sequence of data
   * (scalars; vectors) which had been read from an input file.*/
  SequenceStatisticalOptions(const BaseEnvironment& env,
                                    const std::string&            prefix);

  //! Constructor: with alternative option values.
  /*! In this constructor, the input options are given by \c alternativeOptionsValues, thus, they
   * are not read from an input file.*/
  SequenceStatisticalOptions(const BaseEnvironment& env,
                                    const std::string&            prefix,
                                    const SsOptionsValues& alternativeOptionsValues);
  //! Destructor
  ~SequenceStatisticalOptions();
  //@}

  //! @name Statistical methods
  //@{

  //! Finds the initially discarded portion of the chain. Access to private attribute m_initialDiscardedPortions
  const std::vector<double>&       initialDiscardedPortions() const;

  //! Compute autocorrelation via definition.    Access to private attribute m_autoCorrComputeViaDef
  bool                       autoCorrComputeViaDef() const;

  //! Compute autocorrelation via FFT. Access to private attribute m_autoCorrComputeViaFft
  bool                       autoCorrComputeViaFft() const;

  //! Returns the second lag of the autocorrelation. Access to private attribute m_autoCorrSecondLag
  unsigned int               autoCorrSecondLag    () const;

  //! Returns the spacing of the autocorrelation. Access to private attribute m_autoCorrLagSpacing
  unsigned int               autoCorrLagSpacing   () const;

  //! Returns the number of lags of the autocorrelation. Access to private attribute m_autoCorrNumLags
  unsigned int               autoCorrNumLags      () const;

  //! Displays autocorrelation. Access to private attribute m_autoCorrDisplay
  bool                       autoCorrDisplay      () const;

  //! Writes autocorrelation. Access to private attribute m_autoCorrWrite
  bool                       autoCorrWrite        () const;

  //! Computes KDE. Access to private attribute m_kdeCompute
  bool                       kdeCompute         () const;

  //! Returns number of evaluation positions for KDE. Access to private attribute m_kdeNumEvalPositions
  unsigned int               kdeNumEvalPositions() const;

  //! Finds the covariance matrix. Access to private attribute m_covMatrixCompute
  bool                       covMatrixCompute () const;

  //! Finds the correlation matrix. Access to private attribute m_corrMatrixCompute
  bool                       corrMatrixCompute() const;
  //@}

  //! @name I/O method
  //@{
  //! Prints the initial discarded portion of the chain; and, optionally, other attributes of the chain.
  void                       print(std::ostream& os) const;
  //@}

  //! @name Public attribute
  //@{
  SsOptionsValues     m_ov;
  //@}

#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS

        unsigned int               meanMonitorPeriod() const;

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

        bool                       meanStaccCompute() const;

        bool                       histCompute        () const;
        unsigned int               histNumInternalBins() const;

        bool                       cdfStaccCompute         () const;
        unsigned int               cdfStaccNumEvalPositions() const;
#endif

private:
  //! Defines the options for the chain
  void   defineMyOptions  (boost::program_options::options_description& optionsDesc) const;

  //! Reads the chain options
  void   getMyOptionValues(boost::program_options::options_description& optionsDesc);

  std::string                   m_prefix;
  const BaseEnvironment& m_env;
  boost::program_options::options_description*      m_optionsDesc;

  std::string                   m_option_help;
  std::string                   m_option_initialDiscardedPortions;

  std::string                   m_option_autoCorr_computeViaDef;
  std::string                   m_option_autoCorr_computeViaFft;
  std::string                   m_option_autoCorr_secondLag;
  std::string                   m_option_autoCorr_lagSpacing;
  std::string                   m_option_autoCorr_numLags;
  std::string                   m_option_autoCorr_display;
  std::string                   m_option_autoCorr_write;
  std::string                   m_option_kde_compute;
  std::string                   m_option_kde_numEvalPositions;
  std::string                   m_option_covMatrix_compute;
  std::string                   m_option_corrMatrix_compute;

#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
  std::string                   m_option_mean_monitorPeriod;
  std::string                   m_option_bmm_run;
  std::string                   m_option_bmm_lengths;
  std::string                   m_option_bmm_display;
  std::string                   m_option_bmm_write;
  std::string                   m_option_fft_compute;
  std::string                   m_option_fft_paramId;
  std::string                   m_option_fft_size;
  std::string                   m_option_fft_testInversion;
  std::string                   m_option_fft_write;
  std::string                   m_option_psd_compute;
  std::string                   m_option_psd_numBlocks;
  std::string                   m_option_psd_hopSizeRatio;
  std::string                   m_option_psd_paramId;
  std::string                   m_option_psd_write;
  std::string                   m_option_psdAtZero_compute;
  std::string                   m_option_psdAtZero_numBlocks;
  std::string                   m_option_psdAtZero_hopSizeRatio;
  std::string                   m_option_psdAtZero_display;
  std::string                   m_option_psdAtZero_write;
  std::string                   m_option_geweke_compute;
  std::string                   m_option_geweke_naRatio;
  std::string                   m_option_geweke_nbRatio;
  std::string                   m_option_geweke_display;
  std::string                   m_option_geweke_write;
  std::string                   m_option_meanStacc_compute;
  std::string                   m_option_hist_compute;
  std::string                   m_option_hist_numInternalBins;
  std::string                   m_option_cdfStacc_compute;
  std::string                   m_option_cdfStacc_numEvalPositions;
#endif

};

std::ostream& operator<<(std::ostream& os, const SequenceStatisticalOptions& obj);

#endif // ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS

#endif // UQ_SEQUENCE_STATISTICAL_OPTIONS_H

