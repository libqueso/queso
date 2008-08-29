/* uq/libs/basic/inc/uqRealFft.C
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

#include <uqFft.h>

template <>
uqFftClass<double>::uqFftClass(
  const uqEnvironmentClass& env,
        unsigned int        fftSize)
  :
  m_env    (env),
  m_fftSize(fftSize)
{
  m_realWkSpace = gsl_fft_real_workspace_alloc       (fftSize);
  m_realWvTable = gsl_fft_real_wavetable_alloc       (fftSize);
//m_hcWvTable   = gsl_fft_halfcomplex_wavetable_alloc(fftSize);
}

template <>
uqFftClass<double>::~uqFftClass()
{
//gsl_fft_halfcomplex_wavetable_free(m_hcWvTable);
  gsl_fft_real_wavetable_free       (m_realWvTable);
  gsl_fft_real_workspace_free       (m_realWkSpace);
}

template <>
void
uqFftClass<double>::forward(
  const std::vector<double>&                data, 
        std::vector<std::complex<double> >& result)
{
  std::vector<double> internalData(data);
  result.resize(m_fftSize,std::complex<double >(0.,0.));

  //double sumOfAllTerms = 0.;
  //for (unsigned int j = 0; j < fftSize; ++j) {
  //  sumOfAllTerms += internalData[j];
  //}
  gsl_fft_real_transform(&internalData[0],
                         1,
                         m_fftSize,
                         m_realWvTable,
                         m_realWkSpace);
  //std::cout << "After FFT"
  //          << ", sumOfAllTerms = "          << sumOfAllTerms
  //          << ", sumOfAllTerms - dft[0] = " << sumOfAllTerms - internalData[0]
  //          << std::endl;

  unsigned int halfFFTSize = m_fftSize/2;
  double realPartOfFFT = 0.;
  double imagPartOfFFT = 0.;
  for (unsigned int j = 0; j < internalData.size(); ++j) {
    if (j == 0) {
      realPartOfFFT = internalData[j];
      imagPartOfFFT = 0.;
    }
    else if (j < halfFFTSize) {
      realPartOfFFT = internalData[j];
      imagPartOfFFT = internalData[m_fftSize-j];
    }
    else if (j == halfFFTSize) {
      realPartOfFFT = internalData[j];
      imagPartOfFFT = 0.;
    }
    else {
      realPartOfFFT =  internalData[m_fftSize-j];
      imagPartOfFFT = -internalData[j];
    }
    result[j] = std::complex<double>(realPartOfFFT,imagPartOfFFT);
  }

  return;
}

template <>
void
uqFftClass<double>::inverse(
  const std::vector<double>&                data, 
        std::vector<std::complex<double> >& result)
{
  return;
}

