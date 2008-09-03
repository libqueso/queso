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
uqFftClass<double>::uqFftClass(const uqEnvironmentClass& env)
  :
  m_env    (env),
  m_fftSize(0) // Yes, zero
{
}

template <>
uqFftClass<double>::~uqFftClass()
{
  freeTables();
}

template <>
void
uqFftClass<double>::forward(
  const std::vector<double>&                data, 
        unsigned int                        fftSize,
        std::vector<std::complex<double> >& result)
{
  if (result.size() != fftSize) {
    result.resize(fftSize,std::complex<double>(0.,0.));
    std::vector<std::complex<double> >(result).swap(result);
  }

  std::vector<double> internalData(fftSize,0.);
  unsigned int minSize = std::min((unsigned int) data.size(),fftSize);
  for (unsigned int j = 0; j < minSize; ++j) {
    internalData[j] = data[j];
  }

  //double sumOfAllTerms = 0.;
  //for (unsigned int j = 0; j < fftSize; ++j) {
  //  sumOfAllTerms += internalData[j];
  //}

  //allocTables(fftSize);
  m_realWkSpace = gsl_fft_real_workspace_alloc(fftSize);
  m_realWvTable = gsl_fft_real_wavetable_alloc(fftSize);

  gsl_fft_real_transform(&internalData[0],
                         1,
                         fftSize,
                         m_realWvTable,
                         m_realWkSpace);

  gsl_fft_real_wavetable_free(m_realWvTable);
  gsl_fft_real_workspace_free(m_realWkSpace);
  //freeTables();

  //std::cout << "After FFT"
  //          << ", sumOfAllTerms = "          << sumOfAllTerms
  //          << ", sumOfAllTerms - dft[0] = " << sumOfAllTerms - internalData[0]
  //          << std::endl;

  unsigned int halfFFTSize = fftSize/2;
  double realPartOfFFT = 0.;
  double imagPartOfFFT = 0.;
  for (unsigned int j = 0; j < internalData.size(); ++j) {
    if (j == 0) {
      realPartOfFFT = internalData[j];
      imagPartOfFFT = 0.;
    }
    else if (j < halfFFTSize) {
      realPartOfFFT = internalData[2*j-1];
      imagPartOfFFT = internalData[2*j  ];
    }
    else if (j == halfFFTSize) {
      realPartOfFFT = internalData[2*j-1];
      imagPartOfFFT = 0.;
    }
    else {
      realPartOfFFT =  internalData[2*(fftSize-j)-1];
      imagPartOfFFT = -internalData[2*(fftSize-j)  ];
    }
    result[j] = std::complex<double>(realPartOfFFT,imagPartOfFFT);
  }

  return;
}

template <>
void
uqFftClass<double>::inverse(
  const std::vector<double>&                data, 
        unsigned int                        fftSize,
        std::vector<std::complex<double> >& result)
{
  if (result.size() != fftSize) {
    result.resize(fftSize,std::complex<double>(0.,0.));
    std::vector<std::complex<double> >(result).swap(result);
  }

  std::vector<double> internalData(2*fftSize,0.);                          // Yes, twice the fftSize
  unsigned int minSize = 2 * std::min((unsigned int) data.size(),fftSize); // Yes, 2*
  for (unsigned int j = 0; j < minSize; ++j) {
    internalData[2*j] = data[j];
  }

  //allocTables(fftSize);
  m_complexWkSpace = gsl_fft_complex_workspace_alloc(fftSize);
  m_complexWvTable = gsl_fft_complex_wavetable_alloc(fftSize);

  gsl_fft_complex_inverse(&internalData[0],
                          1,
                          fftSize,
                          m_complexWvTable,
                          m_complexWkSpace);

  gsl_fft_complex_wavetable_free(m_complexWvTable);
  gsl_fft_complex_workspace_free(m_complexWkSpace);
  //freeTables();

  for (unsigned int j = 0; j < fftSize; ++j) {
    result[j] = std::complex<double>(internalData[2*j],internalData[2*j+1]);
  }

  return;
}

template <class T>
void
uqFftClass<T>::allocTables(unsigned int fftSize)
{
  if (m_fftSize != fftSize) {
    if (m_fftSize != 0) freeTables();
    m_fftSize = fftSize;
    m_realWkSpace = gsl_fft_real_workspace_alloc       (fftSize);
    m_realWvTable = gsl_fft_real_wavetable_alloc       (fftSize);
  //m_hcWvTable   = gsl_fft_halfcomplex_wavetable_alloc(fftSize);
  }

  return;
}

template <class T>
void
uqFftClass<T>::freeTables()
{
  if (m_fftSize != 0) {
  //gsl_fft_halfcomplex_wavetable_free(m_hcWvTable);
    gsl_fft_real_wavetable_free       (m_realWvTable);
    gsl_fft_real_workspace_free       (m_realWkSpace);
    m_fftSize = 0;
  }

  return;
}

