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

#include <uqFft.h>
#include <gsl/gsl_fft_complex.h>

template <>
void
uqFftClass<std::complex<double> >::forward(
  const std::vector<std::complex<double> >& data, 
        unsigned int                        fftSize,
        std::vector<std::complex<double> >& forwardResult)
{
  exit(1);
#if 0
  if (forwardResult.size() != fftSize) {
    forwardResult.resize(fftSize,std::complex<double>(0.,0.));
    std::vector<std::complex<double> >(forwardResult).swap(forwardResult);
  }

  std::vector<double> internalData(fftSize,0.);
  unsigned int minSize = std::min((unsigned int) data.size(),fftSize);
  for (unsigned int j = 0; j < minSize; ++j) {
    internalData[j] = data[j];
  }

  gsl_fft_real_workspace* realWkSpace = gsl_fft_real_workspace_alloc(fftSize);
  gsl_fft_real_wavetable* realWvTable = gsl_fft_real_wavetable_alloc(fftSize);

  gsl_fft_real_transform(&internalData[0],
                         1,
                         fftSize,
                         realWvTable,
                         realWkSpace);

  gsl_fft_real_wavetable_free(realWvTable);
  gsl_fft_real_workspace_free(realWkSpace);

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
    forwardResult[j] = std::complex<double>(realPartOfFFT,imagPartOfFFT);
  }
#endif
  return;
}

template <>
void
uqFftClass<std::complex<double> >::inverse(
  const std::vector<std::complex<double> >& data, 
        unsigned int                        fftSize,
        std::vector<std::complex<double> >& inverseResult)
{
  if (inverseResult.size() != fftSize) {
    inverseResult.resize(fftSize,std::complex<double>(0.,0.));
    std::vector<std::complex<double> >(inverseResult).swap(inverseResult);
  }

  std::vector<double> internalData(2*fftSize,0.);                          // Yes, twice the fftSize
  unsigned int minSize = 2 * std::min((unsigned int) data.size(),fftSize); // Yes, 2*
  for (unsigned int j = 0; j < minSize; ++j) {
    internalData[2*j  ] = data[j].real();
    internalData[2*j+1] = data[j].imag();
  }

  gsl_fft_complex_workspace* complexWkSpace = gsl_fft_complex_workspace_alloc(fftSize);
  gsl_fft_complex_wavetable* complexWvTable = gsl_fft_complex_wavetable_alloc(fftSize);

  gsl_fft_complex_inverse(&internalData[0],
                          1,
                          fftSize,
                          complexWvTable,
                          complexWkSpace);

  gsl_fft_complex_wavetable_free(complexWvTable);
  gsl_fft_complex_workspace_free(complexWkSpace);

  for (unsigned int j = 0; j < fftSize; ++j) {
    inverseResult[j] = std::complex<double>(internalData[2*j],internalData[2*j+1]);
  }

  return;
}
