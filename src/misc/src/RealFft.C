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

#include <queso/Fft.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_complex.h>

namespace QUESO {

// Math methods------------------------------------------
template <>
void
Fft<double>::forward(
  const std::vector<double>&                data,
        unsigned int                        fftSize,
        std::vector<std::complex<double> >& forwardResult)
{
  if (forwardResult.size() != fftSize) {
    forwardResult.resize(fftSize,std::complex<double>(0.,0.));
    std::vector<std::complex<double> >(forwardResult).swap(forwardResult);
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
  gsl_fft_real_workspace* realWkSpace = gsl_fft_real_workspace_alloc(fftSize);
  gsl_fft_real_wavetable* realWvTable = gsl_fft_real_wavetable_alloc(fftSize);

  gsl_fft_real_transform(&internalData[0],
                         1,
                         fftSize,
                         realWvTable,
                         realWkSpace);

  gsl_fft_real_wavetable_free(realWvTable);
  gsl_fft_real_workspace_free(realWkSpace);
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
    forwardResult[j] = std::complex<double>(realPartOfFFT,imagPartOfFFT);
  }

  return;
}
//-------------------------------------------------------
template <>
void
Fft<double>::inverse(
  const std::vector<double>&                data,
        unsigned int                        fftSize,
        std::vector<std::complex<double> >& inverseResult)
{
  if (inverseResult.size() != fftSize) {
    inverseResult.resize(fftSize,std::complex<double>(0.,0.));
    std::vector<std::complex<double> >(inverseResult).swap(inverseResult);
  }

  std::vector<double> internalData(2*fftSize,0.); // Yes, twice the fftSize
  unsigned int minSize = std::min((unsigned int) data.size(),fftSize);
  for (unsigned int j = 0; j < minSize; ++j) {
    internalData[2*j] = data[j];
  }

  //if (m_subDisplayFile()) {
  //  *m_subDisplayFile() << "In Fft<double>::inverse()"
  //                     << ": about to call gsl_fft_complex_inverse()"
  //                     << " with fftSize = "         << fftSize
  //                     << "; internalData.size() = " << internalData.size()
  //                     << std::endl;
  //}

  gsl_fft_complex_workspace* complexWkSpace = gsl_fft_complex_workspace_alloc(fftSize);
  gsl_fft_complex_wavetable* complexWvTable = gsl_fft_complex_wavetable_alloc(fftSize);

  gsl_fft_complex_inverse(&internalData[0],
                          1,
                          fftSize,
                          complexWvTable,
                          complexWkSpace);

  gsl_fft_complex_wavetable_free(complexWvTable);
  gsl_fft_complex_workspace_free(complexWkSpace);

  //if (m_subDisplayFile()) {
  //  *m_subDisplayFile() << "In Fft<double>::inverse()"
  //                     << ": returned from gsl_fft_complex_inverse()"
  //                     << " with fftSize = "          << fftSize
  //                     << "; inverseResult.size() = " << inverseResult.size()
  //                     << std::endl;
  //}

  for (unsigned int j = 0; j < fftSize; ++j) {
    inverseResult[j] = std::complex<double>(internalData[2*j],internalData[2*j+1]);
  }

  return;
}
#if 0
template <class T>
void
Fft<T>::allocTables(unsigned int fftSize)
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
Fft<T>::freeTables()
{
  if (m_fftSize != 0) {
  //gsl_fft_halfcomplex_wavetable_free(m_hcWvTable);
    gsl_fft_real_wavetable_free       (m_realWvTable);
    gsl_fft_real_workspace_free       (m_realWkSpace);
    m_fftSize = 0;
  }

  return;
}
#endif

}  // End namespace QUESO
