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

#ifndef UQ_FFT_H
#define UQ_FFT_H

#include <queso/Environment.h>
#include <vector>
#include <complex>

namespace QUESO {

/*! \file Fft.h
    \brief Matrix class.
*/

/*! \class Fft
    \brief Class for a Fast Fourier Transform (FFT) algorithm.

    This class implements a Fast Fourier Transform (FFT) algorithm.
    Fast Fourier Transforms are efficient algorithms for calculating the discrete Fourier
    transform (DFT),
    \f[ x_j = \sum_{k=0}^{N-1} z_k \exp(-2\pi i j k / N) \f]
    and its inverse.
    The DFT usually arises as an approximation to the continuous Fourier transform when
    functions are sampled at discrete intervals in space or time. The naive evaluation of
    the discrete Fourier transform is a matrix-vector multiplication \f$ Ab \f$. A general
    matrix-vector multiplication takes \f$ O(N^2)\f$ operations for \f$ N \f$ data-points.
    Fast Fourier transform algorithms use a divide-and-conquer strategy to factorize the
    matrix \f$ A \f$ into smaller sub-matrices, corresponding to the integer factors of the
    length \f$ N \f$. If \f$ N \f$ can be factorized into a product of integers
    \f$ f_1 f_2 ... f_n \f$ then the DFT can be computed in \f$ O(N \sum f_i) \f$ operations.
    For a radix-2 FFT this gives an operation count of \f$ O(N \log_2 N)\f$.

    \todo: Implement Forward Fourier Transform for Complex data.

*/
/* If the function to be transformed is not
    harmonically related to the sampling frequency, the response of an FFT looks like a sinc
    function (although the integrated power is still correct). Aliasing (also known as leakage)
    can be reduced by apodization using an apodization function. However, aliasing reduction
    is at the expense of broadening the spectral response.*/

template <class T>
class Fft
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default Constructor
  Fft(const BaseEnvironment& env);

  //! Destructor
  ~Fft();
  //@}

  //! @name Mathematical methods
  //@{
  //! Calculates the forward Fourier transform (for real data. TODO: complex data).
  /*! This function uses GSL function 'gsl_fft_real_transform' to compute the FFT of
   * \c data, a real array (of time-ordered real data) of length \c fftSize, using a
   * mixed radix decimation-in-frequency algorithm. There is no restriction on the length
   * \c fftSize. Efficient modules are provided for subtransforms of length 2, 3, 4 and 5.
   * Any remaining factors are computed with a slow, \f$ O(\c fftSize ^2) \f$, general-
   * \c fftSize module. The caller must supply a wavetable containing trigonometric lookup
   * tables and a workspace work.\n The definition of the forward Fourier transform,
   * \f$ x = FFT(z) \f$ of size \f$ N \f$ is:
   * \f[ x_j = \sum_{k=0}^{N-1} z_k \exp(-2\pi i j k / N).  \f]*/
  void forward(const std::vector<T>&                     data,
                     unsigned int                        fftSize,
                     std::vector<std::complex<double> >& result);

  //! Calculates the inverse Fourier transform for real and complex data.
  /*! This function uses GSL function 'gsl_fft_complex_inverse' to compute the FFT of
   * \c data, a real array (of time-ordered real data) of length \c fftSize, using a
   * mixed radix decimation-in-frequency algorithm. There is no restriction on the length
   * \c fftSize. Efficient modules are provided for subtransforms of length 2, 3, 4 and 5.
   * Any remaining factors are computed with a slow, \f$ O(fftSize^2) \f$, general-
   * \c fftSize module. The caller must supply a wavetable containing trigonometric lookup
   * tables and a workspace work.\n The definition of the inverse Fourier transform,
   * \f$ x = IFFT(z)\f$ of size \f$ N \f$ is:
   * \f[ z_j = {1 \over N} \sum_{k=0}^{N-1} x_k \exp(2\pi i j k / N).\f]
   * The factor of \f$ 1/N \f$ makes this a true inverse. */
  void inverse(const std::vector<T>&                     data,
                     unsigned int                        fftSize,
                     std::vector<std::complex<double> >& result);
  //@}
private:
  //void allocTables(unsigned int fftSize);
  //void freeTables ();

  const BaseEnvironment& m_env;
  //unsigned int               m_fftSize;

  //gsl_fft_real_workspace*    m_realWkSpace;
  //gsl_fft_real_wavetable*    m_realWvTable;
  //gsl_fft_complex_workspace* m_complexWkSpace;
  //gsl_fft_complex_wavetable* m_complexWvTable;
};

}  // End namespace QUESO

#endif // UQ_FFT_H
