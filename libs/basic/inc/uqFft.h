/* uq/libs/basic/inc/uqFft.h
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

#ifndef __UQ_FFT_H__
#define __UQ_FFT_H__

#include <uqEnvironment.h>
//#include <gsl/gsl_fft_real.h>
//#include <gsl/gsl_fft_complex.h>
#include <vector>
#include <complex>

template <class T>
class uqFftClass
{
public:
  uqFftClass(const uqEnvironmentClass& env);
 ~uqFftClass();

  void forward(const std::vector<T>&                     data, 
                     unsigned int                        fftSize,
                     std::vector<std::complex<double> >& result);

  void inverse(const std::vector<T>&                     data, 
                     unsigned int                        fftSize,
                     std::vector<std::complex<double> >& result);

private:
  //void allocTables(unsigned int fftSize);
  //void freeTables ();

  const uqEnvironmentClass& m_env;
  //unsigned int               m_fftSize;

  //gsl_fft_real_workspace*    m_realWkSpace;
  //gsl_fft_real_wavetable*    m_realWvTable;
  //gsl_fft_complex_workspace* m_complexWkSpace;
  //gsl_fft_complex_wavetable* m_complexWvTable;
};

template <class T>
uqFftClass<T>::uqFftClass(const uqEnvironmentClass& env)
  :
  m_env(env)
{
}

template <class T>
uqFftClass<T>::~uqFftClass()
{
}

#endif // __UQ_FFT_H__
