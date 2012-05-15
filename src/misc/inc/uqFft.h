//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012 The PECOS Development Team
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
// 
// $Id$
//
//--------------------------------------------------------------------------

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
  uqFftClass(const uqBaseEnvironmentClass& env);
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

  const uqBaseEnvironmentClass& m_env;
  //unsigned int               m_fftSize;

  //gsl_fft_real_workspace*    m_realWkSpace;
  //gsl_fft_real_wavetable*    m_realWvTable;
  //gsl_fft_complex_workspace* m_complexWkSpace;
  //gsl_fft_complex_wavetable* m_complexWvTable;
};

template <class T>
uqFftClass<T>::uqFftClass(const uqBaseEnvironmentClass& env)
  :
  m_env(env)
{
}

template <class T>
uqFftClass<T>::~uqFftClass()
{
}

#endif // __UQ_FFT_H__
