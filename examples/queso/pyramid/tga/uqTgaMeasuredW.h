/* uq/examples/queso/pyramid/uqTgaMeasuredW.h
 *
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_TGA_MEASURED_W_H__
#define __UQ_TGA_MEASURED_W_H__

#include <uq1D1DFunction.h>
#include <uqTgaDefines.h>
#include <uqEnvironment.h>
#include <uqDefines.h>

template<class P_V, class P_M>
class
uqTgaMeasuredWClass
{
public:
  uqTgaMeasuredWClass(const P_V&                 params,
                      const std::vector<double>& times,
                      const std::vector<double>& temps,
                      const std::vector<double>& ws,
                      const std::vector<double>& variances);
 ~uqTgaMeasuredWClass();

  const P_V&                 params   () const;
  const std::vector<double>& times    () const;
  const std::vector<double>& temps    () const;
  const std::vector<double>& ws       () const;
  const std::vector<double>& variances() const;

protected:
  const uqBaseEnvironmentClass& m_env;
        P_V                     m_params;
        std::vector<double>     m_times;
        std::vector<double>     m_temps;
        std::vector<double>     m_ws;
        std::vector<double>     m_variances;
};

template<class P_V, class P_M>
uqTgaMeasuredWClass<P_V,P_M>::uqTgaMeasuredWClass(
  const P_V&                 params,
  const std::vector<double>& times,
  const std::vector<double>& temps,
  const std::vector<double>& ws,
  const std::vector<double>& variances)
  :
  m_env      (params.env()       ),
  m_params   (params             ),
  m_times    (times.size()    ,0.),
  m_temps    (temps.size()    ,0.),
  m_ws       (ws.size()       ,0.),
  m_variances(variances.size(),0.)
{
  if (m_env.rank() == 0) {
    std::cout << "Entering uqTgaMeasuredWClass::constructor()"
              << std::endl;
  }

  unsigned int tmpSize = m_times.size();
  for (unsigned int i = 0; i < tmpSize; ++i) {
    m_times    [i] = times    [i];
    m_temps    [i] = temps    [i];
    m_ws       [i] = ws       [i];
    m_variances[i] = variances[i];
  }

  if (m_env.rank() == 0) {
    std::cout << "Leaving uqTgaMeasuredWClass::constructor()"
              << std::endl;
  }
}

template<class P_V, class P_M>
uqTgaMeasuredWClass<P_V,P_M>::~uqTgaMeasuredWClass()
{
}

template<class P_V, class P_M>
const P_V&
uqTgaMeasuredWClass<P_V,P_M>::params() const
{
  return m_params;
}

template<class P_V, class P_M>
const std::vector<double>&
uqTgaMeasuredWClass<P_V,P_M>::times() const
{
  return m_times;
}

template<class P_V, class P_M>
const std::vector<double>&
uqTgaMeasuredWClass<P_V,P_M>::temps() const
{
  return m_temps;
}

template<class P_V, class P_M>
const std::vector<double>&
uqTgaMeasuredWClass<P_V,P_M>::ws() const
{
  return m_ws;
}

template<class P_V, class P_M>
const std::vector<double>&
uqTgaMeasuredWClass<P_V,P_M>::variances() const
{
  return m_variances;
}
#endif // __UQ_TGA_MEASURED _W_H__
