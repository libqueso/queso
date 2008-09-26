/* uq/libs/mcmc/inc/uqBasicScalarRV.h
 *
 * Copyright (C) 2008 The PECOS Team, http://www.ices.utexas.edu/centers/pecos
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

#ifndef __UQ_BASIC_SCALAR_RV_H__
#define __UQ_BASIC_SCALAR_RV_H__

#include <string>
#include <iostream>
#include <math.h>

class uqBasicScalarRVClass
{
public:
  uqBasicScalarRVClass(const std::string& name,
                       double             initialValue,
                       double             minValue = -INFINITY,
                       double             maxValue = INFINITY,
                       double             priorMu = 0.,
                       double             priorSigma = INFINITY);
 ~uqBasicScalarRVClass();

  std::string name           () const;
  double      initialValue   () const;
  double      minValue       () const;
  double      maxValue       () const;
  double      priorMu        () const;
  double      priorSigma     () const;

  void        setName        (const std::string& name);
  void        setInitialValue(double initialValue);
  void        setMinValue    (double minValue);
  void        setMaxValue    (double maxValue);
  void        setPriorMu     (double priorMu);
  void        setPriorSigma  (double priorSigma);

  void        print          (std::ostream& os) const;

private:
  std::string m_name;
  double      m_initialValue;
  double      m_minValue;
  double      m_maxValue;
  double      m_priorMu;
  double      m_priorSigma;
};

std::ostream& operator<<(std::ostream& os, const uqBasicScalarRVClass& obj);

#endif // __UQ_BASIC_SCALAR_RV_H__
