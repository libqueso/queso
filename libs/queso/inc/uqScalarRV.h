/* uq/libs/queso/inc/uqScalarRV.h
 *
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu/
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

#ifndef __UQ_SCALAR_RV_H__
#define __UQ_SCALAR_RV_H__

#include <string>
#include <iostream>
#include <math.h>

class uqBaseScalarRVClass
{
public:
  uqBaseScalarRVClass(double minValue    = -INFINITY,
                      double maxValue    = INFINITY,
                      double expectValue = 0.,
                      double stdDevValue = INFINITY);
 ~uqBaseScalarRVClass();

  double minValue      () const;
  double maxValue      () const;
  double expectValue   () const;
  double stdDevValue   () const;

  void   setMinValue   (double minValue);
  void   setMaxValue   (double maxValue);
  void   setExpectValue(double expect);
  void   setStdDevValue(double stdDev);

  void   print         (std::ostream& os) const;

private:
  double m_minValue;
  double m_maxValue;
  double m_expectValue;
  double m_stdDevValue;
};

std::ostream& operator<<(std::ostream& os, const uqBaseScalarRVClass& obj);

#endif // __UQ_SCALAR_RV_H__
