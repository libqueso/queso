/* uq/examples/queso/pyramid/uqTgaTemperatureProfiles.h
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

#ifndef __UQ_TGA_TEMPERATURE_PROFILES_H__
#define __UQ_TGA_TEMPERATURE_PROFILES_H__

double temperatureFunction1_Value(double time, const void* info)
{
  double value = 0.;

#if 0
  double a1 = 1.;
  double T1 = 950.;
  double a2 = -1.;
  double b2 = 60.;
  double a3 = -1./54.;
  double b3 = 11./6.;
  double c3 = -60;
#endif
#if 0
  double a1 = 2.;
  double T1 = 1850.;
  double a2 = -2.;
  double b2 = 120.;
  double a3 = -1./54.;
  double b3 = 17./6.;
  double c3 = -120;
#endif
#if 1
  double a1 = 5.;
  double T1 = 4550.;
  double a2 = -5.;
  double b2 = 300.;
  double a3 = -1./54.;
  double b3 = 35./6.;
  double c3 = -300;
#endif

  if (time <= 0.) {
    value = 50.;
  }
  else if (( 0. < time) && (time <= 30. )) {
    value = 50. + a1*time*time;
  }
  else if ((30. < time) && (time <= 90. )) {
    value = T1 + a2*(time-30.)*(time-30.) + b2*(time-30.);
  }
  else if ((90. < time) && (time <= 120.)) {
    value = T1 + a3*(time-90.)*(time-90.)*(time-90.) + b3*(time-90.)*(time-90.) + c3*(time-90.);
  }
  else {
    value = 300.;
  }

  return value;
}

double temperatureFunction1_Deriv(double time, const void* info)
{
  double deriv = 0.;

#if 0
  double a1 = 1.;
  double a2 = -1.;
  double b2 = 60.;
  double a3 = -1./54.;
  double b3 = 11./6.;
  double c3 = -60;
#endif
#if 0
  double a1 = 2.;
  double a2 = -4;
  double b2 = 120.;
  double a3 = -1./54.;
  double b3 = 17./6.;
  double c3 = -120;
#endif
#if 1
  double a1 = 5.;
  double a2 = -5;
  double b2 = 300.;
  double a3 = -1./54.;
  double b3 = 35./6.;
  double c3 = -300;
#endif

  if (time <= 0.) {
    deriv = 0.;
  }
  else if (( 0. < time) && (time <= 30. )) {
    deriv = 2.*a1*time;
  }
  else if ((30. < time) && (time <= 90. )) {
    deriv = 2.*a2*(time-30.) + b2;
  }
  else if ((90. < time) && (time <= 120.)) {
    deriv = 3.*a3*(time-90.)*(time-90.) + 2.*b3*(time-90.) + c3;
  }
  else {
    deriv = 0.;
  }

  return deriv;
}

double temperatureFunction2_Value(double time, const void* info)
{
  double value = 0.;

  double a = 3.;
  double T0 = 288.;
  double T1 = 4140.;
  double T2 = 288.;
  double t1 = (T1-T0)/a;
  double t2 = t1+300.;
  double t3 = t2+(T1-T2)/a;

  if (time <= 0.) {
    value = T0;
  }
  else if (( 0. < time) && (time <= t1 )) {
    value = T0 + a*time;
  }
  else if (( t1 < time) && (time <= t2 )) {
    value = T1;
  }
  else if (( t2 < time) && (time <= t3 )) {
    value = T1 - a*(time-t2);
  }
  else {
    value = T2;
  }

  return value;
}

double temperatureFunction2_Deriv(double time, const void* info)
{
  double deriv = 0.;

  double a = 3.;
  double T0 = 288.;
  double T1 = 4140.;
  double T2 = 288.;
  double t1 = (T1-T0)/a;
  double t2 = t1+300.;
  double t3 = t2+(T1-T2)/a;

  if (time <= 0.) {
    deriv = 0.;
  }
  else if (( 0. < time) && (time <= t1 )) {
    deriv = a;
  }
  else if (( t1 < time) && (time <= t2 )) {
    deriv = 0.;
  }
  else if (( t2 < time) && (time <= t3 )) {
    deriv = -a;
  }
  else {
    deriv = 0.;
  }

  return deriv;
}

double temperatureFunction3_Value(double time, const void* info)
{
  double value = 0.;

  double a = 6.;
  double T0 = 288.;
  double T1 = 5080;
  double T2 = 288.;
  double t1 = (T1-T0)/a;
  double t2 = t1+300.;
  double t3 = t2+(T1-T2)/a;

  if (time <= 0.) {
    value = T0;
  }
  else if (( 0. < time) && (time <= t1 )) {
    value = T0 + a*time;
  }
  else if (( t1 < time) && (time <= t2 )) {
    value = T1;
  }
  else if (( t2 < time) && (time <= t3 )) {
    value = T1 - a*(time-t2);
  }
  else {
    value = T2;
  }

  return value;
}

double temperatureFunction3_Deriv(double time, const void* info)
{
  double deriv = 0.;

  double a = 6.;
  double T0 = 288.;
  double T1 = 5080;
  double T2 = 288.;
  double t1 = (T1-T0)/a;
  double t2 = t1+300.;
  double t3 = t2+(T1-T2)/a;

  if (time <= 0.) {
    deriv = 0.;
  }
  else if (( 0. < time) && (time <= t1 )) {
    deriv = a;
  }
  else if (( t1 < time) && (time <= t2 )) {
    deriv = 0.;
  }
  else if (( t2 < time) && (time <= t3 )) {
    deriv = -a;
  }
  else {
    deriv = 0.;
  }

  return deriv;
}

double temperatureFunction4_Value(double time, const void* info)
{
  double value = 0.;

  double a = 24.;
  double T0 = 288.;
  double T1 = 7200.;
  double T2 = 288.;
  double t1 = (T1-T0)/a;
  double t2 = t1+300.;
  double t3 = t2+(T1-T2)/a;

  if (time <= 0.) {
    value = T0;
  }
  else if (( 0. < time) && (time <= t1 )) {
    value = T0 + a*time;
  }
  else if (( t1 < time) && (time <= t2 )) {
    value = T1;
  }
  else if (( t2 < time) && (time <= t3 )) {
    value = T1 - a*(time-t2);
  }
  else {
    value = T2;
  }

  return value;
}

double temperatureFunction4_Deriv(double time, const void* info)
{
  double deriv = 0.;

  double a = 24.;
  double T0 = 288.;
  double T1 = 7200.;
  double T2 = 288.;
  double t1 = (T1-T0)/a;
  double t2 = t1+300.;
  double t3 = t2+(T1-T2)/a;

  if (time <= 0.) {
    deriv = 0.;
  }
  else if (( 0. < time) && (time <= t1 )) {
    deriv = a;
  }
  else if (( t1 < time) && (time <= t2 )) {
    deriv = 0.;
  }
  else if (( t2 < time) && (time <= t3 )) {
    deriv = -a;
  }
  else {
    deriv = 0.;
  }

  return deriv;
}
#endif // __UQ_TGA_TEMPERATURE_PROFILES_H__
