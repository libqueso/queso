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

#ifndef __UQ_VECTOR_H__
#define __UQ_VECTOR_H__

#include <uqEnvironment.h>
#include <Epetra_Map.h>
#include <gsl/gsl_randist.h>
#include <iostream>

class uqVectorClass
{
public:
           uqVectorClass();
           uqVectorClass(const uqBaseEnvironmentClass& env, const Epetra_Map& map);
           uqVectorClass(const uqVectorClass& rhs);
  virtual ~uqVectorClass();

  uqVectorClass& operator= (const uqVectorClass& rhs);
  uqVectorClass& operator*=(double a);
  uqVectorClass& operator/=(double a);
  uqVectorClass& operator+=(const uqVectorClass& rhs);
  uqVectorClass& operator-=(const uqVectorClass& rhs);

    const uqBaseEnvironmentClass& env                 ()           const;
    const Epetra_Map&             map                 ()           const;
          void                    setPrintHorizontally(bool value) const; // Yes, 'const'
          bool                    getPrintHorizontally()           const;

  virtual unsigned int            size                () const = 0;
  virtual void                    cwSet               (double value) = 0;
  virtual void                    cwSetGaussian       (const gsl_rng* rng, double mean, double stdDev) = 0;
  virtual void                    cwInvert            () = 0;
  virtual void                    sort                () = 0;
  virtual void                    print               (std::ostream& os) const = 0;

protected:
  virtual void                    copy                (const uqVectorClass& src);

  const uqBaseEnvironmentClass& m_env;
  const Epetra_Map&             m_map;
  mutable bool                  m_printHorizontally;
};

#endif // __UQ_VECTOR_H__
