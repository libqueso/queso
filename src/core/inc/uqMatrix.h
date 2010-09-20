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

#ifndef __UQ_MATRIX_H__
#define __UQ_MATRIX_H__

#include <uqEnvironment.h>
#include <uqVector.h>
#include <iostream>

class uqMatrixClass
{
public:
           uqMatrixClass();
           uqMatrixClass(const uqBaseEnvironmentClass& env, const Epetra_Map& map);
           uqMatrixClass(const uqMatrixClass& rhs);
  virtual ~uqMatrixClass();

  uqMatrixClass& operator= (const uqMatrixClass& rhs);
  uqMatrixClass& operator*=(double a);
  uqMatrixClass& operator+=(const uqMatrixClass& rhs);
  uqMatrixClass& operator-=(const uqMatrixClass& rhs);

    const uqBaseEnvironmentClass& env                 ()           const;
    const Epetra_Map&             map                 ()           const;
          unsigned int            numOfProcsForStorage()           const;
          void                    setPrintHorizontally(bool value) const; // Yes, 'const'
          bool                    getPrintHorizontally()           const;
          void                    setInDebugMode      (bool value) const; // Yes, 'const'
          bool                    getInDebugMode      ()           const;

  virtual unsigned int            numRowsLocal        () const = 0;
  virtual unsigned int            numRowsGlobal       () const = 0;
  virtual unsigned int            numCols             () const = 0;
  virtual int                     chol                () = 0;
  virtual void                    zeroLower           (bool includeDiagonal = false) = 0;
  virtual void                    zeroUpper           (bool includeDiagonal = false) = 0;
  virtual void                    print               (std::ostream& os) const = 0;

protected:
  virtual void                    copy                (const uqMatrixClass& src);

  const   uqBaseEnvironmentClass& m_env;
  const   Epetra_Map&             m_map;
  mutable bool                    m_printHorizontally;
  mutable bool                    m_inDebugMode;
};

#endif // __UQ_MATRIX_H__
