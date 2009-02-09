/* libs/basic/inc/uqMatrix.h
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
          void                    setPrintHorizontally(bool value) const; // Yes, 'const'
          bool                    getPrintHorizontally()           const;

  virtual unsigned int            numRows             () const = 0;
  virtual unsigned int            numCols             () const = 0;
  virtual int                     chol                () = 0;
  virtual void                    zeroLower           (bool includeDiagonal = false) = 0;
  virtual void                    zeroUpper           (bool includeDiagonal = false) = 0;
  virtual void                    print               (std::ostream& os) const = 0;

protected:
  virtual void                    copy                (const uqMatrixClass& src);

  const uqBaseEnvironmentClass& m_env;
  const Epetra_Map&             m_map;
  mutable bool                  m_printHorizontally;
};

#endif // __UQ_MATRIX_H__
