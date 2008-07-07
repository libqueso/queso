/* libs/basic/inc/uqTrilinosVector.h
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

#ifndef __UQ_TRILINOS_MATRIX_H__
#define __UQ_TRILINOS_MATRIX_H__

#include <uqMatrix.h>
#include <Epetra_CrsMatrix.h>

class uqTrilinosMatrixClass : public uqMatrixClass
{
public:
  uqTrilinosMatrixClass(const Epetra_Map& map,
                        unsigned int      numCols);
  uqTrilinosMatrixClass(const Epetra_Map& map,
                        unsigned int      dim,
                        double            diagValue); // MATLAB eye
  uqTrilinosMatrixClass(const uqVectorClass* v,
                        bool                 invert = false); // MATLAB diag or diag_inv
  uqTrilinosMatrixClass(const uqVectorClass* v1,
                        const uqVectorClass* v2); // vv2mat
  uqTrilinosMatrixClass(const uqTrilinosMatrixClass& B);
 ~uqTrilinosMatrixClass();

  uqTrilinosMatrixClass& operator= (const uqTrilinosMatrixClass& rhs);
  uqTrilinosMatrixClass& operator*=(double a);
  uqTrilinosMatrixClass& operator/=(double a);
  uqTrilinosMatrixClass& operator+=(const uqTrilinosMatrixClass& rhs);
  uqTrilinosMatrixClass& operator-=(const uqTrilinosMatrixClass& rhs);

  int            rank      () const;
  unsigned int   numRows   () const;
  unsigned int   numCols   () const;
  double         get       (unsigned int i, unsigned int j) const;
  void           set       (unsigned int i, unsigned int j, double value);
  int            chol      ();
  void           zeroLower (bool includeDiagonal = false);
  void           zeroUpper (bool includeDiagonal = false);
  uqMatrixClass* multiply  (const uqMatrixClass* m2) const;
  uqMatrixClass* transpose () const;
  uqVectorClass* multiply  (const uqVectorClass* x) const;

  void           print     (std::ostream& os) const;

private:
  const Epetra_Map& map() const;

  void           copy      (const uqTrilinosMatrixClass* src);
  void           scale     (double alpha);
  void           add       (const uqTrilinosMatrixClass* B);
  void           sub       (const uqTrilinosMatrixClass* B);

  Epetra_CrsMatrix* m_mat;
};

std::ostream& operator<<(std::ostream& os, const uqTrilinosMatrixClass& obj);

#endif // __UQ_TRILINOS_MATRIX_H__
