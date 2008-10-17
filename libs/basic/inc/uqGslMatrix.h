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

#ifndef __UQ_GSL_MATRIX_H__
#define __UQ_GSL_MATRIX_H__

#include <uqMatrix.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <uqGslVector.h>

class uqGslMatrixClass : public uqMatrixClass
{
public:
  uqGslMatrixClass();
  uqGslMatrixClass(const uqEnvironmentClass& env,
                   const Epetra_Map& map,
                   unsigned int numCols);
  uqGslMatrixClass(const uqEnvironmentClass& env,
                   const Epetra_Map& map,
                   double diagValue); // MATLAB eye
  uqGslMatrixClass(const uqGslVectorClass& v,
                   double diagValue); // MATLAB eye
  uqGslMatrixClass(const uqGslVectorClass& v); // MATLAB diag
  uqGslMatrixClass(const uqGslMatrixClass& B);
 ~uqGslMatrixClass();

  uqGslMatrixClass& operator= (const uqGslMatrixClass& rhs);
  uqGslMatrixClass& operator*=(double a);
  uqGslMatrixClass& operator/=(double a);
  uqGslMatrixClass& operator+=(const uqGslMatrixClass& rhs);
  uqGslMatrixClass& operator-=(const uqGslMatrixClass& rhs);
            double& operator()(unsigned int i, unsigned int j);
      const double& operator()(unsigned int i, unsigned int j) const;

  unsigned int      numRows       () const;
  unsigned int      numCols       () const;
  int               chol          ();
  void              zeroLower     (bool includeDiagonal = false);
  void              zeroUpper     (bool includeDiagonal = false);
  uqGslMatrixClass  transpose     () const;
  uqGslVectorClass  multiply      (const uqGslVectorClass& x) const;
  uqGslVectorClass  invertMultiply(const uqGslVectorClass& b) const;
  void              invertMultiply(const uqGslVectorClass& b, uqGslVectorClass& x) const;
  void              print         (std::ostream& os) const;

private:
  void              copy          (const uqGslMatrixClass& src);
  void              multiply      (const uqGslVectorClass& x, uqGslVectorClass& y) const;

          gsl_matrix*      m_mat;
  mutable gsl_matrix*      m_LU;
  mutable gsl_permutation* m_permutation;
};

uqGslMatrixClass operator*    (double a,                    const uqGslMatrixClass& mat);
uqGslVectorClass operator*    (const uqGslMatrixClass& mat, const uqGslVectorClass& vec);
uqGslMatrixClass operator*    (const uqGslMatrixClass& m1,  const uqGslMatrixClass& m2 );
uqGslMatrixClass operator+    (const uqGslMatrixClass& m1,  const uqGslMatrixClass& m2 );
uqGslMatrixClass matrixProduct(const uqGslVectorClass& v1,  const uqGslVectorClass& v2 );
uqGslMatrixClass diagScaling  (const uqGslVectorClass& vec, const uqGslMatrixClass& mat);
std::ostream&    operator<<   (std::ostream& os,            const uqGslMatrixClass& obj);

#endif // __UQ_GSL_MATRIX_H__
