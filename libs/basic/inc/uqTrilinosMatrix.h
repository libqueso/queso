/* libs/basic/inc/uqTrilinosVector.h
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

#ifndef __UQ_TRILINOS_MATRIX_H__
#define __UQ_TRILINOS_MATRIX_H__

#include <uqMatrix.h>
#include <Epetra_SerialDenseMatrix.h>
//#include <Epetra_CrsMatrix.h>
#include <uqTrilinosVector.h>


class uqTrilinosMatrixClass : public uqMatrixClass
{
public:
  uqTrilinosMatrixClass();
  uqTrilinosMatrixClass(const uqBaseEnvironmentClass& env,
                        const Epetra_Map&         map,
                        unsigned int              numCols);

  uqTrilinosMatrixClass(const uqBaseEnvironmentClass& env,
                        const Epetra_Map&         map,
                        unsigned int              numRows,
                        unsigned int              numCols);

  uqTrilinosMatrixClass(const uqBaseEnvironmentClass& env,
                        const Epetra_Map&         map,
                        double                    diagValue); // MATLAB eye
  uqTrilinosMatrixClass(const uqTrilinosVectorClass& v,
                        double                       diagValue); // MATLAB eye
  uqTrilinosMatrixClass(const uqTrilinosVectorClass& v); // MATLAB diag
  uqTrilinosMatrixClass(const uqTrilinosMatrixClass& B);
 ~uqTrilinosMatrixClass();

  uqTrilinosMatrixClass& operator= (const uqTrilinosMatrixClass& rhs);
  uqTrilinosMatrixClass& operator*=(double a);
  uqTrilinosMatrixClass& operator/=(double a);
  uqTrilinosMatrixClass& operator+=(const uqTrilinosMatrixClass& rhs);
  uqTrilinosMatrixClass& operator-=(const uqTrilinosMatrixClass& rhs);
                 double& operator()(unsigned int i, unsigned int j);
           const double& operator()(unsigned int i, unsigned int j) const;

  unsigned int          numRows       () const;
  unsigned int          numCols       () const;
  int                   chol          ();
  void                  zeroLower     (bool includeDiagonal = false);
  void                  zeroUpper     (bool includeDiagonal = false);
  uqTrilinosMatrixClass transpose     () const;
  uqTrilinosVectorClass multiply      (const uqTrilinosVectorClass& x) const;
  uqTrilinosVectorClass invertMultiply(const uqTrilinosVectorClass& b) const;
  void                  invertMultiply(const uqTrilinosVectorClass& b, uqTrilinosVectorClass& x) const;
  void                  print         (std::ostream& os) const;

private:

  void                  copy          (const uqTrilinosMatrixClass& src);
  void                  multiply      (const uqTrilinosVectorClass& x, uqTrilinosVectorClass& y) const;

  const Epetra_Map&         m_map;
  //Epetra_CrsMatrix*       m_mat;
  Epetra_SerialDenseMatrix* m_mat;
};

uqTrilinosMatrixClass operator*    (double a,                         const uqTrilinosMatrixClass& mat);
uqTrilinosVectorClass operator*    (const uqTrilinosMatrixClass& mat, const uqTrilinosVectorClass& vec);
uqTrilinosMatrixClass operator*    (const uqTrilinosMatrixClass& m1,  const uqTrilinosMatrixClass& m2 );
uqTrilinosMatrixClass operator+    (const uqTrilinosMatrixClass& m1,  const uqTrilinosMatrixClass& m2 );
uqTrilinosMatrixClass matrixProduct(const uqTrilinosVectorClass& v1,  const uqTrilinosVectorClass& v2 );
uqTrilinosMatrixClass diagScaling  (const uqTrilinosVectorClass& vec, const uqTrilinosMatrixClass& mat);
std::ostream&         operator<<   (std::ostream& os,                 const uqTrilinosMatrixClass& obj);

#endif // __UQ_TRILINOS_MATRIX_H__
