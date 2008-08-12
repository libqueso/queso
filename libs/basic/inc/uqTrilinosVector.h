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

#ifndef __UQ_TRILINOS_VECTOR_H__
#define __UQ_TRILINOS_VECTOR_H__

#include <uqVector.h>
//#include <Epetra_Vector.h>
#include <Epetra_SerialDenseMatrix.h>

class uqTrilinosVectorClass : public uqVectorClass
{
public:
  uqTrilinosVectorClass();
  uqTrilinosVectorClass(const uqEnvironmentClass& env, const Epetra_Map& map);
  uqTrilinosVectorClass(const uqEnvironmentClass& env, const Epetra_Map& map, double d1, double d2, unsigned int size); // MATLAB linspace
  uqTrilinosVectorClass(const uqTrilinosVectorClass& v,                       double d1, double d2, unsigned int size); // MATLAB linspace
  uqTrilinosVectorClass(const uqTrilinosVectorClass& y);
 ~uqTrilinosVectorClass();

  uqTrilinosVectorClass& operator= (const uqTrilinosVectorClass& rhs);
  uqTrilinosVectorClass& operator*=(double a);
  uqTrilinosVectorClass& operator/=(double a);
  uqTrilinosVectorClass& operator*=(const uqTrilinosVectorClass& rhs);
  uqTrilinosVectorClass& operator/=(const uqTrilinosVectorClass& rhs);
  uqTrilinosVectorClass& operator+=(const uqTrilinosVectorClass& rhs);
  uqTrilinosVectorClass& operator-=(const uqTrilinosVectorClass& rhs);
                 double& operator[](unsigned int i);
           const double& operator[](unsigned int i) const;

  unsigned int   size            () const;
  double         norm2Sq         () const;
  double         norm2           () const;
  double         sumOfComponents () const;
  void           cwSet           (double value);
  void           cwSetGaussian   (gsl_rng* rng, double mean, double stdDev);
  void           cwInvert        ();
  void           sort            ();
  void           print           (std::ostream& os) const;

  bool           atLeastOneComponentSmallerThan(const uqTrilinosVectorClass& rhs) const;
  bool           atLeastOneComponentBiggerThan (const uqTrilinosVectorClass& rhs) const;
  //Epetra_Vector* data          () const;
  Epetra_SerialDenseMatrix* data () const;
  const Epetra_Map& map          () const;

private:

  void           copy            (const uqTrilinosVectorClass& src);

  const Epetra_Map&         m_map;
  //Epetra_Vector*          m_vec;
  Epetra_SerialDenseMatrix* m_vec;
};

uqTrilinosVectorClass operator/    (const double a,                 const uqTrilinosVectorClass& x  );
uqTrilinosVectorClass operator/    (const uqTrilinosVectorClass& x, const uqTrilinosVectorClass& y  );
uqTrilinosVectorClass operator*    (const double a,                 const uqTrilinosVectorClass& x  );
uqTrilinosVectorClass operator*    (const uqTrilinosVectorClass& x, const uqTrilinosVectorClass& y  );
double                scalarProduct(const uqTrilinosVectorClass& x, const uqTrilinosVectorClass& y  );
uqTrilinosVectorClass operator+    (const uqTrilinosVectorClass& x, const uqTrilinosVectorClass& y  );
uqTrilinosVectorClass operator-    (const uqTrilinosVectorClass& x, const uqTrilinosVectorClass& y  );
std::ostream&         operator<<   (std::ostream& os,               const uqTrilinosVectorClass& obj);

#endif // __UQ_TRILINOS_VECTOR_H__
