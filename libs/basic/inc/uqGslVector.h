/* libs/basic/inc/uqGslVector.h
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

#ifndef __UQ_GSL_VECTOR_H__
#define __UQ_GSL_VECTOR_H__

#include <uqVector.h>
#include <gsl/gsl_vector.h>

class uqGslVectorClass : public uqVectorClass
{
public:
  uqGslVectorClass();
  uqGslVectorClass(const uqEnvironmentClass& env, const Epetra_Map& map);
  uqGslVectorClass(const uqEnvironmentClass& env, const Epetra_Map& map, double value);
  uqGslVectorClass(const uqEnvironmentClass& env, double d1, double d2, const Epetra_Map& map); // MATLAB linspace
  uqGslVectorClass(const uqGslVectorClass&     v, double d1, double d2);                        // MATLAB linspace
  uqGslVectorClass(const uqGslVectorClass& y);
 ~uqGslVectorClass();

  uqGslVectorClass& operator= (const uqGslVectorClass& rhs);
  uqGslVectorClass& operator*=(double a);
  uqGslVectorClass& operator/=(double a);
  uqGslVectorClass& operator*=(const uqGslVectorClass& rhs);
  uqGslVectorClass& operator/=(const uqGslVectorClass& rhs);
  uqGslVectorClass& operator+=(const uqGslVectorClass& rhs);
  uqGslVectorClass& operator-=(const uqGslVectorClass& rhs);
            double& operator[](unsigned int i);
      const double& operator[](unsigned int i) const;

  unsigned int size             () const;
  double       norm2Sq          () const;
  double       norm2            () const;
  double       sumOfComponents  () const;
  void         cwSet            (double value);
  void         cwSetGaussian    (const gsl_rng* rng, double mean, double stdDev);
  void         cwSetGaussian    (const gsl_rng* rng, const uqGslVectorClass& meanVec, const uqGslVectorClass& stdDevVec);
  void         cwSetUniform     (const gsl_rng* rng, const uqGslVectorClass& aVec,    const uqGslVectorClass& bVec     );
  void         cwInvert         ();
  void         sort             ();
  void         print            (std::ostream& os) const;

  bool         atLeastOneComponentSmallerThan(const uqGslVectorClass& rhs) const;
  bool         atLeastOneComponentBiggerThan (const uqGslVectorClass& rhs) const;
  gsl_vector*  data             () const; // Necessary for uqGslMatrixClass::invertMultiply()

private:

  void         copy             (const uqGslVectorClass& src);

  gsl_vector* m_vec;
};

uqGslVectorClass operator/    (const double a,            const uqGslVectorClass& x  );
uqGslVectorClass operator/    (const uqGslVectorClass& x, const uqGslVectorClass& y  );
uqGslVectorClass operator*    (const double a,            const uqGslVectorClass& x  );
uqGslVectorClass operator*    (const uqGslVectorClass& x, const uqGslVectorClass& y  );
double           scalarProduct(const uqGslVectorClass& x, const uqGslVectorClass& y  );
uqGslVectorClass operator+    (const uqGslVectorClass& x, const uqGslVectorClass& y  );
uqGslVectorClass operator-    (const uqGslVectorClass& x, const uqGslVectorClass& y  );
std::ostream&    operator<<   (std::ostream& os,          const uqGslVectorClass& obj);

#endif // __UQ_GSL_VECTOR_H__
