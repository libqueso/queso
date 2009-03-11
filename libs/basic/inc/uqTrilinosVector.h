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

#ifndef __UQ_TRILINOS_VECTOR_H__
#define __UQ_TRILINOS_VECTOR_H__

#include <uqVector.h>
//#include <Epetra_Vector.h>
#include <Epetra_SerialDenseMatrix.h>

class uqTrilinosVectorClass : public uqVectorClass
{
public:
  uqTrilinosVectorClass();
  uqTrilinosVectorClass(const uqBaseEnvironmentClass& env, const Epetra_Map& map);
  uqTrilinosVectorClass(const uqBaseEnvironmentClass& env, const Epetra_Map& map, double value);
  uqTrilinosVectorClass(const uqBaseEnvironmentClass& env, const Epetra_Map& map, double d1, double d2, unsigned int size); // MATLAB linspace
  uqTrilinosVectorClass(const uqTrilinosVectorClass&    v,                        double d1, double d2, unsigned int size); // MATLAB linspace
  uqTrilinosVectorClass(const uqTrilinosVectorClass&    y);
 ~uqTrilinosVectorClass();

  unsigned int numberOfProcessorsRequiredForStorage() const;

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
  void           cwSetGaussian   (const gsl_rng* rng, double mean, double stdDev);
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

  const Epetra_Map&       m_map;
  //Epetra_Vector*          m_vec;
  Epetra_SerialDenseMatrix* m_vec;
};

uqTrilinosVectorClass operator/    (      double a,                 const uqTrilinosVectorClass& x  );
uqTrilinosVectorClass operator/    (const uqTrilinosVectorClass& x, const uqTrilinosVectorClass& y  );
uqTrilinosVectorClass operator*    (      double a,                 const uqTrilinosVectorClass& x  );
uqTrilinosVectorClass operator*    (const uqTrilinosVectorClass& x, const uqTrilinosVectorClass& y  );
double                scalarProduct(const uqTrilinosVectorClass& x, const uqTrilinosVectorClass& y  );
uqTrilinosVectorClass operator+    (const uqTrilinosVectorClass& x, const uqTrilinosVectorClass& y  );
uqTrilinosVectorClass operator-    (const uqTrilinosVectorClass& x, const uqTrilinosVectorClass& y  );
std::ostream&         operator<<   (std::ostream& os,               const uqTrilinosVectorClass& obj);

#endif // __UQ_TRILINOS_VECTOR_H__
