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
