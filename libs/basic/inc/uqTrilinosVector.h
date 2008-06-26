#ifndef __UQ_TRILINOS_VECTOR_H__
#define __UQ_TRILINOS_VECTOR_H__

#include <uqVector.h>
#include <Epetra_Vector.h>

class uqTrilinosVectorClass : public uqVectorClass
{
public:
  uqTrilinosVectorClass(const Epetra_Map& map);
  uqTrilinosVectorClass(const Epetra_Map& map, double d1, double d2); // MATLAB linspace
  uqTrilinosVectorClass(const uqTrilinosVectorClass& y);
 ~uqTrilinosVectorClass();

  uqTrilinosVectorClass& operator= (const uqTrilinosVectorClass& rhs);
  uqTrilinosVectorClass& operator*=(double a);
  uqTrilinosVectorClass& operator/=(double a);
  uqTrilinosVectorClass& operator+=(const uqTrilinosVectorClass& rhs);
  uqTrilinosVectorClass& operator-=(const uqTrilinosVectorClass& rhs);

  int            rank            () const;
  unsigned int   size            () const;
  double         tranposeMultiply(const uqVectorClass& y) const;
  double         get             (unsigned int i) const;
  void           set             (unsigned int i, double value);
  void           set             (double value);
  void           invert          ();
  void           sort            ();
  void           print           (std::ostream& os) const;

  bool           atLeastOneComponentSmallerThan(const uqTrilinosVectorClass& rhs) const;
  bool           atLeastOneComponentBiggerThan (const uqTrilinosVectorClass& rhs) const;
  Epetra_Vector* data            () const;

private:
  const Epetra_Map& map          () const;

  void           copy            (const uqTrilinosVectorClass& src);
  void           scale           (double a);
  void           add             (const uqTrilinosVectorClass& y);
  void           sub             (const uqTrilinosVectorClass& y);

  Epetra_Vector* m_vec;
};

std::ostream& operator<<(std::ostream& os, const uqTrilinosVectorClass& obj);

#endif // __UQ_TRILINOS_VECTOR_H__
