#include <uqTrilinosMatrix.h>
#include <uqTrilinosVector.h>
#include <Epetra_MpiComm.h>
#include <uqDefines.h>

uqTrilinosMatrixClass::uqTrilinosMatrixClass(
  const Epetra_Map& map,
  unsigned int      numCols)
  :
  uqMatrixClass(),
  m_mat(new Epetra_CrsMatrix(Copy,map,numCols))
{
  UQ_FATAL_TEST_MACRO((m_mat == NULL),
                      m_env.rank(),
                      "uqTrilinosMatrixClass::constructor()",
                      "null matrix generated");
}
 
uqTrilinosMatrixClass::uqTrilinosMatrixClass(
  const Epetra_Map& map,
  unsigned int      dim,
  double            diagValue)
{
}

uqTrilinosMatrixClass::uqTrilinosMatrixClass(
  const uqVectorClass* v,
  bool                 invert)
{
}

uqTrilinosMatrixClass::uqTrilinosMatrixClass(
  const uqVectorClass* v1,
  const uqVectorClass* v2)
{
}
 
//uqTrilinosMatrixClass::uqTrilinosMatrixClass(
//  const Epetra_Map&    map,
//  double               condNumber,
//  const uqVectorClass* direction,
//  uqMatrixClass*       precMatrix)
//{
//}

uqTrilinosMatrixClass::uqTrilinosMatrixClass(const uqTrilinosMatrixClass& B)
{
}

uqTrilinosMatrixClass::~uqTrilinosMatrixClass()
{
  if (m_mat) delete m_mat;
}

uqTrilinosMatrixClass&
uqTrilinosMatrixClass::operator=(const uqTrilinosMatrixClass& obj)
{
  this->copy(&obj);
  return *this;
}

uqTrilinosMatrixClass&
uqTrilinosMatrixClass::operator*=(double a)
{
  this->scale(a);
  return *this;
}

uqTrilinosMatrixClass&
uqTrilinosMatrixClass::operator/=(double a)
{
  this->scale(1./a);
  return *this;
}

uqTrilinosMatrixClass&
uqTrilinosMatrixClass::operator+=(const uqTrilinosMatrixClass& rhs)
{
  this->add(&rhs);
  return *this;
}

uqTrilinosMatrixClass&
uqTrilinosMatrixClass::operator-=(const uqTrilinosMatrixClass& rhs)
{
  this->sub(&rhs);
  return *this;
}

void
uqTrilinosMatrixClass::copy(const uqTrilinosMatrixClass* src)
{
  return;
}

int
uqTrilinosMatrixClass::rank() const
{
  return this->map().Comm().MyPID();
}

unsigned int
uqTrilinosMatrixClass::numRows() const
{
  return 0;
}

unsigned int
uqTrilinosMatrixClass::numCols() const
{
  return 0;
}

double
uqTrilinosMatrixClass::get(unsigned int i, unsigned int j) const
{
  return 0.;
}

void
uqTrilinosMatrixClass::set(unsigned int i, unsigned int j, double value)
{
  return;
}

void
uqTrilinosMatrixClass::scale(double alpha)
{
  return;
}

void
uqTrilinosMatrixClass::add(const uqTrilinosMatrixClass* B)
{
  return;
}

void
uqTrilinosMatrixClass::sub(const uqTrilinosMatrixClass* B)
{
  return;
}

int
uqTrilinosMatrixClass::chol()
{
  return 0;
}

void
uqTrilinosMatrixClass::zeroLower(bool includeDiagonal)
{
  return;
}

void
uqTrilinosMatrixClass::zeroUpper(bool includeDiagonal)
{
  return;
}

uqMatrixClass*
uqTrilinosMatrixClass::multiply(const uqMatrixClass* m2) const
{
  return NULL;
}

uqMatrixClass*
uqTrilinosMatrixClass::transpose () const
{
  return NULL;
}

uqVectorClass*
uqTrilinosMatrixClass::multiply(
  const uqVectorClass* x) const
{
  UQ_FATAL_TEST_MACRO((this->numCols() != x->size()),
                      m_env.rank(),
                      "uqTrilinosMatrixClass::multiply(), vector",
                      "matrix and x have incompatible sizes");

  const uqTrilinosVectorClass* triX = dynamic_cast<const uqTrilinosVectorClass*>(x);
  UQ_FATAL_TEST_MACRO((triX == NULL),
                      m_env.rank(),
                      "uqTrilinosMatrixClass::multiply(), vector",
                      "x is not of type uqTrilinosVectorClass");

  uqTrilinosVectorClass* y = new uqTrilinosVectorClass(this->map());
  m_mat->Multiply(false,(*triX->data()),*(y->data()));

  return y;
}

const Epetra_Map&
uqTrilinosMatrixClass::map() const
{
  return (const Epetra_Map&) (m_mat->Map());
}

void
uqTrilinosMatrixClass::print(std::ostream& os) const
{
  return;
}

std::ostream&
operator<<(std::ostream& os, const uqTrilinosMatrixClass& obj)
{
  obj.print(os);

  return os;
}
