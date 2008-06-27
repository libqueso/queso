#include <uqTrilinosVector.h>
#include <Epetra_MpiComm.h>
#include <uqDefines.h>

uqTrilinosVectorClass::uqTrilinosVectorClass(const Epetra_Map& map)
  :
  uqVectorClass(),
  m_vec(new Epetra_Vector(map,true))
{
  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.rank(),
                      "uqTrilinosVectorClass::constructor()",
                      "null vector generated");
}

uqTrilinosVectorClass::uqTrilinosVectorClass(const Epetra_Map& map, double d1, double d2)
  :
  uqVectorClass(),
  m_vec(new Epetra_Vector(map,true))
{
  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.rank(),
                      "uqTrilinosVectorClass::constructor(), linspace",
                      "null vector generated");

  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.rank(),
                    "uqTrilinosVectorClass::constructor(), linspace",
                    "failed");
}

uqTrilinosVectorClass::uqTrilinosVectorClass(const uqTrilinosVectorClass& v)
  :
  uqVectorClass(),
  m_vec(new Epetra_Vector(v.map(),true))
{
  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.rank(),
                      "uqTrilinosVectorClass::constructor(), copy",
                      "null vector generated");
  this->uqVectorClass::copy(v);
  this->copy(v);
}

uqTrilinosVectorClass::~uqTrilinosVectorClass()
{
  if (m_vec) delete m_vec;
}

uqTrilinosVectorClass&
uqTrilinosVectorClass::operator=(const uqTrilinosVectorClass& obj)
{
  this->copy(obj);
  return *this;
}

uqTrilinosVectorClass&
uqTrilinosVectorClass::operator*=(double a)
{
  this->scale(a);
  return *this;
}

uqTrilinosVectorClass&
uqTrilinosVectorClass::operator/=(double a)
{
  this->scale(1./a);
  return *this;
}

uqTrilinosVectorClass&
uqTrilinosVectorClass::operator+=(const uqTrilinosVectorClass& rhs)
{
  this->add(rhs);
  return *this;
}

uqTrilinosVectorClass&
uqTrilinosVectorClass::operator-=(const uqTrilinosVectorClass& rhs)
{
  this->sub(rhs);
  return *this;
}

void
uqTrilinosVectorClass::copy(const uqTrilinosVectorClass& src)
{
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.rank(),
                    "uqTrilinosVectorClass::copy()",
                    "failed");

  return;
}

int
uqTrilinosVectorClass::rank() const
{
  return m_vec->Map().Comm().MyPID();
}

unsigned int
uqTrilinosVectorClass::size() const
{
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.rank(),
                    "uqTrilinosVectorClass::size()",
                    "failed");

  return 0;
}

double
uqTrilinosVectorClass::tranposeMultiply(const uqVectorClass& y) const
{
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.rank(),
                    "uqTrilinosVectorClass::tranposeMultiply()",
                    "failed");
  return 0.;
}

double
uqTrilinosVectorClass::get(unsigned int i) const
{
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.rank(),
                    "uqTrilinosVectorClass::get()",
                    "failed");
  return 0.;
}

void
uqTrilinosVectorClass::set(unsigned int i, double value)
{
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.rank(),
                    "uqTrilinosVectorClass::set()",
                    "failed");
  return;
}

void
uqTrilinosVectorClass::set(double value)
{
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.rank(),
                    "uqTrilinosVectorClass::set()",
                    "failed");
  return;
}

void
uqTrilinosVectorClass::scale(double a)
{
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.rank(),
                    "uqTrilinosVectorClass::scale()",
                    "failed");
  return;
}

void
uqTrilinosVectorClass::add(const uqTrilinosVectorClass& y)
{
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.rank(),
                    "uqTrilinosVectorClass::add()",
                    "failed");
  return;
}

void
uqTrilinosVectorClass::sub(const uqTrilinosVectorClass& y)
{
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.rank(),
                    "uqTrilinosVectorClass::sub()",
                    "failed");
  return;
}

void
uqTrilinosVectorClass::invert()
{
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.rank(),
                    "uqTrilinosVectorClass::invert()",
                    "failed");
  return;
}

void
uqTrilinosVectorClass::sort()
{
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.rank(),
                    "uqTrilinosVectorClass::sort()",
                    "failed");
  return;
}

void
uqTrilinosVectorClass::print(std::ostream& os) const
{
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.rank(),
                    "uqTrilinosVectorClass::print()",
                    "failed");
  return;
}

const Epetra_Map&
uqTrilinosVectorClass::map() const
{
  return (const Epetra_Map&) (m_vec->Map());
}

Epetra_Vector*
uqTrilinosVectorClass::data() const
{
  return m_vec;
}

bool
uqTrilinosVectorClass::atLeastOneComponentSmallerThan(const uqTrilinosVectorClass& rhs) const
{
  UQ_FATAL_TEST_MACRO((this->size() != rhs.size()),
                      m_env.rank(),
                      "uqTrilinosVectorClass::atLeastOneComponentSmallerThan()",
                      "vectors have different sizes");

  bool result = false;

  return result;
}

bool
uqTrilinosVectorClass::atLeastOneComponentBiggerThan (const uqTrilinosVectorClass& rhs) const
{
  UQ_FATAL_TEST_MACRO((this->size() != rhs.size()),
                      m_env.rank(),
                      "uqTrilinosVectorClass::atLeastOneComponentBiggerThan()",
                      "vectors have different sizes");

  bool result = false;

  return result;
}

std::ostream&
operator<<(std::ostream& os, const uqTrilinosVectorClass& obj)
{
  obj.print(os);

  return os;
}
