#include <uq1DNode.h>

uq1DNodeClass::uq1DNodeClass(
  unsigned int       globalId,
  double             x,
  uqNodePositionEnum nodePositionType,
  uqBCEnum           bcType,
  double             bcValue)
  :
  m_myGlobalId    (globalId),
  m_myX           (x),
  m_myPositionType(nodePositionType),
  m_myBCType      (bcType),
  m_myBCValue     (bcValue)
{
}

uq1DNodeClass::uq1DNodeClass(const uq1DNodeClass& obj)
{
  this->copy(obj);
}

uq1DNodeClass&
uq1DNodeClass::operator=(const uq1DNodeClass& rhs)
{
  this->copy(rhs);

  return *this;
}

void
uq1DNodeClass::copy(const uq1DNodeClass& src)
{
  m_myGlobalId     = src.m_myGlobalId;
  m_myX            = src.m_myX;
  m_myPositionType = src.m_myPositionType;
  m_myBCType       = src.m_myBCType;
  m_myBCValue      = src.m_myBCValue;

  return;
}

uq1DNodeClass::~uq1DNodeClass()
{
}

unsigned int
uq1DNodeClass::globalId() const
{
  return m_myGlobalId;
}

double
uq1DNodeClass::x() const
{
  return m_myX;
}

void
uq1DNodeClass::print(std::ostream& os) const
{
  return;
}
