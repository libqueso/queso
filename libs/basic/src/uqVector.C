#include <uqVector.h>
#include <uqDefines.h>

uqVectorClass::uqVectorClass()
  :
  m_env(*(new uqEnvironmentClass()))
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqVectorClass::constructor(), default",
                      "should not be used by user");
}

uqVectorClass::uqVectorClass(const uqEnvironmentClass& env)
  :
  m_env(env),
  m_printHorizontally(true)
{
}

uqVectorClass::uqVectorClass(const uqVectorClass& rhs)
  :
  m_env(rhs.m_env)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.rank(),
                      "uqVectorClass::constructor(), copy",
                      "code should not execute through here");
}

uqVectorClass::~uqVectorClass()
{
}

void
uqVectorClass::copy(const uqVectorClass& src)
{
  //m_env               = env;
  m_printHorizontally = src.m_printHorizontally;

  return;
}

uqVectorClass&
uqVectorClass::operator= (const uqVectorClass& rhs)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.rank(),
                      "uqVectorClass::operator=()",
                      "code should not execute through here");
  return *this;
}

uqVectorClass&
uqVectorClass::operator*=(double a)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.rank(),
                      "uqVectorClass::operator*=()",
                      "code should not execute through here");
  return *this;
}

uqVectorClass&
uqVectorClass::operator/=(double a)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.rank(),
                      "uqVectorClass::operator/=()",
                      "code should not execute through here");
  return *this;
}

uqVectorClass&
uqVectorClass::operator+=(const uqVectorClass& rhs)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.rank(),
                      "uqVectorClass::operator+=()",
                      "code should not execute through here");
  return *this;
}

uqVectorClass&
uqVectorClass::operator-=(const uqVectorClass& rhs)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.rank(),
                      "uqVectorClass::operator-=()",
                      "code should not execute through here");
  return *this;
}

const uqEnvironmentClass&
uqVectorClass::env() const
{
  return m_env;
}

void
uqVectorClass::setPrintHorizontally(bool value) const
{
  m_printHorizontally = value;
  return;
}

bool
uqVectorClass::getPrintHorizontally() const
{
  return m_printHorizontally;
}
