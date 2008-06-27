#include <uqMatrix.h>
#include <uqDefines.h>

uqMatrixClass::uqMatrixClass()
  :
  m_env(*(new uqEnvironmentClass()))
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqMatrixClass::constructor(), default",
                      "should not be used by user");
}

uqMatrixClass::uqMatrixClass(const uqEnvironmentClass& env)
  :
  m_env(env)
{
}

uqMatrixClass::uqMatrixClass(const uqMatrixClass& rhs)
  :
  m_env(rhs.m_env)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.rank(),
                      "uqMatrixClass::constructor(), copy",
                      "code should not execute through here");
}

uqMatrixClass::~uqMatrixClass()
{
}

void
uqMatrixClass::copy(const uqMatrixClass& src)
{
  //m_env = env;

  return;
}

uqMatrixClass&
uqMatrixClass::operator= (const uqMatrixClass& rhs)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.rank(),
                      "uqMatrixClass::operator=()",
                      "code should not execute through here");
  return *this;
}

uqMatrixClass&
uqMatrixClass::operator*=(double a)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.rank(),
                      "uqMatrixClass::operator*=()",
                      "code should not execute through here");
  return *this;
}

uqMatrixClass&
uqMatrixClass::operator+=(const uqMatrixClass& rhs)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.rank(),
                      "uqMatrixClass::operator+=()",
                      "code should not execute through here");
  return *this;
}

uqMatrixClass&
uqMatrixClass::operator-=(const uqMatrixClass& rhs)
{
  UQ_FATAL_TEST_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                      m_env.rank(),
                      "uqMatrixClass::operator-=()",
                      "code should not execute through here");
  return *this;
}

const uqEnvironmentClass&
uqMatrixClass::env() const
{
  return m_env;
}
