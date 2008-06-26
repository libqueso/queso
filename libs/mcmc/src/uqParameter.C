#include <uqParameter.h>

uqParameterClass::uqParameterClass(
  const std::string& name,
  double             initialValue,
  double             minValue,
  double             maxValue,
  double             priorMu,
  double             priorSigma)
  :
  m_name        (name),
  m_initialValue(initialValue),
  m_minValue    (minValue),
  m_maxValue    (maxValue),
  m_priorMu     (priorMu),
  m_priorSigma  (priorSigma)
{
}

uqParameterClass::~uqParameterClass()
{
}

std::string
uqParameterClass::name() const
{
  return m_name;
}

double
uqParameterClass::initialValue() const
{
  return m_initialValue;
}

double
uqParameterClass::minValue() const
{
  return m_minValue;
}

double
uqParameterClass::maxValue() const
{
  return m_maxValue;
}

double
uqParameterClass::priorMu() const
{
  return m_priorMu;
}

double
uqParameterClass::priorSigma() const
{
  return m_priorSigma;
}

void
uqParameterClass::setName(const std::string& name)
{
  m_name = name;
  return;
}

void
uqParameterClass::setInitialValue(double initialValue)
{
  m_initialValue = initialValue;
  return;
}

void
uqParameterClass::setMinValue(double minValue)
{
  m_minValue = minValue;
  return;
}

void
uqParameterClass::setMaxValue(double maxValue)
{
  m_maxValue = maxValue;
  return;
}

void
uqParameterClass::setPriorMu(double priorMu)
{
  m_priorMu = priorMu;
  return;
}

void
uqParameterClass::setPriorSigma(double priorSigma)
{
  m_priorSigma = priorSigma;
  return;
}

void
uqParameterClass::print(std::ostream& os) const
{
  os << this->m_name         << " "
     << this->m_initialValue << " "
     << this->m_minValue     << " "
     << this->m_maxValue     << " "
     << this->m_priorMu      << " "
     << this->m_priorSigma;

  return;
}

std::ostream&
operator<<(std::ostream& os, const uqParameterClass& param)
{
  param.print(os);

  return os;
}
