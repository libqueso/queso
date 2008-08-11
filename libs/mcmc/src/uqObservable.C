#include <uqObservable.h>

uqObservableClass::uqObservableClass(
  const std::string& name,
  unsigned int       numberOfObservations,
  double             priorVariance,
  double             varianceAccuracy)
  :
  m_name                (name),
  m_numberOfObservations(numberOfObservations),
  m_priorVariance       (priorVariance),
  m_varianceAccuracy    (varianceAccuracy)
{
}

uqObservableClass::~uqObservableClass()
{
}

std::string
uqObservableClass::name() const
{
  return m_name;
}

unsigned int
uqObservableClass::numberOfObservations() const
{
  return m_numberOfObservations;
}

double
uqObservableClass::priorVariance() const
{
  return m_priorVariance;
}

double
uqObservableClass::varianceAccuracy() const
{
  return m_varianceAccuracy;
}

void
uqObservableClass::setName(const std::string& name)
{
  m_name = name;
  return;
}

void
uqObservableClass::setNumberOfObservations(unsigned int numberOfObservations)
{
  m_numberOfObservations = numberOfObservations;
  return;
}

void
uqObservableClass::setPriorVariance(double priorVariance)
{
  m_priorVariance = priorVariance;
  return;
}

void
uqObservableClass::setVarianceAccuracy(double varianceAccuracy)
{
  m_varianceAccuracy = varianceAccuracy;
  return;
}

void
uqObservableClass::print(std::ostream& os) const
{
  os << this->m_name                 << " "
     << this->m_numberOfObservations << " "
     << this->m_priorVariance        << " "
     << this->m_varianceAccuracy;

  return;
}

std::ostream&
operator<<(std::ostream& os, const uqObservableClass& param)
{
  param.print(os);

  return os;
}
