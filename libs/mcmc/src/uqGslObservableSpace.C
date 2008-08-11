#include <uqObservableSpace.h>
#include <uqGslMatrix.h>

template<>
void
uqObservableSpaceClass<uqGslVectorClass,uqGslMatrixClass>::createNumbersOfObservations() const
{
  m_numbersOfObservations = new uqGslVectorClass(m_env,this->dim());

  for (unsigned int i = 0; i < m_observables.size(); ++i) {
    if (m_observables[i]) (*m_numbersOfObservations)[i] = (double) m_observables[i]->numberOfObservations();
  }

  return;
}

template<>
void
uqObservableSpaceClass<uqGslVectorClass,uqGslMatrixClass>::createPriorVariances() const
{
  m_priorVariances = new uqGslVectorClass(m_env,this->dim());

  for (unsigned int i = 0; i < m_observables.size(); ++i) {
    if (m_observables[i]) (*m_priorVariances)[i] = m_observables[i]->priorVariance();
  }

  return; 
}

template<>
void
uqObservableSpaceClass<uqGslVectorClass,uqGslMatrixClass>::createVarianceAccuracies() const
{
  m_varianceAccuracies = new uqGslVectorClass(m_env,this->dim());

  for (unsigned int i = 0; i < m_observables.size(); ++i) {
    if (m_observables[i]) (*m_varianceAccuracies)[i] = m_observables[i]->varianceAccuracy();
  }

  return;
}

template<>
void
uqObservableSpaceClass<uqGslVectorClass,uqGslMatrixClass>::printObservableNames(std::ostream& os, bool printHorizontally) const
{
  if (printHorizontally) { 
    for (unsigned int i = 0; i < this->dim(); ++i) {
      os << "'" << this->observable(i).name() << "'"
         << " ";
    }
  }
  else {
    for (unsigned int i = 0; i < this->dim(); ++i) {
      os << "'" << this->observable(i).name() << "'"
         << std::endl;
    }
  }

  return;
}
