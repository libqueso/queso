#include <uqParamSpace.h>
#include <uqGslMatrix.h>

template<>
void
uqParamSpaceClass<uqGslVectorClass,uqGslMatrixClass>::createInitialValues() const
{
  m_initialValues = new uqGslVectorClass(m_env,this->dim());

  for (unsigned int i = 0; i < m_parameters.size(); ++i) {
    if (m_parameters[i]) (*m_initialValues)[i] = m_parameters[i]->initialValue();
  }

  return;
}

template<>
void
uqParamSpaceClass<uqGslVectorClass,uqGslMatrixClass>::createMinValues() const
{
  m_minValues = new uqGslVectorClass(m_env,this->dim());

  for (unsigned int i = 0; i < m_parameters.size(); ++i) {
    if (m_parameters[i]) (*m_minValues)[i] = m_parameters[i]->minValue();
  }

  return; 
}

template<>
void
uqParamSpaceClass<uqGslVectorClass,uqGslMatrixClass>::createMaxValues() const
{
  m_maxValues = new uqGslVectorClass(m_env,this->dim());

  for (unsigned int i = 0; i < m_parameters.size(); ++i) {
    if (m_parameters[i]) (*m_maxValues)[i] = m_parameters[i]->maxValue();
  }

  return;
}

template<>
void
uqParamSpaceClass<uqGslVectorClass,uqGslMatrixClass>::createPriorMuValues() const
{
  m_priorMuValues = new uqGslVectorClass(m_env,this->dim());

  for (unsigned int i = 0; i < m_parameters.size(); ++i) {
    if (m_parameters[i]) (*m_priorMuValues)[i] = m_parameters[i]->priorMu();
  }

  return; 
}

template<>
void
uqParamSpaceClass<uqGslVectorClass,uqGslMatrixClass>::createPriorSigmaValues() const
{
  m_priorSigmaValues = new uqGslVectorClass(m_env,this->dim());

  for (unsigned int i = 0; i < m_parameters.size(); ++i) {
    if (m_parameters[i]) (*m_priorSigmaValues)[i] = m_parameters[i]->priorSigma();
  }

  return;
}

template<>
void
uqParamSpaceClass<uqGslVectorClass,uqGslMatrixClass>::printParameterNames(std::ostream& os, bool printHorizontally) const
{
  if (printHorizontally) { 
    for (unsigned int i = 0; i < this->dim(); ++i) {
      os << "'" << this->parameter(i).name() << "'"
         << " ";
    }
  }
  else {
    for (unsigned int i = 0; i < this->dim(); ++i) {
      os << "'" << this->parameter(i).name() << "'"
         << std::endl;
    }
  }

  return;
}
