#ifndef __UQ_CHAIN_POSITION_H__
#define __UQ_CHAIN_POSITION_H__

#include <uqParameter.h>
#include <uqEnvironment.h>

template <class V>
class uqChainPositionClass
{
public:
  uqChainPositionClass(const uqEnvironmentClass& env);
  uqChainPositionClass(const uqEnvironmentClass& env,
                       const V& paramValues,
                       bool     outOfBounds,
                       double   m2lPrior,
                       const V& m2lLikelihoodResults,
                       const V& lrSigma2,
                       double   logPosterior);
  uqChainPositionClass(const uqChainPositionClass<V>& rhs);
 ~uqChainPositionClass();

  uqChainPositionClass<V>& operator= (const uqChainPositionClass<V>& rhs);

  const V& paramValues  () const;
  bool     outOfBounds  () const;
  double   m2lPrior     () const;
  const V& m2lLikelihood() const;
  const V& lrSigma2     () const;
  double   logPosterior () const;

  void     set          (const V& paramValues,
                         bool     outOfBounds,
                         double   m2lPrior,
                         const V& m2lLikelihoodResults,
                         const V& lrSigma2,
                         double   logPosterior);

  void     print        (std::ostream& os) const;

private:
  const uqEnvironmentClass& m_env;
  V*     m_paramValues;
  bool   m_outOfBounds;
  double m_m2lPrior;
  V*     m_m2lLikelihoodResults;
  V*     m_lrSigma2;
  double m_logPosterior;
};

template <class V>
uqChainPositionClass<V>::uqChainPositionClass(const uqEnvironmentClass& env)
  :
  m_env                 (env),
  m_paramValues         (NULL),
  m_outOfBounds         (false),
  m_m2lPrior            (0.),
  m_m2lLikelihoodResults(NULL),
  m_lrSigma2            (NULL),
  m_logPosterior        (0.)
{
}

template <class V>
uqChainPositionClass<V>::uqChainPositionClass(
  const uqEnvironmentClass& env,
  const V& paramValues,
  bool     outOfBounds,
  double   m2lPrior,
  const V& m2lLikelihoodResults,
  const V& lrSigma2,
  double   logPosterior)
  :
  m_env                 (env),
  m_paramValues         (new V(paramValues)),
  m_outOfBounds         (outOfBounds),
  m_m2lPrior            (m2lPrior),
  m_m2lLikelihoodResults(new V(m2lLikelihoodResults)),
  m_lrSigma2            (new V(lrSigma2)),
  m_logPosterior        (logPosterior)
{
}

template <class V>
uqChainPositionClass<V>::uqChainPositionClass(const uqChainPositionClass<V>& rhs)
  :
  m_env                 (rhs.m_env                         ),
  m_paramValues         (new V(*rhs.m_paramValues)         ),
  m_outOfBounds         (rhs.m_outOfBounds                 ),
  m_m2lPrior            (rhs.m_m2lPrior                    ),
  m_m2lLikelihoodResults(new V(*rhs.m_m2lLikelihoodResults)),
  m_lrSigma2            (new V(*rhs.m_lrSigma2)            ),
  m_logPosterior        (rhs.m_logPosterior                )
{
}

template <class V>
uqChainPositionClass<V>::~uqChainPositionClass()
{
  if (m_lrSigma2)             delete m_lrSigma2;
  if (m_m2lLikelihoodResults) delete m_m2lLikelihoodResults;
  if (m_paramValues)          delete m_paramValues;
}

template <class V>
uqChainPositionClass<V>&
uqChainPositionClass<V>::operator=(const uqChainPositionClass<V>& rhs)
{
  if (m_paramValues == NULL) m_paramValues = new V(*rhs.m_paramValues);
  else                      *m_paramValues = *rhs.m_paramValues;
  m_outOfBounds   = rhs.m_outOfBounds;
  m_m2lPrior      = rhs.m_m2lPrior;
  if (m_m2lLikelihoodResults == NULL) m_m2lLikelihoodResults = new V(*rhs.m_m2lLikelihoodResults);
  else                               *m_m2lLikelihoodResults = *rhs.m_m2lLikelihoodResults;
  if (m_lrSigma2             == NULL) m_lrSigma2             = new V(*rhs.m_lrSigma2);
  else                               *m_lrSigma2             = *rhs.m_lrSigma2;
  m_logPosterior  = rhs.m_logPosterior;

  return *this;
}

template <class V>
const V&
uqChainPositionClass<V>::paramValues() const
{
  UQ_FATAL_TEST_MACRO((m_paramValues == NULL),
                      m_env.rank(),
                      "uqChainPositionClass<V>::paramValues()",
                      "m_paramValues is NULL");
  return *m_paramValues;
}

template <class V>
bool
uqChainPositionClass<V>::outOfBounds() const
{
  return m_outOfBounds;
}

template <class V>
double
uqChainPositionClass<V>::m2lPrior() const
{
  return m_m2lPrior;
}

template <class V>
const V&    
uqChainPositionClass<V>::m2lLikelihood() const
{
  return *m_m2lLikelihoodResults;
}

template <class V>
const V&
uqChainPositionClass<V>::lrSigma2() const
{
  return *m_lrSigma2;
}

template <class V>
double
uqChainPositionClass<V>::logPosterior() const
{
  return m_logPosterior;
}

template <class V>
void
uqChainPositionClass<V>::set(
  const V& paramValues,
  bool     outOfBounds,
  double   m2lPrior,
  const V& m2lLikelihoodResults,
  const V& lrSigma2,
  double   logPosterior)
{
  if (m_paramValues == NULL) m_paramValues = new V(paramValues);
  else                      *m_paramValues = paramValues;
  m_outOfBounds   = outOfBounds;
  m_m2lPrior      = m2lPrior;
  if (m_m2lLikelihoodResults == NULL) m_m2lLikelihoodResults = new V(m2lLikelihoodResults);
  else                               *m_m2lLikelihoodResults = m2lLikelihoodResults;
  if (m_lrSigma2             == NULL) m_lrSigma2             = new V(lrSigma2);
  else                               *m_lrSigma2             = lrSigma2;
  m_logPosterior  = logPosterior;

  return;
}

template <class V>
std::ostream& operator<<(std::ostream& os, const uqChainPositionClass<V>& obj)
{
  obj.print(os);

  return os;
}

#endif // __UQ_CHAIN_POSITION_H__
