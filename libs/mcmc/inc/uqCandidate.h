#ifndef __UQ_CANDIDATE_H__
#define __UQ_CANDIDATE_H__

#include <uqParameter.h>
#include <uqEnvironment.h>

template <class V>
class uqCandidateClass
{
public:
  uqCandidateClass(const uqEnvironmentClass& env);
  uqCandidateClass(const uqEnvironmentClass& env,
                   const V& paramValues,
                   bool     outOfBounds,
                   double   m2lPrior,
                   const V& m2lLikelihoodResults,
                   const V& lrSigma2,
                   double   logPosterior);
  uqCandidateClass(const uqCandidateClass<V>& rhs);
 ~uqCandidateClass();

  uqCandidateClass<V>& operator= (const uqCandidateClass<V>& rhs);

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
uqCandidateClass<V>::uqCandidateClass(const uqEnvironmentClass& env)
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
uqCandidateClass<V>::uqCandidateClass(
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
uqCandidateClass<V>::uqCandidateClass(const uqCandidateClass<V>& rhs)
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
uqCandidateClass<V>::~uqCandidateClass()
{
  if (m_lrSigma2)             delete m_lrSigma2;
  if (m_m2lLikelihoodResults) delete m_m2lLikelihoodResults;
  if (m_paramValues)          delete m_paramValues;
}

template <class V>
uqCandidateClass<V>&
uqCandidateClass<V>::operator=(const uqCandidateClass<V>& rhs)
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
uqCandidateClass<V>::paramValues() const
{
  UQ_FATAL_TEST_MACRO((m_paramValues == NULL),
                      m_env.rank(),
                      "uqCandidateClass<V>::paramValues()",
                      "m_paramValues is NULL");
  return *m_paramValues;
}

template <class V>
bool
uqCandidateClass<V>::outOfBounds() const
{
  return m_outOfBounds;
}

template <class V>
double
uqCandidateClass<V>::m2lPrior() const
{
  return m_m2lPrior;
}

template <class V>
const V&    
uqCandidateClass<V>::m2lLikelihood() const
{
  return *m_m2lLikelihoodResults;
}

template <class V>
const V&
uqCandidateClass<V>::lrSigma2() const
{
  return *m_lrSigma2;
}

template <class V>
double
uqCandidateClass<V>::logPosterior() const
{
  return m_logPosterior;
}

template <class V>
void
uqCandidateClass<V>::set(
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
std::ostream& operator<<(std::ostream& os, const uqCandidateClass<V>& obj)
{
  obj.print(os);

  return os;
}

#endif // __UQ_CANDIDATE_H__
