#ifndef __UQ_1D_PROBLEM__
#define __UQ_1D_PROBLEM__

#include <uq1DScalarElement.h>
#include <uq1DNode.h>
#include <uqEnvironment.h>

class uq1DProblemClass
{
public:
  uq1DProblemClass(const uqEnvironmentClass& env,
                   double a,
                   double b,
                   double c,
                   double f,
                   double g,
                   double h);
 ~uq1DProblemClass();

  void   femSolve(unsigned int numDiv,
                  unsigned int order);
  double uValue  (double x) const;
  double uExact  (double x) const;
  double a       () const;
  double b       () const;
  double c       () const;
  double f       () const;
  double g       () const;
  double h       () const;

protected:
  const uqEnvironmentClass& m_env;
  double m_a;
  double m_b;
  double m_c;
  double m_f;
  double m_g;
  double m_h;
  double m_coefA;
  double m_coefB;
  double m_coefC;

  unsigned int m_numDiv;
  unsigned int m_order;

  unsigned int                         m_numNodes;
  unsigned int                         m_numElements;
  std::vector<uq1DNodeClass*>          m_nodes;
  std::vector<uq1DScalarElementClass*> m_elements;
  uqGslVectorClass*                    m_globalU;

};

#endif // __UQ_1D_PROBLEM__
