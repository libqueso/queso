#include <uq1DProblem.h>
#include <uqGslMatrix.h>
#include <uqGslVector.h>
#include <uqDefines.h>

uq1DProblemClass::uq1DProblemClass(
  const uqEnvironmentClass& env,
  double a,
  double b,
  double c,
  double f,
  double g,
  double h)
  :
  m_env        (env),
  m_a          (a),
  m_b          (b),
  m_c          (c),
  m_f          (f),
  m_g          (g),
  m_h          (h),
  m_coefA      (0.5*m_f/m_c),
  m_coefB      (m_h - 0.5*m_f/m_c - m_g),
  m_coefC      (m_g),
  m_numDiv     (0),
  m_order      (0),
  m_numNodes   (0),
  m_numElements(0),
  m_nodes      (0),
  m_elements   (0),
  m_globalU    (NULL)
{
}

uq1DProblemClass::~uq1DProblemClass()
{
  if (m_globalU)    delete m_globalU;
  for (unsigned int i = 0; i < m_elements.size(); ++i) {
    if (m_elements[i]) delete m_elements[i];
  }
  for (unsigned int i = 0; i < m_nodes.size(); ++i) {
    if (m_nodes[i]) delete m_nodes[i];
  }
}

void
uq1DProblemClass::femSolve(
  unsigned int numDiv,
  unsigned int order)
{
  m_numDiv = numDiv;
  m_order  = order;

  //////////////////////////////////////////////////
  // Initialize some important quantities
  //////////////////////////////////////////////////
  unsigned int m_numNodes    = m_numDiv+1;
  unsigned int m_numElements = m_numDiv;

  //////////////////////////////////////////////////
  // Instantiate nodes
  //////////////////////////////////////////////////
  m_nodes.resize(m_numNodes,NULL);
  m_nodes[0] = new uq1DNodeClass(0,
                                 m_a,
                                 UQ_AT_BOUNDARY_NODE_NODE_POS,
                                 UQ_DIRICHLET_BC,
                                 m_g);
  for (unsigned int i = 1; i < m_numNodes-1; ++i) {
    double x = ( ((double) i)*m_b + ((double) (m_numDiv - i))*m_a )/(double) m_numDiv;
    m_nodes[i] = new uq1DNodeClass(i,
                                   x,
                                   UQ_INTERNAL_NODE_POS,
                                   UQ_NONE_BC,
                                   0.);
  }
  m_nodes[m_numNodes-1] = new uq1DNodeClass(m_numNodes-1,
                                            m_b,
                                            UQ_AT_BOUNDARY_NODE_NODE_POS,
                                            UQ_DIRICHLET_BC,
                                            m_h);

  std::set<unsigned int> m_setOfDirichletNodeIds;
  m_setOfDirichletNodeIds.insert(0);
  m_setOfDirichletNodeIds.insert(m_numNodes-1);

  //////////////////////////////////////////////////
  // Instantiate elements
  //////////////////////////////////////////////////
  m_elements.resize(m_numElements,NULL);
  for (unsigned int e = 0; e < m_numElements; ++e) {
    m_elements[e] = new uq1DScalarElementClass(m_env,*(m_nodes[e]),*(m_nodes[e+1]),m_order);
  }

  //////////////////////////////////////////////////
  // Determine global number of dofs
  //////////////////////////////////////////////////
  unsigned int m_globalNumDofs = 1;
  for (unsigned int e = 0; e < m_numElements; ++e) {
    m_globalNumDofs += (m_elements[e]->numDofs()-1);
  }

  //////////////////////////////////////////////////
  // Globally number all dofs in all elements
  //////////////////////////////////////////////////
  unsigned globalDofId = 0;
  for (unsigned int e = 0; e < m_numElements; ++e) {
    uq1DScalarElementClass* elem = m_elements[e];
    for (unsigned int i = 0; i < elem->numDofs(); ++i) {
      elem->setGlobalDofId(i,globalDofId);
      globalDofId++;
    }
    globalDofId--;
  }

  //////////////////////////////////////////////////
  // Assemble matrices
  //////////////////////////////////////////////////
  uqGslMatrixClass m_globalMatrix(m_env,m_globalNumDofs,m_globalNumDofs);
  for (unsigned int e = 0; e < m_numElements; ++e) {
    uq1DScalarElementClass* elem = m_elements[e];
    elem->updateStiffnessMatrix(m_globalMatrix,
                                m_c,
                                m_setOfDirichletNodeIds);
  }
  //std::cout << "m_globalMatrix = " << m_globalMatrix
  //          << std::endl;

  //////////////////////////////////////////////////
  // Assemble rhs
  //////////////////////////////////////////////////
  uqGslVectorClass m_globalRhs(m_env,m_globalNumDofs);
  for (unsigned int e = 0; e < m_numElements; ++e) {
    uq1DScalarElementClass* elem = m_elements[e];
    elem->updateRhs(m_globalRhs,
                    m_f,
                    m_setOfDirichletNodeIds);
  }
  uqGslVectorClass m_DirichletRhs(m_env,m_globalNumDofs);
  m_DirichletRhs[0                ] = m_g;
  m_DirichletRhs[m_globalNumDofs-1] = m_h;
  m_globalRhs += m_DirichletRhs;
  //std::cout << "final m_globalRhs = " << m_globalRhs
  //          << std::endl;

  //////////////////////////////////////////////////
  // Solve system of linear equations
  //////////////////////////////////////////////////
  m_globalU = new uqGslVectorClass(m_env,m_globalNumDofs);
  m_globalMatrix.invertMultiply(m_globalRhs, *m_globalU);
  //std::cout << "m_globalU = " << *m_globalU
  //          << std::endl;

#if 0
  //////////////////////////////////////////////////
  // Compute exact solution
  //////////////////////////////////////////////////
  double coefA = (0.5*m_f/m_c);
  double coefB = (m_h - 0.5*m_f/m_c - m_g);
  double coefC = m_g;
  uqGslVectorClass m_x     (m_env,m_globalNumDofs);
  uqGslVectorClass m_exactU(m_env,m_globalNumDofs);
  for (unsigned int e = 0; e < m_numElements; ++e) {
    uq1DScalarElementClass* elem = m_elements[e];
    for (unsigned int i = 0; i < elem->numDofs(); ++i) {
      unsigned int ii = elem->localDof(i).globalId();
      double x        = elem->xOfLocalDof(i);
      double value    = coefA*x*x + coefB*x + coefC;
      m_x     [ii] = x;
      m_exactU[ii] = value;
    }
  }
  //std::cout << "m_x = "            << m_x
  //          << std::endl;
  //std::cout << "Exact solution = " << coefA
  //          << "*x^2 + "           << coefB
  //          << "*x + "             << coefC
  //          << std::endl;
  //std::cout << "m_exactU = "       << m_exactU
  //          << std::endl;

  double exactNorm2 = m_exactU.norm2();
  m_exactU -= *m_globalU;
  double diffNorm2 = m_exactU.norm2();
  std::cout << "||diff||_2/||m_exactU||_2 = " << diffNorm2
            << " / "                          << exactNorm2
            << "  = "                         << m_exactU.norm2()/exactNorm2
            << std::endl;
#endif

  return;
}

double
uq1DProblemClass::uValue(double x) const
{
  if ((x < m_a) || (x > m_b)) {
    std::cerr << "In uq1DProblemClass::uValue()"
              << ": requested x ("             << x
              << ") is out of the intervarl (" << m_a
              << ", "                          << m_b
              << ")"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(((x < m_a) || (x > m_b)),
                      m_env.rank(),
                      "uq1DProblemClass::uValue()",
                      "x out of range");

  unsigned int e = 0;
  bool found = false;
  while ((found == false) && (e < m_elements.size())) {
    if ((m_elements[e]->node(0).x() <= x) &&
        (x <= m_elements[e]->node(1).x())) {
      found = true;
    }
    else {
      e++;
    }
  }
  UQ_FATAL_TEST_MACRO(((x < m_a) || (x > m_b)),
                      m_env.rank(),
                      "uq1DProblemClass::uValue()",
                      "element not found");

  uq1DScalarElementClass* elem = m_elements[e];
  double result = 0.;
  for (unsigned int i = 0; i < elem->numDofs(); ++i) {
    unsigned int ii = elem->localDof(i).globalId();
    double bcc = elem->bcc(x);
    result += (*m_globalU)[ii] * elem->phi(i,bcc);
  }

  return result;
}

double
uq1DProblemClass::uExact(double x) const
{
  if ((x < m_a) || (x > m_b)) {
    std::cerr << "In uq1DProblemClass::uExact()"
              << ": requested x ("             << x
              << ") is out of the intervarl (" << m_a
              << ", "                          << m_b
              << ")"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(((x < m_a) || (x > m_b)),
                      m_env.rank(),
                      "uq1DProblemClass::uExact()",
                      "x out of range");

  double value = m_coefA*x*x + m_coefB*x + m_coefC;

  return value;
}

double
uq1DProblemClass::a() const
{
  return m_a;
}

double
uq1DProblemClass::b() const
{
  return m_b;
}

double
uq1DProblemClass::c() const
{
  return m_c;
}

double
uq1DProblemClass::f() const
{
  return m_f;
}

double
uq1DProblemClass::g() const
{
  return m_g;
}

double
uq1DProblemClass::h() const
{
  return m_h;
}
