#ifndef __UQ_1D_SCALAR_ELEMENT_H__
#define __UQ_1D_SCALAR_ELEMENT_H__

#include <uq1DNode.h>
#include <uq1DLocalDof.h>
#include <uqEnvironment.h>
#include <uqGslMatrix.h>
#include <uqGslVector.h>
#include <uqDefines.h>
#include <iostream>

class uq1DScalarElementClass
{
public:
  uq1DScalarElementClass(const uqEnvironmentClass& env,
                         const uq1DNodeClass&      node0,
                         const uq1DNodeClass&      node1,
                               unsigned int        order);
 ~uq1DScalarElementClass();

  template<class M>
  void                     updateStiffnessMatrix(M&     globalMatrix,
                                                 double cTerm,
                                                 const std::set<unsigned int>& setOfDirichletNodeIds);
  template<class V>
  void                     updateRhs            (V&     globalRhs,
                                                 double fTerm,
                                                 const std::set<unsigned int>& setOfDirichletNodeIds);
  void                     setGlobalDofId       (unsigned int i, unsigned int globalDofId);
  unsigned int             numDofs              () const;
  const uq1DLocalDofClass& localDof             (unsigned int i);
  //double                   xOfLocalDof          (unsigned int i);
  const uq1DNodeClass&     node                 (unsigned int i);
  double                   bcc                  (double x);
  double                   phi                  (unsigned int i, double bcc);
  void                     print                (std::ostream& os) const;

protected:
  double                   magnitude            ();
  double                   gradPhi              (unsigned int i, double bcc);
  unsigned int             numIntegrationPoints ();

  const uqEnvironmentClass&       m_env;
  std::vector<uq1DNodeClass*>     m_nodes;
  unsigned int                    m_order;
  std::vector<uq1DLocalDofClass*> m_localDofs;
  unsigned int                    m_integrationOrder;
  std::vector<double>             m_w;
  std::vector<double>             m_bcc;
};

template<class M>
void
uq1DScalarElementClass::updateStiffnessMatrix(
  M&     globalMatrix,
  double cTerm,
  const std::set<unsigned int>& setOfDirichletNodeIds)
{
  // The matrix below is always Gsl, since it is local
  uqGslMatrixClass localMatrix(globalMatrix.env(),numDofs(),numDofs());

  // Compute local matrix
  for (unsigned int i = 0; i < numDofs(); ++i) {
    for (unsigned int j = 0; j < numDofs(); ++j) {
      double value = 0.;
      for (unsigned int k = 0; k < numIntegrationPoints(); ++k) {
        double bcc = m_bcc[k];
        value += -m_w[k]*cTerm*gradPhi(i,bcc)*gradPhi(j,bcc)/magnitude();
      }
      localMatrix(i,j) = value;
    }
  }

  // Add newly computed values to their respective global positions
  for (unsigned int i = 0; i < numDofs(); ++i) {
    unsigned int ii = m_localDofs[i]->globalId();
    unsigned int respi = m_localDofs[i]->globalIdOfRespectiveNode();
    if ((respi                             != UQ_INVALID_NODE_ID         ) &&
        (setOfDirichletNodeIds.find(respi) != setOfDirichletNodeIds.end())) {
      globalMatrix(ii,ii) = 1.;
    }
    else {
      for (unsigned int j = 0; j < numDofs(); ++j) {
        unsigned int jj = m_localDofs[j]->globalId();
        //unsigned int respj = m_localDofs[j]->globalIdOfRespectiveNode();
        //if ((respj                             != UQ_INVALID_NODE_ID         ) &&
        //    (setOfDirichletNodeIds.find(respj) != setOfDirichletNodeIds.end())) {
        //  // Do nothing especial
        //}
        globalMatrix(ii,jj) += localMatrix(i,j);
      }
    }
  }

  return;
}

template<class V>
void
uq1DScalarElementClass::updateRhs(
  V&     globalRhs,
  double fTerm,
  const std::set<unsigned int>& setOfDirichletNodeIds)
{
  // The vector below is always Gsl, since it is local
  uqGslVectorClass localRhs(globalRhs.env(),numDofs());

  // Compute local rhs
  for (unsigned int i = 0; i < numDofs(); ++i) {
    double value = 0.;
    for (unsigned int k = 0; k < numIntegrationPoints(); ++k) {
      double bcc = m_bcc[k];
      value += m_w[k]*fTerm*phi(i,bcc)*magnitude();
    }
    localRhs[i] = value;
  }

  for (unsigned int i = 0; i < numDofs(); ++i) {
    unsigned int ii    = m_localDofs[i]->globalId();
    unsigned int respi = m_localDofs[i]->globalIdOfRespectiveNode();
    if ((respi                             != UQ_INVALID_NODE_ID         ) &&
        (setOfDirichletNodeIds.find(respi) != setOfDirichletNodeIds.end())) {
      // Keep global position of vector with value 0.
    }
    else {
      globalRhs[ii] += localRhs[i];
    }
  }

  return;
}
#endif // __UQ_1D_SCALAR_ELEMENT_H__
