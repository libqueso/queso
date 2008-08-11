/* libs/fp/src/uq1DScalarElement.C
 * 
 * Copyright (C) 2008 The PECOS Team, http://www.ices.utexas.edu/centers/pecos
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <uq1DScalarElement.h>
#include <math.h>

uq1DScalarElementClass::uq1DScalarElementClass(
  const uqEnvironmentClass& env,
  const uq1DNodeClass&      node0,
  const uq1DNodeClass&      node1,
        unsigned int        order)
  :
  m_env             (env),
  m_nodes           (2),//,NULL),
  m_order           (order),
  m_localDofs       (numDofs()),//,NULL),
  m_integrationOrder(2),
  m_w               (numIntegrationPoints(),0.),
  m_bcc             (numIntegrationPoints(),0.)
{
  m_nodes[0] = new uq1DNodeClass(node0);
  m_nodes[1] = new uq1DNodeClass(node1);

  m_localDofs[0] = new uq1DLocalDofClass(0.,UQ_AT_NODE_LOCAL_DOF_POS,node0.globalId());
  for (unsigned int i = 1; i < m_localDofs.size()-1; ++i) {
    m_localDofs[i] = new uq1DLocalDofClass((double)i/(double)(m_localDofs.size()-1.),UQ_INTERNAL_LOCAL_DOF_POS,UQ_INVALID_NODE_ID);
  }  
  m_localDofs[m_localDofs.size()-1] = new uq1DLocalDofClass(1.,UQ_AT_NODE_LOCAL_DOF_POS,node1.globalId());

  switch (m_integrationOrder) {
    case 2:
      m_w[0] = 5./18.;
      m_w[1] = 8./18.;
      m_w[2] = 5./18.;
      m_bcc[0] = 0.5*(-sqrt(3./5.)+1);
      m_bcc[1] = 0.5*(          0.+1);
      m_bcc[2] = 0.5*( sqrt(3./5.)+1);
    break;

    default:
      UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                        m_env.rank(),
                        "uq1DScalarElementClass::constructor()",
                        "no code for this integration order yet");
    break;
  }
}

uq1DScalarElementClass::~uq1DScalarElementClass()
{
  for (unsigned int i = 0; i < m_localDofs.size(); ++i) {
    if (m_localDofs[i]) delete m_localDofs[i];
  }
  delete m_nodes[1];
  delete m_nodes[0];
}

void
uq1DScalarElementClass::setGlobalDofId(
  unsigned int i,
  unsigned int globalDofId)
{
  m_localDofs[i]->setGlobalId(globalDofId);

  return;
}

unsigned int
uq1DScalarElementClass::numDofs() const
{
  unsigned int result = 0;

  switch (m_order) {
    case 2:
      result = 3;
    break;

    default:
      UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                        m_env.rank(),
                        "uq1DScalarElementClass::numDofs()",
                        "no code for this order yet");
    break;
  }

  return result;
}

const uq1DLocalDofClass&
uq1DScalarElementClass::localDof(unsigned int i)
{
  return *(m_localDofs[i]);
}

//double
//uq1DScalarElementClass::xOfLocalDof(unsigned int i)
//{
//  double x = m_nodes[0]->x() + m_localDofs[i]->bcc()*magnitude();
//
//  return x;
//}

double
uq1DScalarElementClass::magnitude()
{
  return (m_nodes[1]->x() - m_nodes[0]->x());
}

const uq1DNodeClass&
uq1DScalarElementClass::node(unsigned int i)
{
  UQ_FATAL_TEST_MACRO((i > 1),
                      m_env.rank(),
                      "uq1DScalarElementClass::node()",
                      "i out of range");

  return *(m_nodes[i]);
}

double
uq1DScalarElementClass::bcc(double x)
{
  UQ_FATAL_TEST_MACRO((x < m_nodes[0]->x()) || (m_nodes[1]->x() < x),
                      m_env.rank(),
                      "uq1DScalarElementClass::bcc()",
                      "x out of range");

  double result = (x - m_nodes[0]->x())/magnitude();

  return result;
}

double
uq1DScalarElementClass::phi(unsigned int i, double bcc)
{
  UQ_FATAL_TEST_MACRO((bcc < 0.) || (1. < bcc),
                      m_env.rank(),
                      "uq1DScalarElementClass::phi()",
                      "bcc out of range");

  double result = 0.;

  switch (m_order) {
    case 2:
      switch (i) {
        case 0:
          result = 2.*(bcc-0.5)*(bcc-1.);
        break;

        case 1:
          result = -4.*bcc*(bcc-1.);
        break;

        case 2:
          result = 2.*bcc*(bcc-0.5);
        break;

        default:
          UQ_FATAL_RC_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                            m_env.rank(),
                            "uq1DScalarElementClass::phi()",
                            "invalid i");
        break;
      }
    break;

    default:
      UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                        m_env.rank(),
                        "uq1DScalarElementClass::phi()",
                        "no code for this order yet");
    break;
  }

  return result;
}

double
uq1DScalarElementClass::gradPhi(unsigned int i, double bcc)
{
  UQ_FATAL_TEST_MACRO((bcc < 0.) || (1. < bcc),
                      m_env.rank(),
                      "uq1DScalarElementClass::gradPhi()",
                      "bcc out of range");

  double result = 0.;

  switch (m_order) {
    case 2:
      switch (i) {
        case 0:
          result = 4.*bcc-3.;
        break;

        case 1:
          result = -4.*(2.*bcc-1.);
        break;

        case 2:
          result = 4.*bcc-1.;
        break;

        default:
          UQ_FATAL_RC_MACRO(UQ_INVALID_INTERNAL_STATE_RC,
                            m_env.rank(),
                            "uq1DScalarElementClass::gradPhi()",
                            "invalid i");
        break;
      }
    break;

    default:
      UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                        m_env.rank(),
                        "uq1DScalarElementClass::gradPhi()",
                        "no code for this order yet");
    break;
  }

  return result;
}

unsigned int
uq1DScalarElementClass::numIntegrationPoints()
{
  unsigned int result = 0;

  switch (m_integrationOrder) {
    case 2:
      result = 3;
    break;

    default:
      UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                        m_env.rank(),
                        "uq1DScalarElementClass::numIntegrationPoints()",
                        "no code for this integration order yet");
    break;
  }

  return result;
}

void
uq1DScalarElementClass::print(std::ostream& os) const
{
  return;
}
