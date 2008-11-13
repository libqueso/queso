/* libs/fp/inc/uq1DProblem.h
 * 
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu/
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_1D_PROBLEM__
#define __UQ_1D_PROBLEM__

#include <uq1DScalarElement.h>
#include <uq1DNode.h>
#include <uqEnvironment.h>

class uq1DProblemClass
{
public:
  uq1DProblemClass(const uqBaseEnvironmentClass& env,
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
  const uqBaseEnvironmentClass& m_env;
  const Epetra_Map*         m_map;
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
