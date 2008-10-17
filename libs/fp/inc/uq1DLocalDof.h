/* libs/fp/inc/uq1DLocalDof.h
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

#ifndef __UQ_1D_LOCAL_DOFS_H__
#define __UQ_1D_LOCAL_DOFS_H__

#include <iostream>
#include <uqDefines.h>

enum uqLocalDofPositionEnum {
  UQ_INTERNAL_LOCAL_DOF_POS = 0,
  UQ_AT_NODE_LOCAL_DOF_POS,
  UQ_AT_EDGE_LOCAL_DOF_POS,
  UQ_AT_FACE_LOCAL_DOF_POS
};

class uq1DLocalDofClass
{
public:
  uq1DLocalDofClass(double                 bcc,
                    uqLocalDofPositionEnum localPositionType,
                    unsigned int           globalIdOfRespectiveNode);
  uq1DLocalDofClass(const uq1DLocalDofClass& obj);
 ~uq1DLocalDofClass();

  uq1DLocalDofClass& operator=(const uq1DLocalDofClass& rhs);

  double       bcc                     () const;
  unsigned int globalId                () const;
  void         setGlobalId             (unsigned int globalId);
  unsigned int globalIdOfRespectiveNode() const;
  void         print                   (std::ostream& os) const;

protected:
  void         copy                    (const uq1DLocalDofClass& src);

  double                 m_bcc;
  uqLocalDofPositionEnum m_myLocalPositionType;
  unsigned int           m_globalIdOfRespectiveNode;
  unsigned int           m_myGlobalId;
};

#endif // __UQ_1D_LOCAL_DOFS_H__
