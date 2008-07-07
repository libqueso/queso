/* uq/libs/mcmc/inc/uqTrilinosParamSpace.h
 *
 * Copyright (C) 2008 The PECOS Team, http://www.ices.utexas.edu/centers/pecos

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

#ifndef __UQ_TRILINOS_PARAM_SPACE_H__
#define __UQ_TRILINOS_PARAM_SPACE_H__

#include <uqTrilinosVector.h>
#include <uqTrilinosMatrix.h>
#include <uqParamSpace.h>
#include <Epetra_MpiComm.h>

class uqTrilinosParamSpaceClass : public uqParamSpaceClass<uqTrilinosVectorClass,uqTrilinosMatrixClass>
{
public:
  uqTrilinosParamSpaceClass(const Epetra_MpiComm& comm, unsigned int dimension);
 ~uqTrilinosParamSpaceClass();

  uqTrilinosVectorClass* newVector             () const;
  uqTrilinosMatrixClass* newMatrix             () const;

  void                   print                 (std::ostream& os) const;

  const Epetra_Map&      map                   () const;

protected:
  void                   createInitialValues   () const;
  void                   createMinValues       () const;
  void                   createMaxValues       () const;
  void                   createPriorMuValues   () const;
  void                   createPriorSigmaValues() const;

private:
  const Epetra_Map* m_map;
};

std::ostream& operator<<(std::ostream& os, const uqTrilinosParamSpaceClass& space);

#endif // __UQ_TRILINOS_PARAM_SPACE_H__
