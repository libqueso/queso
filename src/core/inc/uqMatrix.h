//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, 
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
// 
// $Id$
//
//--------------------------------------------------------------------------

#ifndef __UQ_MATRIX_H__
#define __UQ_MATRIX_H__

#include <uqEnvironment.h>
#include <uqVector.h>
#include <iostream>

class uqMatrixClass
{
public:
           uqMatrixClass();
           uqMatrixClass(const uqBaseEnvironmentClass& env, const uqMapClass& map);
           uqMatrixClass(const uqMatrixClass& rhs);
  virtual ~uqMatrixClass();

  uqMatrixClass& operator= (const uqMatrixClass& rhs);
  uqMatrixClass& operator*=(double a);
  uqMatrixClass& operator+=(const uqMatrixClass& rhs);
  uqMatrixClass& operator-=(const uqMatrixClass& rhs);

    const uqBaseEnvironmentClass& env                 ()           const;
    const uqMapClass&             map                 ()           const;
          unsigned int            numOfProcsForStorage()           const;
          void                    setPrintHorizontally(bool value) const; // Yes, 'const'
          bool                    getPrintHorizontally()           const;
          void                    setInDebugMode      (bool value) const; // Yes, 'const'
          bool                    getInDebugMode      ()           const;

  virtual unsigned int            numRowsLocal        () const = 0;
  virtual unsigned int            numRowsGlobal       () const = 0;
  virtual unsigned int            numCols             () const = 0;
  virtual int                     chol                () = 0;
  virtual void                    zeroLower           (bool includeDiagonal = false) = 0;
  virtual void                    zeroUpper           (bool includeDiagonal = false) = 0;
  virtual void                    print               (std::ostream& os) const = 0;

protected:
  virtual void                    copy                (const uqMatrixClass& src);

  const   uqBaseEnvironmentClass& m_env;
#ifdef QUESO_CLASSES_INSTANTIATE_NEW_MAPS
  const   uqMapClass              m_map;
#else
  const   uqMapClass&             m_map;
#endif
  mutable bool                    m_printHorizontally;
  mutable bool                    m_inDebugMode;
};

#endif // __UQ_MATRIX_H__
