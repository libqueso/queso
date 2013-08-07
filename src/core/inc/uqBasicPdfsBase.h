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

#ifndef __UQ_BASIC_PDFS_BASE_H__
#define __UQ_BASIC_PDFS_BASE_H__

#include <uqDefines.h>
#include <iostream>

/*! \file uqBasicPdfsBaseClass.h
    \brief Class for Basic PDFs.
*/

/*! \class uqBasicPdfsBaseClass
    \brief TODO: Base class for basic PDFs (via either GSL or Boost). 
    
    \todo This class \b will acommodate the definition of a Joint PDF using either GSL or Boost 
    distributions. It will ultimately be called by uqBaseJointPdfClass and/or its
    derived classes (via m_env.basicPdfs()) during the construction of Joint PDFs.
*/
class uqBasicPdfsBaseClass
{
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  uqBasicPdfsBaseClass();
  
  //! Constructor.
  uqBasicPdfsBaseClass(int worldRank);
  
  //! Virtual destructor.
  virtual ~uqBasicPdfsBaseClass();
  //@}

  //! @name Mathematical methods
  //@{  
  //! TODO: Actual value of the Beta PDF (calculated via either Boost or GSL libraries). See template specialization.
  virtual double betaPdfActualValue (double x, double alpha, double beta) const = 0;
  
  //! TODO: Actual value of the Gamma PDF (calculated via either Boost or GSL libraries). See template specialization.
  virtual double gammaPdfActualValue(double x, double a,     double b   ) const = 0;
  //@}
protected:
  int m_worldRank;
};

#endif // __UQ_BASIC_PDFS_BASE_H__
