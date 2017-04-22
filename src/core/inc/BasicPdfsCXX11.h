//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2017 The PECOS Development Team
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

#ifndef QUESO_BASIC_PDFS_CXX11_H
#define QUESO_BASIC_PDFS_CXX11_H

#include <queso/Defines.h>

#ifdef QUESO_HAVE_CXX11

#include <queso/BasicPdfsBase.h>

namespace QUESO {

/*!
 * \file BasicPdfsCXX11.h
 * \brief Class for Basic PDFs using C++11 math functions
 */

/*!
 * \class BasicPdfsCXX11
 * \brief Base class for basic PDFs using C++ math functions
 *
 * \todo This class \b will acommodate the definition of a Joint PDF using
 * C++11 math functions to bootstrap our own PDF evaluations since std::random
 * doesn't have PDF evaluations.  It will ultimately be called by BaseJointPdf
 * and/or its derived classes (via m_env.basicPdfs()) during the construction
 * of Joint PDFs.
 */
class BasicPdfsCXX11 : public BasicPdfsBase
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor.
  BasicPdfsCXX11(int worldRank);

  //! Destructor.
  ~BasicPdfsCXX11();
  //@}

  //! @name Mathematical methods
  //@{
  //! Actual value of the Beta PDF.
  double betaPdfActualValue (double x, double alpha, double beta) const;

  //! Actual value of the Gamma PDF.
  double gammaPdfActualValue(double x, double a,     double b   ) const;
  //@}
protected:

private:
  //! Default constructor.
  BasicPdfsCXX11();
};

}  // End namespace QUESO

#endif  // QUESO_HAVE_CXX11

#endif  // QUESO_BASIC_PDFS_CXX11_H
