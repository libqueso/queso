//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
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

#ifndef UQ_GSL_BLOCK_MATRIX_H
#define UQ_GSL_BLOCK_MATRIX_H

/*!
 * \file GslBlockMatrix.h
 * \brief QUESO matrix class using GSL.
 */

#include <queso/Matrix.h>

namespace QUESO {

/*!
 * \class GslBlockMatrix
 * \brief Class for matrix operations using GSL library.
 *
 * This class creates and provides basic support for matrices of templated
 * type as a specialization of Matrix using GSL matrices, which are defined
 * by an encapsulated gsl_matrix structure.
 */

class GslBlockMatrix : public Matrix
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Shaped Constructor: creates a square matrix with size \c v.sizeLocal() and diagonal values all equal to \c diagValue.
  GslBlockMatrix(const std::vector<unsigned int> & blockSizes,
      double diagValue);

  //! Destructor
  ~GslBlockMatrix();
  //@}

  //! Return block \c i in the block diagonal matrix
  GslBlockMatrix & getBlock(unsigned int i) const;

  //! Return the number of blocks in the block diagonal matrix
  unsigned int numBlocks() const;

  //! This function calculates the inverse of \c this matrix, multiplies it with vector \c b and stores the result in vector \c x.
  /*!
   * It checks for a previous LU decomposition of \c this matrix and does not
   * recompute it if m_MU != NULL.
   */
  void invertMultiply (const GslVector & b, GslVector & x) const;

  //! @name I/O methods
  //@{
  //! Print method. Defines the behavior of operator<< inherited from the Object class.
  void print (std::ostream & os) const;
  //@}

private:

  //! Default Constructor
  /*! Creates an empty matrix vector of no dimension. It should not be used by user.*/
  GslBlockMatrix();
};

std::ostream & operator<<(std::ostream & os, const GslBlockMatrix & obj);

}  // End namespace QUESO

#endif // UQ_GSL_BLOCK_MATRIX_H
