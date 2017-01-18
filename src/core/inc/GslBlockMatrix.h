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

#ifndef UQ_GSL_BLOCK_MATRIX_H
#define UQ_GSL_BLOCK_MATRIX_H

/*!
 * \file GslBlockMatrix.h
 * \brief QUESO block matrix class using GSL.
 */

#include <vector>
#include <queso/Environment.h>
#include <queso/Matrix.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorSpace.h>

namespace QUESO {

/*!
 * \class GslBlockMatrix
 * \brief Class for representing block matrices using GSL library.
 *
 * This class provides basic 'invertMultiply' support for matrices of block
 * diagonal structure.  Each block is implemented as a GslMatrix object.
 */

class GslBlockMatrix : Matrix
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Creates a square matrix with size defined by \c v and diagonal values all equal to \c diagValue.
  /*!
   * The \c blockSizes array determines the sizes of each (square) block.
   */
  GslBlockMatrix(const std::vector<unsigned int> & blockSizes,
      const GslVector & v, double diagValue);

  //! Destructor
  ~GslBlockMatrix();
  //@}

  //! Not implemented yet
  virtual unsigned int numRowsLocal() const;

  //! Not implemented yet
  virtual unsigned int numRowsGlobal() const;

  //! Not implemented yet
  virtual unsigned int numCols() const;

  //! Not implemented yet
  virtual int chol();

  //! Not implemented yet
  virtual void zeroLower(bool includeDiagonal=false);

  //! Not implemented yet
  virtual void zeroUpper(bool includeDiagonal=false);

  //! Return block \c i in the block diagonal matrix
  GslMatrix & getBlock(unsigned int i) const;

  //! Return the number of blocks in the block diagonal matrix
  unsigned int numBlocks() const;

  //! This function calculates the inverse of \c this matrix, multiplies it with vector \c b and stores the result in vector \c x.
  /*!
   * It checks for a previous LU decomposition of each block matrix and does
   * not recompute it if m_MU != NULL for each block.
   */
  void invertMultiply(const GslVector & b, GslVector & x) const;

  //! @name I/O methods
  //@{
  //! Print method. Defines the behavior of operator<< inherited from the Object class.
  virtual void print (std::ostream & os) const;
  //@}

private:
  std::vector<VectorSpace<GslVector, GslMatrix> *> m_vectorSpaces;
  std::vector<GslMatrix *> m_blocks;
};

std::ostream & operator<<(std::ostream & os, const GslBlockMatrix & obj);

}  // End namespace QUESO

#endif // UQ_GSL_BLOCK_MATRIX_H
