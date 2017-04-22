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

#include <queso/GslBlockMatrix.h>

namespace QUESO {

GslBlockMatrix::GslBlockMatrix(const std::vector<unsigned int> & blockSizes,
    const GslVector & v, double diagValue)
  : Matrix(v.env(), v.map()),
    m_vectorSpaces(blockSizes.size()),
    m_blocks(blockSizes.size())
{
  for (unsigned int i = 0; i < this->m_vectorSpaces.size(); i++) {
    this->m_vectorSpaces[i] = new VectorSpace<GslVector, GslMatrix>(m_env,
        "block_param_", blockSizes[i], NULL);
    this->m_blocks[i] = new GslMatrix(this->m_vectorSpaces[i]->zeroVector(),
        diagValue);
  }
}

GslBlockMatrix::~GslBlockMatrix()
{
  for (unsigned int i = 0; i < this->m_vectorSpaces.size(); i++) {
    delete this->m_blocks[i];
    delete this->m_vectorSpaces[i];
  }
}


unsigned int
GslBlockMatrix::numRowsLocal() const
{
  queso_not_implemented();
}

unsigned int
GslBlockMatrix::numRowsGlobal() const
{
  queso_not_implemented();
  return 0;
}

unsigned int
GslBlockMatrix::numCols() const
{
  queso_not_implemented();
  return 0;
}

int
GslBlockMatrix::chol()
{
  queso_not_implemented();
  return 0;
}

void
GslBlockMatrix::zeroLower(bool includeDiagonal)
{
  queso_not_implemented();
}

void
GslBlockMatrix::zeroUpper(bool includeDiagonal)
{
  queso_not_implemented();
}


GslMatrix &
GslBlockMatrix::getBlock(unsigned int i) const
{
  return *(this->m_blocks[i]);
}

unsigned int
GslBlockMatrix::numBlocks() const
{
  return this->m_blocks.size();
}

void
GslBlockMatrix::invertMultiply(const GslVector & b, GslVector & x) const
{
  unsigned int totalCols = 0;

  for (unsigned int i = 0; i < this->m_blocks.size(); i++) {
    totalCols += this->m_blocks[i]->numCols();
  }

  if (totalCols != b.sizeLocal()) {
    queso_error_msg("block matrix and rhs have incompatible sizes");
  }

  if (x.sizeLocal() != b.sizeLocal()) {
    queso_error_msg("solution and rhs have incompatible sizes");
  }

  unsigned int blockOffset = 0;

  // Do an invertMultiply for each block
  for (unsigned int i = 0; i < this->m_blocks.size(); i++) {
    GslVector blockRHS(this->m_vectorSpaces[i]->zeroVector());
    GslVector blockSol(this->m_vectorSpaces[i]->zeroVector());

    // Be sure to copy over the RHS to the right sized vector
    for (unsigned int j = 0; j < this->m_blocks[i]->numCols(); j++) {
      blockRHS[j] = b[blockOffset + j];
    }

    // Solve
    this->m_blocks[i]->invertMultiply(blockRHS, blockSol);

    // Be sure to copy the block solution back to the global solution vector
    for (unsigned int j = 0; j < this->m_blocks[i]->numCols(); j++) {
      x[blockOffset + j] = blockSol[j];
    }

    // Remember to increment the offset so we don't lose our place for the next
    // block
    blockOffset += this->m_blocks[i]->numCols();
  }
}

void
GslBlockMatrix::print(std::ostream& os) const
{
  for (unsigned int i = 0; i < this->numBlocks(); i++) {
    this->getBlock(i).print(os);
  }
}

std::ostream&
operator<<(std::ostream& os, const GslBlockMatrix & obj)
{
  obj.print(os);

  return os;
}

}  // End namespace QUESO
