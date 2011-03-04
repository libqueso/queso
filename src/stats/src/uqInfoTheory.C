//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011 The PECOS Development Team
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
// $Id:$
//
//--------------------------------------------------------------------------

#include <uqInfoTheory.h>

#ifdef QUESO_HAS_ANN

//*****************************************************
// Function: distANN_XY
//*****************************************************
void distANN_XY( const ANNpointArray dataX, const ANNpointArray dataY, 
		 double* distsXY, 
		 unsigned int dimX, unsigned int dimY, 
		 unsigned int xN, unsigned int yN, 
		 unsigned int k, double eps ) 
{

  ANNkd_tree* kdTree;
  ANNdistArray nnDist;
  ANNidxArray nnIdx;
	
  // Allocate memory
  nnIdx = new ANNidx[ k+1 ];
  nnDist = new ANNdist[ k+1 ];
  kdTree = new ANNkd_tree( dataY, yN, dimY );
  
  // Get the distances to all the points
  for( unsigned int i = 0; i < xN ; i++ ) {
    kdTree->annkSearch( dataX[ i ], k+1, nnIdx, nnDist, eps );
    distsXY[ i ] = nnDist[ k ];
  }
	
  // Deallocate memory
  delete [] nnIdx;
  delete [] nnDist;
  delete kdTree;
  annClose();

  return;
}

//*****************************************************
// Function: normalizeData
//*****************************************************
void normalizeANN_XY( ANNpointArray dataXY, unsigned int dimXY,
		      ANNpointArray dataX, unsigned int dimX,
		      ANNpointArray dataY, unsigned int dimY,
		      unsigned int N ) 
{

  ANNpoint meanXY, stdXY;

  meanXY = annAllocPt(dimXY); // is init with 0
  stdXY = annAllocPt(dimXY); // is init with 0

  // get mean
  for( unsigned int i = 0; i < N; i++ ) {
    for( unsigned int j = 0; j < dimX; j++ ) {
      meanXY[ j ] += dataXY[ i ][ j ];
    }
    for( unsigned int j = 0; j < dimY; j++ ) {
      meanXY[ dimX + j ] += dataXY[ i ][ dimX + j ];
    }
  }
  for( unsigned int j = 0; j < dimXY; j++ ) {
    meanXY[ j ] = meanXY[ j ] / (double)N;
  }

  // get standard deviation
  for( unsigned int i = 0; i < N; i++ ) {
    for( unsigned int j = 0; j < dimXY; j++ ) {
      stdXY[ j ] += pow( dataXY[ i ][ j ] - meanXY[ j ], 2.0 );
    }
  }
  for( unsigned int j = 0; j < dimXY; j++ ) {
    stdXY[ j ] = sqrt( stdXY[ j ] / ((double)N-1.0) );
  }

  /*
  std::cout << "Mean RV - ";
  annPrintPt( meanXY, dimXY, std::cout );
  std::cout << std::endl << "Std RV - ";
  annPrintPt( stdXY, dimXY, std::cout );
  std::cout << std::endl;
  */

  // get normalized points and populate marginals
  for( unsigned int i = 0; i < N; i++ ) {
    // normalize joint
    for( unsigned int j = 0; j < dimXY; j++ ) {
      dataXY[ i ][ j ] = ( dataXY[ i ][ j ] - meanXY[ j ] ) / stdXY[ j ];
    }

    // populate marginals
    for( unsigned int j = 0; j < dimX; j++ ) {
      dataX[ i ][ j ] = dataXY[ i ][ j ];
    }
    for( unsigned int j = 0; j < dimY; j++ ) {
      dataY[ i ][ j ] = dataXY[ i ][ dimX + j ];
    }
  }

}

#endif // QUESO_HAS_ANN
