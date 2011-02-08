//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010 The PECOS Development Team
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

#include <config.h>
#ifdef HAVE_ANN

#include <uqInfoTheory.h>

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

#endif // HAVE_ANN
