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
// $Id: uqInfoTheory.h 14612 2011-01-02 12:39:58Z gabriel $
//
//--------------------------------------------------------------------------

#ifndef __UQ_INFO_THEORY_H__
#define __UQ_INFO_THEORY_H__

#include <uqDefines.h>
#ifdef QUESO_HAS_ANN

#include <ANN/ANN.h>
#include <ANN/ANNx.h>

#define UQ_INFTH_ANN_NO_SMP        10000
#define UQ_INFTH_ANN_EPS           0.0
#define UQ_INFTH_ANN_KNN           6


void distANN_XY( const ANNpointArray dataX, const ANNpointArray dataY, 
		 double* distsXY, 
		 unsigned int dimX, unsigned int dimY, 
		 unsigned int xN, unsigned int yN, 
		 unsigned int k, double eps );

#endif // QUESO_HAS_ANN

#endif // __UQ_INFO_THEORY_H__
