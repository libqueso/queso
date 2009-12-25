/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * $Id$
 *
 * Brief description of this file: 
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include <example_compute.h>
#include <uqGslMatrix.h>
#include <uqVectorRV.h>

void compute(const uqFullEnvironmentClass& env) {
  // Step 1 of 9: Instantiate the parameter space
  uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>
    paramSpace(env, "param_", 2, NULL);

  // Step 2 of 9: Instantiate the parameter domain
  uqGslVectorClass paramMins(paramSpace.zeroVector());
  paramMins.cwSet(-INFINITY);
  uqGslVectorClass paramMaxs(paramSpace.zeroVector());
  paramMaxs.cwSet( INFINITY);
  uqBoxSubsetClass<uqGslVectorClass,uqGslMatrixClass>
    paramDomain("param_",paramSpace,paramMins,paramMaxs);

  // Step 3 of 9: Instantiate the vector RV
  uqGslVectorClass meanVector(paramSpace.zeroVector());
  meanVector[0] = -1;
  meanVector[1] =  2;
  uqGslMatrixClass* covMatrix = paramSpace.newMatrix();
  (*covMatrix)(0,0) = 4.; (*covMatrix)(0,1) = 0.;
  (*covMatrix)(1,0) = 0.; (*covMatrix)(1,1) = 7.;
  uqGaussianVectorRVClass<uqGslVectorClass,uqGslMatrixClass>
    auxRv("", paramDomain,meanVector,*covMatrix);

  // Step 4 of 9: Instantiate the vector sequence
  const int num_samples = 1000;

  uqSequenceOfVectorsClass<uqGslVectorClass,uqGslMatrixClass>
    auxSeq(paramSpace,num_samples,"aux_seq");

  //std::vector<double> convMeasure(num_samples,0.0);
  
  // Step 5 of 9: Populate the vector sequence
  uqGslVectorClass auxVec(paramSpace.zeroVector());
  for (unsigned int i = 0; i < auxSeq.subSequenceSize(); ++i) {
    auxRv.realizer().realization(auxVec);
    auxSeq.setPositionValues(i,auxVec);
//<<<<<<< .mine
    //if (( i >= 1 ) && (env.numSubEnvironments() > 1))
      //{
	//convMeasure[i] = auxSeq.estimateConvBrooksGelman( 0, i );
      //}
//=======
  //  if ((i >= 1) && (env.numSubEnvironments() > 1))
    //  {
	//convMeasure[i] = auxSeq.estimateConvBrooksGelman( 0, i );
      //}
//>>>>>>> .r5927
  }
  
//std::set<unsigned int> auxSet;
  //auxSet.insert(0);
  //auxSet.insert(1);
  //auxSeq.subWriteContents("anyname",auxSet);

  //if( env.inter0Rank() == 0 )
    //{
      //std::ofstream dataout( "convergence.m", std::ios::out );
      //dataout << "clear all" << std::endl;
      //dataout << "close all" << std::endl;
      //dataout << "index = [ " << 1 << std::endl;
      //for( int i = 2; i < num_samples-1; i++)
	//{
	  //dataout << i << std::endl;
	//}
      //dataout << num_samples-1 << " ];" << std::endl;

      //dataout << " data = [ " << convMeasure[1] << std::endl;
      //for( int i = 2; i < num_samples-1; i++)
	//{
	  //dataout << convMeasure[i] << std::endl;
	//}
      //dataout << convMeasure[ num_samples-1] << "];" << std::endl;
      //dataout << "plot( index, data, 'b-', 'LineWidth', 2 )" << std::endl;
      //dataout << "xlabel( 'Iteration', 'FontSize', 16 )" << std::endl;
      //dataout << "ylabel( 'BG-Convergence', 'FontSize', 16 )" << std::endl;
      //dataout << "title( 'Brooks-Gelman Convergence, Gaussian RV, " 
	//      << env.numSubEnvironments() << " Sequences', 'FontSize', 16 )" << std::endl;
      //dataout << "print -depsc BGConv" << env.numSubEnvironments() << ".eps" << std::endl;
      //std::cout <<"convMeasure = " << convMeasure[num_samples-1] << std::endl;
    //}

  // Step 6 of 9: Compute min, max, mean, covariance and correlation matrices
  uqGslVectorClass minVec (paramSpace.zeroVector());
  uqGslVectorClass maxVec (paramSpace.zeroVector());
  auxSeq.unifiedMinMax(0,auxSeq.subSequenceSize(),minVec,maxVec);

  uqGslVectorClass meanVec(paramSpace.zeroVector());
  auxSeq.unifiedMean(0,auxSeq.subSequenceSize(),meanVec);

  uqGslMatrixClass* covarianceMatrix  = paramSpace.newMatrix();
  uqGslMatrixClass* correlationMatrix = paramSpace.newMatrix();
  uqComputeCovCorrMatricesBetweenVectorSequences(auxSeq,
                                                 auxSeq,
                                                 auxSeq.subSequenceSize(),
                                                 *covarianceMatrix,
                                                 *correlationMatrix);
  
  // Step 7 of 9: Compute cdf accuracy
  unsigned int auxSize = 1000;
  //uqGslVectorClass deltaVec(maxVec-minVec);
  //deltaVec *= (1./(double) (auxSize-1));
  //std::vector<uqGslVectorClass*> evalPositionsVecs(auxSize,NULL);
  std::vector<uqGslVectorClass*> cdfStaccVecs(auxSize,NULL);
  std::vector<uqGslVectorClass*> cdfStaccVecsup(auxSize,NULL);
  std::vector<uqGslVectorClass*> cdfStaccVecslow(auxSize,NULL);
  std::vector<uqGslVectorClass*> xdataVecs(auxSize,NULL);
  std::ofstream outdata;
  std::ofstream outdata2;
  std::ofstream outdata3;
  std::ofstream outdata4;
  outdata.open("./cdf.dat");
  outdata2.open("./cdfup.dat");
  outdata3.open("./cdfdown.dat");
  outdata4.open("./xdata.dat");
  auxSeq.subCdfStacc(0,cdfStaccVecs,cdfStaccVecsup,cdfStaccVecslow,xdataVecs);
//if (env.fullRank() == 0) {
  for (unsigned int i = 0; i < 1000; ++i) {
    //evalPositionsVecs[i] = new uqGslVectorClass(paramSpace.zeroVector());
    //*(evalPositionsVecs[i]) = minVec + ((double) i)*deltaVec;
        //std::cout <<*cdfStaccVecs[i]<< std::endl;
	//std::cout <<*cdfStaccVecsup[i]<< std::endl;
	//std::cout <<*cdfStaccVecslow[i]<< std::endl;
        //std::cout <<*xdataVecs[i]<< std::endl; 
   outdata <<*cdfStaccVecs[i]<< std::endl;
   outdata2 <<*cdfStaccVecsup[i]<< std::endl;
   outdata3 <<*cdfStaccVecslow[i]<< std::endl;
   outdata4 <<*xdataVecs[i]<< std::endl;
    // std::cout <<auxVec[i]<<std::endl;
    //cdfStaccVecs     [i] = new uqGslVectorClass(paramSpace.zeroVector());
  }
 outdata.close();
 outdata2.close();
 outdata3.close();
 outdata4.close();
//}
  //uqGslVectorClass cdfStaccVecs(paramSpace.zeroVector());
  
//if (env.fullRank() == 0) {
  //  std::cout << "\n minVec = "  << minVec
    //          << "\n maxVec = "  << maxVec
      //        << "\n meanVec = " << meanVec
        //      << "\n covMat = "  << *covarianceMatrix
          //    << "\n corrMat = " << *correlationMatrix
            //  << std::endl;
  //}

 //std::cout << "\n cdfmean = "  << cdfStaccVecs<< std::endl;
  
  // Return
  
for (unsigned int i = 0; i < auxSize; ++i) {
  //  delete evalPositionsVecs[i];
    delete cdfStaccVecs     [i];
  }
  delete correlationMatrix;
  delete covarianceMatrix;
  delete covMatrix;

  return;
}
