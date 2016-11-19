#include <cmath>
#include <fstream>
#include <iomanip> // for setprecision

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <sensitivity_mc.h>

template<class P_V, class P_M, class Q_V, class Q_M>
Qoi_mc<P_V, P_M, Q_V, Q_M>::Qoi_mc(const char * prefix,
    const QUESO::VectorSet<P_V, P_M> & domainSet,
    const QUESO::VectorSet<Q_V, Q_M> & imageSet)
  : QUESO::BaseVectorFunction<P_V, P_M, Q_V, Q_M>(prefix, domainSet, imageSet),
    x_loc(3)
{
}

template<class P_V, class P_M, class Q_V, class Q_M>
Qoi_mc<P_V, P_M, Q_V, Q_M>::~Qoi_mc()
{
  // Deconstruct here
}

template<class P_V, class P_M, class Q_V, class Q_M>
void
Qoi_mc<P_V, P_M, Q_V, Q_M>::compute(const P_V & domainVector,
    const P_V * domainDirection,
    Q_V & imageVector, QUESO::DistArray<P_V *> * gradVectors,
    QUESO::DistArray<P_M *> * hessianMatrices,
    QUESO::DistArray<P_V *> * hessianEffects) const
{
  if (domainVector.sizeLocal() != 2) {
    queso_error_msg("domainVector does not have size 2");
  }
  if (imageVector.sizeLocal() != 1) {
    queso_error_msg("imageVector does not have size 1");
  }

  qoi_samples.open ("c_qoi_samplesAi.txt", std::fstream::in | std::fstream::out | std::fstream::app);

 // ----Generate Qoi using samples from a text files -------------
 count++;
 samples.open ("./files_sense/c_samples_Ai.txt", std::fstream::in | std::fstream::out | std::fstream::app);
 
 int cou = 1;

 while (samples >> mf >> cf){
       if (cou == count){
       		m = mf;
 		c = cf;
       	break;
       }
       cou++;
 }

// ------- Generate Qoi using pseudo-random MC samples ------------
//  std::cout << "m = " << domainVector[0] << std::endl; 
//  m = domainVector[0];  // Sample of the RV 'line slope'
//  c = domainVector[1];  // Sample of the RV 'y-intercept'
  double y_obs = 0.0;
  y_obs = m*x_loc + c;

  qoi_samples << std::setprecision(4) << y_obs << "\t\t" << m << "\t\t" << c << std::endl;

  imageVector[0] = y_obs;
  qoi_samples.close();
  samples.close();
}

template class Qoi_mc<QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector,
                   QUESO::GslMatrix>;
