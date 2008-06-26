#ifndef __UQ_COV_COND_H__
#define __UQ_COV_COND_H__

#include <iostream>

template <class V, class M>
void
uqCovCond(
  double   condNumber,
  const V& direction,
  M&       covMatrix,
  M&       precMatrix)
{
  //std::cout << "Entering uqCovCond()"
  //          << std::endl;

  V v1(direction);
  unsigned int size1 = v1.size();
  //std::cout << "In uqCovCond(), v1 contents are:"
  //          << std::endl
  //          << v1
  //          << std::endl;

  V v2(direction.env(),condNumber,1.0,size1); // MATLAB linspace
  v2.cwInvert();
  v2.sort();
  //std::cout << "In uqCovCond(), v2 contents are:"
  //          << std::endl
  //          << v2
  //          << std::endl;

  double v1Norm2 = v1.norm2();
  if (v1[0] >=0) v1[0] += v1Norm2;
  else           v1[0] -= v1Norm2;
  double v1Norm2Sq = v1.norm2Sq();

  M Z(direction.env(),size1,1.0);
  Z -= (2./v1Norm2Sq) * matrixProduct(v1,v1);
  //std::cout << "In uqCovCond(), Z contents are:"
  //          << std::endl
  //          << Z
  //          << std::endl;

  M Zt(Z.transpose());
  covMatrix  = Z * diagScaling(v2,   Zt);
  precMatrix = Z * diagScaling(1./v2,Zt);

  //std::cout << "Leaving uqCovCond()"
  //          << std::endl;
}
#endif // __UQ_COV_COND_H__
