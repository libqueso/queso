#ifndef __UQ_1D_NODE_H__
#define __UQ_1D_NODE_H__

#include <uqEnvironment.h>
#include <iostream>

enum uqNodePositionEnum {
  UQ_INTERNAL_NODE_POS = 0,
  UQ_AT_BOUNDARY_NODE_NODE_POS,
  UQ_AT_BOUNDARY_EDGE_NODE_POS,
  UQ_AT_BOUNDARY_FACE_NODE_POS
};

enum uqBCEnum {
  UQ_NONE_BC = 0,
  UQ_DIRICHLET_BC
};

class uq1DNodeClass
{
public:
  uq1DNodeClass(unsigned int       globalId,
                double             x,
                uqNodePositionEnum nodePositionType,
                uqBCEnum           bcType,
                double             bcValue);
  uq1DNodeClass(const uq1DNodeClass& obj);
 ~uq1DNodeClass();

  uq1DNodeClass& operator=(const uq1DNodeClass& rhs);

  unsigned int globalId() const;
  double       x() const;
  void         print(std::ostream& os) const;

protected:
  void copy(const uq1DNodeClass& src);

  unsigned int       m_myGlobalId;
  double             m_myX;
  uqNodePositionEnum m_myPositionType;
  uqBCEnum           m_myBCType;
  double             m_myBCValue;
};

#endif // __UQ_1D_NODE_H__
