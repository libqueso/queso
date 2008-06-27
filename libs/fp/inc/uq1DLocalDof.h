#ifndef __UQ_1D_LOCAL_DOFS_H__
#define __UQ_1D_LOCAL_DOFS_H__

#include <iostream>
#include <uqDefines.h>

enum uqLocalDofPositionEnum {
  UQ_INTERNAL_LOCAL_DOF_POS = 0,
  UQ_AT_NODE_LOCAL_DOF_POS,
  UQ_AT_EDGE_LOCAL_DOF_POS,
  UQ_AT_FACE_LOCAL_DOF_POS
};

class uq1DLocalDofClass
{
public:
  uq1DLocalDofClass(double                 bcc,
                    uqLocalDofPositionEnum localPositionType,
                    unsigned int           globalIdOfRespectiveNode);
  uq1DLocalDofClass(const uq1DLocalDofClass& obj);
 ~uq1DLocalDofClass();

  uq1DLocalDofClass& operator=(const uq1DLocalDofClass& rhs);

  double       bcc                     () const;
  unsigned int globalId                () const;
  void         setGlobalId             (unsigned int globalId);
  unsigned int globalIdOfRespectiveNode() const;
  void         print                   (std::ostream& os) const;

protected:
  void         copy                    (const uq1DLocalDofClass& src);

  double                 m_bcc;
  uqLocalDofPositionEnum m_myLocalPositionType;
  unsigned int           m_globalIdOfRespectiveNode;
  unsigned int           m_myGlobalId;
};

#endif // __UQ_1D_LOCAL_DOFS_H__
