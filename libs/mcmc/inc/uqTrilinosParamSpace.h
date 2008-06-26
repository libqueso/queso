#ifndef __UQ_TRILINOS_PARAM_SPACE_H__
#define __UQ_TRILINOS_PARAM_SPACE_H__

#include <uqTrilinosVector.h>
#include <uqTrilinosMatrix.h>
#include <uqParamSpace.h>
#include <Epetra_MpiComm.h>

class uqTrilinosParamSpaceClass : public uqParamSpaceClass<uqTrilinosVectorClass,uqTrilinosMatrixClass>
{
public:
  uqTrilinosParamSpaceClass(const Epetra_MpiComm& comm, unsigned int dimension);
 ~uqTrilinosParamSpaceClass();

  uqTrilinosVectorClass* newVector             () const;
  uqTrilinosMatrixClass* newMatrix             () const;

  void                   print                 (std::ostream& os) const;

  const Epetra_Map&      map                   () const;

protected:
  void                   createInitialValues   () const;
  void                   createMinValues       () const;
  void                   createMaxValues       () const;
  void                   createPriorMuValues   () const;
  void                   createPriorSigmaValues() const;

private:
  const Epetra_Map* m_map;
};

std::ostream& operator<<(std::ostream& os, const uqTrilinosParamSpaceClass& space);

#endif // __UQ_TRILINOS_PARAM_SPACE_H__
