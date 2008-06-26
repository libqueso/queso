#include <uqTrilinosParamSpace.h>

uqTrilinosParamSpaceClass::uqTrilinosParamSpaceClass(const Epetra_MpiComm& comm, unsigned int dimension)
  :
  uqParamSpaceClass<uqTrilinosVectorClass,uqTrilinosMatrixClass>(dimension),
  //m_rank(comm.MyPID()), // compiler complains FIXME
  m_map(new Epetra_Map(dimension,0,comm))
{
  //std::cout << "Entering uqTrilinosParamSpaceClass::constructor()"
  //          << std::endl;

  m_rank = comm.MyPID();

  //std::cout << "Leaving uqTrilinosParamSpaceClass::constructor()"
  //          << std::endl;
}
  
uqTrilinosParamSpaceClass::~uqTrilinosParamSpaceClass()
{
  std::cout << "Entering uqTrilinosParamSpaceClass::destructor()"
            << std::endl;

  if (m_map) delete m_map;

  std::cout << "Leaving uqTrilinosParamSpaceClass::destructor()"
            << std::endl;
}

uqTrilinosVectorClass*
uqTrilinosParamSpaceClass::newVector() const
{
  return new uqTrilinosVectorClass(*m_map);
}

uqTrilinosMatrixClass*
uqTrilinosParamSpaceClass::newMatrix() const
{
  return new uqTrilinosMatrixClass(*m_map,this->dim());
}

void
uqTrilinosParamSpaceClass::createInitialValues() const
{
  return;
}

void
uqTrilinosParamSpaceClass::createMinValues() const
{
  return;
}

void
uqTrilinosParamSpaceClass::createMaxValues() const
{
  return;
}

void
uqTrilinosParamSpaceClass::createPriorMuValues() const
{
  return;
}

void
uqTrilinosParamSpaceClass::createPriorSigmaValues() const
{
  return;
}

void
uqTrilinosParamSpaceClass::print(std::ostream& os) const
{
  return;
}

const Epetra_Map&
uqTrilinosParamSpaceClass::map() const
{
 return *m_map;
}

std::ostream&
operator<<(std::ostream& os, const uqTrilinosParamSpaceClass& space)
{
  space.uqParamSpaceClass<uqTrilinosVectorClass,uqTrilinosMatrixClass>::print(os);
  space.print(os);

  return os;
}
