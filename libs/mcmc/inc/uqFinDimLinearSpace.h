#ifndef __UQ_FIN_DIM_LINEAR_SPACE_H__
#define __UQ_FIN_DIM_LINEAR_SPACE_H__

#include <uqEnvironment.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <uqDefines.h>

template <class V, class M>
class uqFinDimLinearSpaceClass
{
public:
           uqFinDimLinearSpaceClass();
           uqFinDimLinearSpaceClass(const uqEnvironmentClass& env);
  virtual ~uqFinDimLinearSpaceClass();

  virtual unsigned int      dim                       ()                 const = 0;
          V*                newVector                 ()                 const; // See template specialization
          V*                newVector                 (const V& v)       const;
          M*                newMatrix                 ()                 const; // See template specialization
          M*                newDiagMatrix             (const V& v)       const;
          M*                newDiagMatrix             (double diagValue) const; // See template specialization

  virtual void              print                     (std::ostream& os) const;

protected:
  const uqEnvironmentClass& m_env;
  unsigned int              m_dim;
};

template <class V, class M>
uqFinDimLinearSpaceClass<V,M>::uqFinDimLinearSpaceClass()
  :
  m_env(*(new uqEnvironmentClass()))
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqFinDimLinearSpaceClass<V,M>::constructor(), default",
                      "should not be used by user");
}

template <class V, class M>
uqFinDimLinearSpaceClass<V,M>::uqFinDimLinearSpaceClass(
  const uqEnvironmentClass& env)
  :
  m_env(env),
  m_dim(0)
{
  //std::cout << "Entering uqFinDimLinearSpaceClass<V,M>::constructor()"
  //          << std::endl;

  //std::cout << "Leaving uqFinDimLinearSpaceClass<V,M>::constructor()"
  //          << std::endl;
}

template <class V, class M>
uqFinDimLinearSpaceClass<V,M>::~uqFinDimLinearSpaceClass()
{
  //std::cout << "Entering uqFinDimLinearSpaceClass<V,M>::destructor()"
  //          << std::endl;

  //std::cout << "Leaving uqFinDimLinearSpaceClass<V,M>::destructor()"
  //          << std::endl;
}

template <class V, class M>
unsigned int
uqFinDimLinearSpaceClass<V,M>::dim() const
{
  return m_dim;
}

template <class V, class M>
V*
uqFinDimLinearSpaceClass<V,M>::newVector(const V& v) const
{
  if (v.size() != m_dim) return NULL;

  return new V(v);
}

template <class V, class M>
M*
  uqFinDimLinearSpaceClass<V,M>::newDiagMatrix(const V& v) const
{
  if (v.size() != m_dim) return NULL;

  return new M(v);
}

template <class V, class M>
void
uqFinDimLinearSpaceClass<V,M>::print(std::ostream& os) const
{
  return;
}
#endif // __UQ_FIN_DIM_LINEAR_SPACE_H__

