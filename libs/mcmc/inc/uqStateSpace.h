#ifndef __UQ_STATE_SPACE_H__
#define __UQ_STATE_SPACE_H__

#include <uqFinDimLinearSpace.h>

template <class V, class M>
class uqStateSpaceClass : public uqFinDimLinearSpaceClass<V,M>
{
public:
           uqStateSpaceClass();
           uqStateSpaceClass(const uqEnvironmentClass& env);
  virtual ~uqStateSpaceClass();

          unsigned int            dim                       () const;
  virtual void                    print                     (std::ostream& os) const;

protected:
          void                    defineMyOptions           (po::options_description& optionsDesc) const;
          void                    getMyOptionValues         (po::options_description& optionsDesc);

  po::options_description* m_optionsDesc;

  using uqFinDimLinearSpaceClass<V,M>::m_env;
  using uqFinDimLinearSpaceClass<V,M>::m_dim;
};

template <class V, class M>
uqStateSpaceClass<V,M>::uqStateSpaceClass()
  :
  uqFinDimLinearSpaceClass<V,M>()
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqStateSpaceClass<V,M>::constructor(), default",
                      "should not be used by user");
}

template <class V, class M>
uqStateSpaceClass<V,M>::uqStateSpaceClass(
  const uqEnvironmentClass& env)
  :
  uqFinDimLinearSpaceClass<V,M>(env),
  m_optionsDesc                (NULL)
{
  //std::cout << "Entering uqStateSpaceClass<V,M>::constructor()"
  //          << std::endl;

  m_optionsDesc = new po::options_description("State space options");
  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.rank() == 0) std::cout << "After getting option values, state of uqStateSpaceClass object is:"
                                   << "\ndimension = " << m_dim
                                   << "\n"
                                   << std::endl;

  //std::cout << "Leaving uqStateSpaceClass<V,M>::constructor()"
  //          << std::endl;
}

template <class V, class M>
uqStateSpaceClass<V,M>::~uqStateSpaceClass()
{
  //std::cout << "Entering uqStateSpaceClass<V,M>::destructor()"
  //          << std::endl;

  if (m_optionsDesc) delete m_optionsDesc;

  //std::cout << "Leaving uqStateSpaceClass<V,M>::destructor()"
  //          << std::endl;
}

template <class V, class M>
void
uqStateSpaceClass<V,M>::defineMyOptions(
  po::options_description& optionsDesc) const
{
  m_optionsDesc->add_options()
    ("uqStateSpace_help",                                              "produce help message for UQ PS")
    ("uqStateSpace_dim",  po::value<unsigned int>()->default_value(0), "Space dimension"               )
  ;

  return;
}

template <class V, class M>
void
uqStateSpaceClass<V,M>::getMyOptionValues(po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count("uqStateSpace_help")) {
    std::cout << optionsDesc
              << std::endl;
  }

  if (m_env.allOptionsMap().count("uqStateSpace_dim")) {
    const po::variables_map& tmpMap = m_env.allOptionsMap();
    m_dim = tmpMap["uqStateSpace_dim"].as<unsigned int>();
  }

  return;
}

template <class V, class M>
unsigned int
uqStateSpaceClass<V,M>::dim() const
{
  return m_dim;
}

template <class V, class M>
void
uqStateSpaceClass<V,M>::print(std::ostream& os) const
{
  os << "m_dim = " << m_dim
     << std::endl;

  return;
}

template<class V, class M>
std::ostream&
operator<<(std::ostream& os, const uqStateSpaceClass<V,M>& space)
{
  space.print(os);

  return os;
}
#endif // __UQ_STATE_SPACE_H__

