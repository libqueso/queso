#ifndef __UQ_ENVIRONMENT_H__
#define __UQ_ENVIRONMENT_H__

#undef UQ_USES_COMMAND_LINE_OPTIONS

#define UQ_ENV_VERBOSITY_DEFAULT_VALUE 0

#ifdef __UQ_USES_TRILINOS__
#include <Epetra_MpiComm.h>
#endif

#include <gsl/gsl_rng.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <iostream>
#include <fstream>

struct uqEnvOptionsStruct {
  uqEnvOptionsStruct(unsigned int verbosity);
 ~uqEnvOptionsStruct();

  unsigned int m_verbosity;
};

class uqEnvironmentClass {
public:
  uqEnvironmentClass();
  uqEnvironmentClass(int& argc, char** &argv);
  uqEnvironmentClass(const uqEnvOptionsStruct& options);
  uqEnvironmentClass(const uqEnvironmentClass& obj);
 ~uqEnvironmentClass();

  uqEnvironmentClass& operator= (const uqEnvironmentClass& rhs);

        int                      rank                     () const;
        void                     barrier                  () const;
#ifdef UQ_USES_COMMAND_LINE_OPTIONS
  const po::options_description& allOptionsDesc           () const;
#endif
        po::variables_map&       allOptionsMap            () const;
        void                     scanInputFileForMyOptions(const po::options_description& optionsDesc) const;
        unsigned int             verbosity                () const;
        gsl_rng*                 rng                      () const;
        bool                     isThereInputFile         () const;
        void                     print                    (std::ostream& os) const;

private:
        void                     commonConstructor        ();
        void                     readEventualInputFile    ();
        void                     defineMyOptions          (po::options_description& optionsDesc) const;
        void                     getMyOptionValues        (po::options_description& optionsDesc);

  int                      m_argc;
  char**                   m_argv;
#ifdef __UQ_USES_TRILINOS__
  Epetra_MpiComm*          m_comm;
#endif
  int                      m_rank;
  int                      m_commSize;
  bool                     m_argsWereProvided;
  bool                     m_thereIsInputFile;
  std::string              m_inputFileName;
  po::options_description* m_allOptionsDesc;
  po::options_description* m_envOptionsDesc;
  po::variables_map*       m_allOptionsMap;
  unsigned int             m_verbosity;
  gsl_rng*                 m_rng;
};

std::ostream& operator<<(std::ostream& os, const uqEnvironmentClass& obj);

#endif // __UQ_ENVIRONMENT_H__
