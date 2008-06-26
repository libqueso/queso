#include <uq1DFemDakotaHandler.h>
#include <iostream>
#include <fstream>

uq1DFemDakotaHandlerClass::uq1DFemDakotaHandlerClass(
  uq1DFemDakotaRequestTypeEnum drt,
  std::string&                 inputFileName)
  :
  m_drt                   (drt),
  m_c                     (0.),
  m_f                     (0.),
  m_g                     (0.),
  m_h                     (0.),
  m_numDiv                (0),
  m_xValues               (0),
  m_uValues               (0),
  m_asv                   (0)
{
  int iRC;

  if (drt == UQ_1D_FEM_FP_DAKOTA_REQUEST_TYPE) {
    iRC = fpRead(inputFileName);
    if (iRC) {
      std::cerr << "Wrong dakota fp information passef from Dakota to " << "uq1DFem_gsl"
                << std::endl;
      exit(-1);
    }
  }
  else if (drt == UQ_1D_FEM_PE_DAKOTA_REQUEST_TYPE) {
    iRC = peRead(inputFileName);
    if (iRC) {
      std::cerr << "Wrong dakota pe information passef from Dakota to " << "uq1DFem_gsl"
                << std::endl;
      exit(-1);
    }
  }
  else {
    std::cerr << "Wrong dakota request (" << drt
              << ") inside uq1DFemDakotaHandlerClass::constructor()"
              << std::endl;
    exit(-1);
  }

  if (m_asv.size() < 1) {
    std::cerr << "Wrong number ("                 << m_asv.size()
              << ") of functions from Dakota to " << "uq1DFem_gsl"
              << std::endl;
    exit(-1);
  }
}

uq1DFemDakotaHandlerClass::~uq1DFemDakotaHandlerClass()
{
}

int
uq1DFemDakotaHandlerClass::fpRead(std::string& inputFileName)
{
  int iRC = commonRead(inputFileName);
  if (iRC) return iRC;

  if (m_xValues.size() != m_asv.size()) {
    std::cerr << "Number of x ("                             << m_xValues.size()
              << ") values and of functions ("               << m_asv.size()
              << ") to be evaluated, passed from Dakota to " << "uq1DFem_gsl"
              << ", are different"
              << std::endl;
    return -1;
  }

  return 0;
}

int
uq1DFemDakotaHandlerClass::peRead(std::string& inputFileName)
{
  int iRC = commonRead(inputFileName);
  if (iRC) return iRC;

  if (m_xValues.size() != m_uValues.size()) {
    std::cerr << "Number of x ("                    << m_xValues.size()
              << ") values and of u ("              << m_uValues.size()
              << ") values, passed from Dakota to " << "uq1DFem_gsl"
              << ", are different"
              << std::endl;
    return -1;
  }

  if (m_asv.size() != 1) {
    std::cerr << "Number of functions ("     << m_asv.size()
              << "), passed from Dakota to " << "uq1DFem_gsl"
              << ", is invalid"
              << std::endl;
    return -1;
  }

  return 0;
}

int
uq1DFemDakotaHandlerClass::commonRead(std::string& inputFileName)
{
  std::ifstream ifs(inputFileName.c_str());
  if (!ifs) {
    std::cerr << "Wrong input file name ('" << inputFileName
              << "') from Dakota to "       << "uq1DFem_gsl"
              << std::endl;
    return -1;
  }

  // Get the parameter vector
  unsigned int num_vars;
  std::string vars_text;
  ifs >> num_vars >> vars_text;

  if (num_vars < 5) {
    std::cerr << "Wrong number ("                  << num_vars
              << ") of parameters from Dakota to " << "uq1DFem_gsl"
              << std::endl;
    return -1;
  }

  std::vector<double>      paramValues(num_vars);
  std::vector<std::string> paramNames (num_vars);
  for (unsigned int i = 0; i < paramValues.size(); ++i) {
    ifs >> paramValues[i] >> paramNames[i];
  }

  unsigned int num_xValues = 0;
  unsigned int num_uValues = 0;
  for (unsigned int i = 0; i < paramValues.size(); ++i) {
    if      (strncmp(paramNames[i].c_str(),"X",1) == 0) num_xValues++;
    else if (strncmp(paramNames[i].c_str(),"U",1) == 0) num_uValues++;
  }
  m_xValues.resize(num_xValues,0.);
  m_uValues.resize(num_uValues,0.);

  for (unsigned int i = 0; i < paramValues.size(); ++i) {
    if (strcmp(paramNames[i].c_str(),"cTerm") == 0) {
      m_c = paramValues[i];
    }
    else if (strcmp(paramNames[i].c_str(),"fTerm") == 0) {
      m_f = paramValues[i];
    }
    else if (strcmp(paramNames[i].c_str(),"gTerm") == 0) {
      m_g = paramValues[i];
    }
    else if (strcmp(paramNames[i].c_str(),"hTerm") == 0) {
      m_h = paramValues[i];
    }
    else if (strcmp(paramNames[i].c_str(),"numDiv") == 0) {
      m_numDiv = (unsigned int) paramValues[i];
    }
    else if (strncmp(paramNames[i].c_str(),"X",1) == 0) {
      const char* justDigits = paramNames[i].c_str();
      double index = strtod(&justDigits[1],NULL);
      //std::cout << "In uq1DFemDakotaHandlerClass::commonRead()"
      //          << ", param = " << paramNames[i]
      //          << ", index = " << index
      //          << std::endl;
      m_xValues[(unsigned int) index] = paramValues[i];
    }
    else if (strncmp(paramNames[i].c_str(),"U",1) == 0) {
      const char* justDigits = paramNames[i].c_str();
      double index = strtod(&justDigits[1],NULL);
      //std::cout << "In uq1DFemDakotaHandlerClass::commonRead()"
      //          << ", param = " << paramNames[i]
      //          << ", index = " << index
      //          << std::endl;
      m_uValues[(unsigned int) index] = paramValues[i];
    }
    else {
      std::cerr << "Wrong name ('"                   << paramNames[i]
                << "') of parameter from Dakota to " << "uq1DFem_gsl"
                << std::endl;
      return -1;
    }
  }

  std::cout << "\nAfter reading from Dakota:";
  std::cout << "\n m_c = " << m_c;
  std::cout << "\n m_f = " << m_f;
  std::cout << "\n m_g = " << m_g;
  std::cout << "\n m_h = " << m_h;
  std::cout << "\n m_numDiv = " << m_numDiv;
  std::cout << "\nm_xValues =";
  for (unsigned int i = 0; i < m_xValues.size(); ++i) {
    std::cout << " " << m_xValues[i];
  }
  std::cout << std::endl;
  std::cout << "\nm_uValues =";
  for (unsigned int i = 0; i < m_uValues.size(); ++i) {
    std::cout << " " << m_uValues[i];
  }
  std::cout << std::endl;

  // Get the ASV vector and ignore the labels
  unsigned int num_fns;
  std::string fns_text;
  ifs >> num_fns >> fns_text;

  m_asv.resize(num_fns,0);
  for (unsigned int i = 0; i < m_asv.size(); ++i) {
    ifs >> m_asv[i];
    ifs.ignore(256, '\n');
  }

  return 0;
}

void
uq1DFemDakotaHandlerClass::writeResults(
  std::string&            outputFileName,
  const uq1DProblemClass& problem)
{
  std::ofstream ofs(outputFileName.c_str());

  if (m_drt == UQ_1D_FEM_FP_DAKOTA_REQUEST_TYPE) {
    for (unsigned int i = 0; i < m_asv.size(); ++i) {
      if (m_asv[i] & 1) {
        ofs << " " << problem.uValue(m_xValues[i]) << std::endl;
      }
      if (m_asv[i] & 2) {
        std::cerr << "In case of a forward propagation request from Dakota"
                  << ", uq1DFemDakotaHandlerClass::writeResults()"
                  << " does not return derivatives"
                  << std::endl;
        exit(-1);
      }
      if (m_asv[i] & 4) {
        std::cerr << "In case of a forward propagation request from Dakota"
                  << ", uq1DFemDakotaHandlerClass::writeResults()"
                  << " does not return Hessians"
                  << std::endl;
        exit(-1);
      }
    }
    std::cout << "\nInfo for checking data written by uq1DFemDakotaHandlerClass::writeResults():";
    std::cout << "\nm_xValues =";
    for (unsigned int i = 0; i < m_xValues.size(); ++i) {
      std::cout << " " << m_xValues[i];
    }
    std::cout << "\nm_uExacts =";
    for (unsigned int i = 0; i < m_xValues.size(); ++i) {
      std::cout << " " << problem.uExact(m_xValues[i]);
    }
  }
  else if (m_drt == UQ_1D_FEM_PE_DAKOTA_REQUEST_TYPE) {
    double cumulativeSquaredError = 0.;
    for (unsigned int i = 0; i < m_xValues.size(); ++i) {
      double diff = m_uValues[i] - problem.uValue(m_xValues[i]);
      cumulativeSquaredError += diff*diff;
    }
    for (unsigned int i = 0; i < m_asv.size(); ++i) {
      if (m_asv[i] & 1) {
        ofs << " " << cumulativeSquaredError << std::endl;
      }
      if (m_asv[i] & 2) {
        std::cerr << "In case of a parameter estimation request from Dakota"
                  << ", uq1DFemDakotaHandlerClass::writeResults()"
                  << " does not return derivatives"
                  << std::endl;
        exit(-1);
      }
      if (m_asv[i] & 4) {
        std::cerr << "In case of a parameter estimation request from Dakota"
                  << ", uq1DFemDakotaHandlerClass::writeResults()"
                  << " does not return Hessians"
                  << std::endl;
        exit(-1);
      }
    }
  }
  else {
    std::cerr << "Wrong dakota request (" << m_drt
              << ") inside uq1DFemDakotaHandlerClass::writeResults()"
              << std::endl;
    exit(-1);
  }

  return;
}

double
uq1DFemDakotaHandlerClass::c() const
{
  return m_c;
}

double
uq1DFemDakotaHandlerClass::f() const
{
  return m_f;
}

double
uq1DFemDakotaHandlerClass::g() const
{
  return m_g;
}

double
uq1DFemDakotaHandlerClass::h() const
{
  return m_h;
}
