#include <uqEnvironment.h>
#include <uq1DFemDakotaHandler.h>
#include <string>

void defineMyOptions  (po::options_description& optionsDesc);
void getMyOptionValues(const uqEnvironmentClass& env,
                       po::options_description&  optionsDesc,
                       double*                   c,
                       double*                   f,
                       double*                   g,
                       double*                   h,
                       unsigned int*             numDiv);

int main(int argc, char* argv[])
{
  double       a = 0.;
  double       b = 1.;
  unsigned int order = 2;

  double       c = 100.;
  double       f = 1000.;
  double       g = 0.;
  double       h = 100.;
  unsigned int numDiv = 20;

  uqEnvironmentClass env;
  po::options_description optionsDesc("1D FEM problem options");

  //////////////////////////////////////////////////
  // Read data eventually passed through input files
  //////////////////////////////////////////////////
  uq1DFemDakotaHandlerClass* dakotaHandler = NULL;
  if ((argc                               == 4) &&
      (strcmp(argv[1],"--fpDakotaDriven") == 0)) {
    std::string inputFileName(argv[1+UQ_1D_FEM_CMD_LINE_SHIFT_VALUE]);
    dakotaHandler = new uq1DFemDakotaHandlerClass(UQ_1D_FEM_FP_DAKOTA_REQUEST_TYPE, inputFileName);
  }
  else if ((argc                               == 4) &&
           (strcmp(argv[1],"--peDakotaDriven") == 0)) {
    std::string inputFileName(argv[1+UQ_1D_FEM_CMD_LINE_SHIFT_VALUE]);
    dakotaHandler = new uq1DFemDakotaHandlerClass(UQ_1D_FEM_PE_DAKOTA_REQUEST_TYPE, inputFileName);
  }

  if (dakotaHandler) {
    c = dakotaHandler->c();
    f = dakotaHandler->f();
    g = dakotaHandler->g();
    h = dakotaHandler->h();
  }
  else {
    defineMyOptions              (optionsDesc);
    env.scanInputFileForMyOptions(optionsDesc);
    getMyOptionValues            (env,optionsDesc,&c,&f,&g,&h,&numDiv);
  }

  //////////////////////////////////////////////////
  // Define problem
  //////////////////////////////////////////////////
  uq1DProblemClass problem(env,
                           a,
                           b,
                           c,
                           f,
                           g,
                           h);

  //////////////////////////////////////////////////
  // Solve problem with FEM  
  //////////////////////////////////////////////////
  problem.femSolve(numDiv, order);

  //////////////////////////////////////////////////
  // Dump results to output stream
  //////////////////////////////////////////////////
  if (dakotaHandler) {
    std::string outputFileName(argv[2+UQ_1D_FEM_CMD_LINE_SHIFT_VALUE]);
    dakotaHandler->writeResults(outputFileName, problem);
  }
  else {
    double uValueAtLeft   = problem.uValue(a);
    double uValueAtCenter = problem.uValue(.5*(a+b));
    double uValueAtRight  = problem.uValue(b);
    std::cout << "values = " << uValueAtLeft
              << ", "        << uValueAtCenter
              << ", "        << uValueAtRight
              << std::endl;
  }

  return 0;
}

void defineMyOptions(po::options_description& optionsDesc)
{
  optionsDesc.add_options()
    ("uq1DFem_help",                                                    "produce help message for UQ 1D FEM")
    ("uq1DFem_c",      po::value<double      >()->default_value( 100.), "diffusion coefficient"             )
    ("uq1DFem_f",      po::value<double      >()->default_value(1000.), "forcing term"                      )
    ("uq1DFem_g",      po::value<double      >()->default_value(   0.), "left Dirichlet boundary condition" )
    ("uq1DFem_h",      po::value<double      >()->default_value( 100.), "right Dirichlet boundary condition")
    ("uq1DFem_numDiv", po::value<unsigned int>()->default_value(   20), "number of divisions for FEM"       )
  ;

  return;
}

void getMyOptionValues(
  const uqEnvironmentClass& env,
  po::options_description&  optionsDesc,
  double*                   c,
  double*                   f,
  double*                   g,
  double*                   h,
  unsigned int*             numDiv)
{
  if (env.allOptionsMap().count("uq1DFem_help")) {
    std::cout << optionsDesc
              << std::endl;
  }

  if (env.allOptionsMap().count("uq1DFem_c")) {
    *c = env.allOptionsMap()["uq1DFem_c"].as<double>();
  }

  if (env.allOptionsMap().count("uq1DFem_f")) {
    *f = env.allOptionsMap()["uq1DFem_f"].as<double>();
  }

  if (env.allOptionsMap().count("uq1DFem_g")) {
    *g = env.allOptionsMap()["uq1DFem_g"].as<double>();
  }

  if (env.allOptionsMap().count("uq1DFem_h")) {
    *h = env.allOptionsMap()["uq1DFem_h"].as<double>();
  }

  if (env.allOptionsMap().count("uq1DFem_numDiv")) {
    *numDiv = env.allOptionsMap()["uq1DFem_numDiv"].as<unsigned int>();
  }

  return;
}
