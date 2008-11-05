/* uq/examples/fp/heat/uqHeatEx.C
 *
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

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
    ("uqHeatEx_help",                                                    "produce help message for UQ 1D FEM")
    ("uqHeatEx_c",      po::value<double      >()->default_value( 100.), "diffusion coefficient"             )
    ("uqHeatEx_f",      po::value<double      >()->default_value(1000.), "forcing term"                      )
    ("uqHeatEx_g",      po::value<double      >()->default_value(   0.), "left Dirichlet boundary condition" )
    ("uqHeatEx_h",      po::value<double      >()->default_value( 100.), "right Dirichlet boundary condition")
    ("uqHeatEx_numDiv", po::value<unsigned int>()->default_value(   20), "number of divisions for FEM"       )
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
  if (env.allOptionsMap().count("uqHeatEx_help")) {
    std::cout << optionsDesc
              << std::endl;
  }

  if (env.allOptionsMap().count("uqHeatEx_c")) {
    *c = env.allOptionsMap()["uqHeatEx_c"].as<double>();
  }

  if (env.allOptionsMap().count("uqHeatEx_f")) {
    *f = env.allOptionsMap()["uqHeatEx_f"].as<double>();
  }

  if (env.allOptionsMap().count("uqHeatEx_g")) {
    *g = env.allOptionsMap()["uqHeatEx_g"].as<double>();
  }

  if (env.allOptionsMap().count("uqHeatEx_h")) {
    *h = env.allOptionsMap()["uqHeatEx_h"].as<double>();
  }

  if (env.allOptionsMap().count("uqHeatEx_numDiv")) {
    *numDiv = env.allOptionsMap()["uqHeatEx_numDiv"].as<unsigned int>();
  }

  return;
}
