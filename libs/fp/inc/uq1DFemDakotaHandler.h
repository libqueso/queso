#ifndef __UQ_1D_FEM_DAKOTA_HANDLER__
#define __UQ_1D_FEM_DAKOTA_HANDLER__

#include <uq1DProblem.h>
#include <vector>
#include <iostream>
#include <string>

#define UQ_1D_FEM_CMD_LINE_SHIFT_VALUE 1

enum uq1DFemDakotaRequestTypeEnum {
  UQ_1D_FEM_FP_DAKOTA_REQUEST_TYPE = 0, // forward propagation
  UQ_1D_FEM_PE_DAKOTA_REQUEST_TYPE      // parameter estimation
};

class uq1DFemDakotaHandlerClass
{
public:
  uq1DFemDakotaHandlerClass(uq1DFemDakotaRequestTypeEnum drt,
                            std::string& inputFileName);
 ~uq1DFemDakotaHandlerClass();

  void writeResults(std::string&            outputFileName,
                    const uq1DProblemClass& problem);
  double c() const;
  double f() const;
  double g() const;
  double h() const;

protected:
  int fpRead    (std::string& inputFileName);
  int peRead    (std::string& inputFileName);
  int commonRead(std::string& inputFileName);

  uq1DFemDakotaRequestTypeEnum m_drt;
  double                       m_c;
  double                       m_f;
  double                       m_g;
  double                       m_h;
  unsigned int                 m_numDiv;
  std::vector<double>          m_xValues;
  std::vector<double>          m_uValues;
  std::vector<int>             m_asv;
};

#endif // __UQ_1D_FEM_DAKOTA_HANDLER__
