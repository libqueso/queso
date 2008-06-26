#ifndef __UQ_GSL_ODE_SOLVER_H__
#define __UQ_GSL_ODE_SOLVER_H__

#include <uqGslVector.h>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

class uqGslOdeSolverClass
{
public:
  uqGslOdeSolverClass(const std::string& solverType,
                      unsigned int       sizeOfStateVector);
 ~uqGslOdeSolverClass();

  void solveODE(int (*stateDot   )(double t, const double currentState[], double answer[],             void* infoForComputingStateDot),
                int (*stateDotJac)(double t, const double currentState[], double* dfdy, double dfdt[], void* infoForComputingStateDot),
                const std::vector<double>&      instants,
                const uqGslVectorClass&         state0,
                void*                           infoForComputingStateDot,
                double                          suggestedDeltaT,
                std::vector<uqGslVectorClass*>& states);
protected:
  unsigned int               m_sizeOfStateVector;
  const gsl_odeiv_step_type* m_solver;
  gsl_odeiv_step*            m_s;
  gsl_odeiv_control*         m_c;
  gsl_odeiv_evolve*          m_e;
};

#endif // __UQ_GSL_ODE_SOLVER_H__
