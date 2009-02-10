/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * $Id$
 *
 * Brief description of this file: 
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include <uqGslOdeSolver.h>
#include <uqDefines.h>

uqGslOdeSolverClass::uqGslOdeSolverClass(
  const std::string& solverType,
  unsigned int       sizeOfStateVector)
  :
  m_sizeOfStateVector(sizeOfStateVector)
{
  if (solverType == "45") {
    m_solver = gsl_odeiv_step_rkf45;
  }
  if (solverType == "gear2") {
    m_solver = gsl_odeiv_step_gear2;
  }
  else {
    m_solver = gsl_odeiv_step_rkf45;
  }
  m_s = gsl_odeiv_step_alloc   (m_solver,m_sizeOfStateVector);
  m_c = gsl_odeiv_control_y_new(1.e-6,0.);
  m_e = gsl_odeiv_evolve_alloc (m_sizeOfStateVector);
}

uqGslOdeSolverClass::~uqGslOdeSolverClass()
{
  gsl_odeiv_evolve_free (m_e);
  gsl_odeiv_control_free(m_c);
  gsl_odeiv_step_free   (m_s);
}

void
uqGslOdeSolverClass::solveODE(
  int (*stateDot   )(double t, const double currentState[], double answer[],             void* infoForComputingStateDot),
  int (*stateDotJac)(double t, const double currentState[], double* dfdy, double dfdt[], void* infoForComputingStateDot),
  const std::vector<double>&      instants,
  const uqGslVectorClass&         state0,
  void*                           infoForComputingStateDot,
  double                          suggestedDeltaInstant,
  std::vector<uqGslVectorClass*>& states)
{
  //std::cout << "Entering uqGslVectorClass::solveODE()"
  //          << std::endl;

  //if (state0.env().rank() == 0) {
  //  std::cout << "In uqGslVectorClass::solveODE()"
  //            << ": instants.size() = " << instants.size()
  //            << std::endl;
  //}

  int iRC;
  iRC = gsl_odeiv_evolve_reset(m_e);
  iRC = gsl_odeiv_step_reset  (m_s);

  UQ_FATAL_TEST_MACRO(state0.size() != m_sizeOfStateVector,
                      state0.env().rank(),
                      "uqGslVectorClass::solveODE()",
                      "state0 has size different than state size expected by this solver instance");

  UQ_FATAL_TEST_MACRO(instants[0] < 0.,
                      state0.env().rank(),
                      "uqGslVectorClass::solveODE()",
                      "routine does not accept negative instants");

  gsl_odeiv_system sys = {stateDot, stateDotJac, m_sizeOfStateVector, infoForComputingStateDot};

  double t            = 0.;
  double finalInstant = instants[instants.size()-1];
  double deltaInstant = suggestedDeltaInstant;
  double y[m_sizeOfStateVector];
  for (unsigned int j = 0; j < m_sizeOfStateVector; ++j) {
    y[j] = state0[j];
  }

  //**********************************************************
  // Reset the answer of this routine
  //**********************************************************
  for (unsigned int i = 0; i < states.size(); ++i) {
    if (states[i] != NULL) {
      delete states[i];
      states[i] = NULL;
    }
  }
  states.resize(instants.size(),NULL);

  //**********************************************************
  // Protection against roundoff errors when 't' is returned from 'gsl_odeiv_evolve_apply()'
  //**********************************************************
  bool nextInstantShouldNotPassARequestedInstant = false;

  //**********************************************************
  // Loop
  //**********************************************************
  unsigned int idOfRequestedInstant = 0;
  while (t < finalInstant) {
    //if (state0.env().rank() == 0) {
    //  std::cout << "In uqGslVectorClass::solveODE()"
    //            << ": beginning while loop"
    //            << ", idOfRequestedInstant = " << idOfRequestedInstant
    //            << ", requested instant = "    << instants[idOfRequestedInstant]
    //            << ", t = "                    << t
    //            << std::endl;
    //}
    //state0.env().barrier();

    bool bRC = nextInstantShouldNotPassARequestedInstant && (t > instants[idOfRequestedInstant]);
    if (bRC) {
      if (state0.env().rank() == 0) {
        std::cout << "In uqGslVectorClass::solveODE()"
                  << ": idOfRequestedInstant = " << idOfRequestedInstant
                  << ", requested instant = "    << instants[idOfRequestedInstant]
                  << ", t = "                    << t
                  << std::endl;
      }
      state0.env().barrier();
    }
    UQ_FATAL_TEST_MACRO(bRC,
                        state0.env().rank(),
                        "uqGslVectorClass::solveODE()",
                        "t should not pass current requested instant");
    nextInstantShouldNotPassARequestedInstant = false;

    //**********************************************************
    // Before going to the next instant,
    // check if a copy of the state is requested at this instant.
    //**********************************************************
    while (t == instants[idOfRequestedInstant]) {
      states[idOfRequestedInstant] = new uqGslVectorClass(state0);
      for (unsigned int j = 0; j < m_sizeOfStateVector; ++j) {
        (*states[idOfRequestedInstant])[j] = y[j];
      }
      idOfRequestedInstant++;

      //**********************************************************
      // State at instant 'finalInstant' should be recorded out of the while loop.
      // That is, 'idOfRequestedInstant' should value at most 'instants.size()-1' inside the while loop.
      //**********************************************************
      UQ_FATAL_TEST_MACRO(idOfRequestedInstant >= instants.size(),
                          state0.env().rank(),
                          "uqGslVectorClass::solveODE()",
                          "'idOfRequestedInstant' is too big inside the while loop");

      UQ_FATAL_TEST_MACRO(instants[idOfRequestedInstant] < instants[idOfRequestedInstant-1],
                          state0.env().rank(),
                          "uqGslVectorClass::solveODE()",
                          "routine does not accept decreasing instants");
    }

    //**********************************************************
    // Before going to the next instant,
    // check if a copy of the state is requested at an intermediate instant.
    // If so, adapt deltaInstant accordingly.
    // Use '<=' instead of just '<', for protection agains roundoffs.
    //**********************************************************
    if (instants[idOfRequestedInstant] <= (t+deltaInstant)) {
      deltaInstant = instants[idOfRequestedInstant]-t;
      nextInstantShouldNotPassARequestedInstant = true;
    }
    //else if (deltaInstant != suggestedDeltaInstant) {
    //  deltaInstant = suggestedDeltaInstant;
    //}
    
    //**********************************************************
    // Finally, go to the next instant.
    //**********************************************************
    //if (state0.env().rank() == 0) {
    //  std::cout << "In uqGslVectorClass::solveODE()"
    //            << ": just about to call gsl_odeiv_evolve_apply()"
    //            << ", idOfRequestedInstant = " << idOfRequestedInstant
    //            << ", requested instant = "    << instants[idOfRequestedInstant]
    //            << ", t = "                    << t
    //            << ", deltaInstant = "         << deltaInstant
    //            << std::endl;
    //}
    //state0.env().barrier();

    double currentInstant = t;
    iRC = gsl_odeiv_evolve_apply(m_e,m_c,m_s,&sys,&t,instants[idOfRequestedInstant],&deltaInstant,y);

    //if (state0.env().rank() == 0) {
    //  std::cout << "In uqGslVectorClass::solveODE()"
    //            << ": just returned from gsl_odeiv_evolve_apply()"
    //            << ", idOfRequestedInstant = " << idOfRequestedInstant
    //            << ", requested instant = "    << instants[idOfRequestedInstant]
    //            << ", t = "                    << t
    //            << ", deltaInstant = "         << deltaInstant
    //            << std::endl;
    //}
    //state0.env().barrier();

    UQ_FATAL_TEST_MACRO(iRC != GSL_SUCCESS,
                        state0.env().rank(),
                        "uqGslVectorClass::solveODE()",
                        "invalid return from gsl_odeiv_evolve_apply()");
    UQ_FATAL_TEST_MACRO(t <= currentInstant,
                        state0.env().rank(),
                        "uqGslVectorClass::solveODE()",
                        "t should always increase");
    UQ_FATAL_TEST_MACRO(deltaInstant <= 0.,
                        state0.env().rank(),
                        "uqGslVectorClass::solveODE()",
                        "deltaInstant should always be positive");
  }

  //**********************************************************
  // 'finalInstant' should be the only remaining requested instant.
  //**********************************************************
  UQ_FATAL_TEST_MACRO(finalInstant != instants[idOfRequestedInstant],
                      state0.env().rank(),
                      "uqGslVectorClass::solveODE()",
                      "'finalInstant' should be the only remainding requested instant after the while loop");
  UQ_FATAL_TEST_MACRO(nextInstantShouldNotPassARequestedInstant == false,
                      state0.env().rank(),
                      "uqGslVectorClass::solveODE()",
                      "'finalInstant' should be a requested instant at the end of the while loop");

  for (unsigned int i = idOfRequestedInstant; i < instants.size(); ++i) {
    states[i] = new uqGslVectorClass(state0);
    for (unsigned int j = 0; j < m_sizeOfStateVector; ++j) {
      (*states[i])[j] = y[j];
    }
  }

  //if (state0.env().rank() == 0) {
  //  std::cout << "In uqGslVectorClass::solveODE()" << std::endl;
  //  for (unsigned int i = 0; i < instants.size(); ++i) {
  //    std::cout << "  for t = " << instants[i]
  //              << ", state = " << *(states[i])
  //              << std::endl;
  //  }
  //}

  //std::cout << "Leaving uqGslVectorClass::solveODE()"
  //          << std::endl;

  return;
}

