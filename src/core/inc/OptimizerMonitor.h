//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2017 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

#ifndef UQ_OPTIMIZER_MONITOR_H
#define UQ_OPTIMIZER_MONITOR_H

#include <vector>
#include <iostream>

namespace QUESO
{
  class BaseEnvironment;

  //! Object to monitor convergence of optimizers
  class OptimizerMonitor
  {
  public:

    //! Constructor
    /*!
     * Input paramter n_iters will reserve space for monitor data
     * for n_iters. This is for efficiency - if the number of iterations
     * goes beyond n_iters, there will not be an error, just less performance
     * by the monitor. Defaults to reserving space for 100 iterations
     */
    OptimizerMonitor( const BaseEnvironment& env, unsigned int n_iters = 100 );

    ~OptimizerMonitor();

    //! Enabling output to std::cout everytime append is called
    /*!
     * Helpful when wanting to monitor convergence progress in real time.
     * The print_xmin argument controls whether or not the current estimate
     * of the minimizer is printed. Not recommended if you have more than a
     * handful of dimensions over which you are optimizing.
     */
    void set_display_output( bool enable_output, bool print_xmin );

    //! Add current value of minimizer, objective, and error norm
    void append( std::vector<double>& x_min, double objective, double norm );

    //! Clears internal datastructures and resets display variables to false.
    void reset();

    //! Prints entire convergence history.
    /*!
     *  The print_xmin argument controls whether or not the estimate
     *  of the minimizer (all components) is printed.
     */
    void print( std::ostream& output, bool print_xmin = false ) const;

  private:

    const BaseEnvironment& m_env;

    bool m_display_conv;
    bool m_print_xmin;

    std::vector<std::vector<double> > m_minimizer_hist;
    std::vector<double> m_objective_hist;
    std::vector<double> m_norm_hist;

    void print_header( std::ostream& output, bool print_xmin ) const;
    void print_iteration( unsigned int iter, std::ostream& output,
                          bool print_xmin ) const;

  };

} // end namespace QUESO

#endif // UQ_OPTIMIZER_MONITOR_H
