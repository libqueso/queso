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

#include <string>
#include <iomanip>
#include <queso/asserts.h>
#include <queso/Environment.h>
#include <queso/OptimizerMonitor.h>

namespace QUESO
{

  OptimizerMonitor::OptimizerMonitor( const BaseEnvironment& env, unsigned int n_iters )
    : m_env(env),
      m_display_conv(false),
      m_print_xmin(false)
  {
    m_minimizer_hist.reserve(n_iters);
    m_objective_hist.reserve(n_iters);
    m_norm_hist.reserve(n_iters);
  }

  OptimizerMonitor::~OptimizerMonitor()
  {}

  void OptimizerMonitor::set_display_output( bool enable_output, bool print_xmin )
  {
    m_display_conv = enable_output;
    m_print_xmin = print_xmin;
  }

  void OptimizerMonitor::append( std::vector<double>& x_min, double objective, double norm )
  {
    // This needs to be done before trying to print
    m_minimizer_hist.push_back(x_min);
    m_objective_hist.push_back(objective);
    m_norm_hist.push_back(norm);

    // Print out to screen if the user set the option
    if( m_display_conv )
      {
        // If we are appending the first entry, print a nice header
        if( m_minimizer_hist.size() == 1 )
          {
            this->print_header(std::cout, m_print_xmin);
          }

        /* We're assuming here the size of the array is the current iteration.
           Shift by -1 since we count from 0 in C++. */
        this->print_iteration(m_norm_hist.size()-1,std::cout,m_print_xmin);
      }

  }

  void OptimizerMonitor::print_header( std::ostream& output, bool print_xmin ) const
  {
    unsigned int width = 37;

    if( print_xmin) width += (m_minimizer_hist[0]).size()*15;

    // Only print output on processor 0
    if( m_env.fullRank() == 0 )
      {
        output.width(5);
        output << "i";

        if( print_xmin)
          {
            for( unsigned int i = 0; i < m_minimizer_hist[0].size(); i++ )
              {
                output.width(9);
                output << "x" << i << std::string(5,' ');
              }
          }

        output.width(9);
        output << "f" << std::string(5,' ');

        output.width(12);
        output << "norm" << std::endl;
        output << std::string(width,'-') << std::endl;
      }
  }

  void OptimizerMonitor::print_iteration( unsigned int iter, std::ostream& output,
                                          bool print_xmin ) const
  {
    if( m_env.fullRank() == 0 )
      {
        output.width(5);
        output << iter;

        if( print_xmin)
          {
            for( unsigned int i = 0; i < m_minimizer_hist[iter].size(); i++ )
              {
                output.width(2);
                output << "  ";
                output.width(13);
                output << std::scientific << m_minimizer_hist[iter][i];
              }
          }

        output.width(2);
        output << "  ";
        output.width(13);
        output << std::scientific << m_objective_hist[iter];

        output.width(2);
        output << "  ";
        output.width(13);
        output << m_norm_hist[iter] << std::endl;
      }
  }

  void OptimizerMonitor::reset()
  {
    m_display_conv = false;
    m_print_xmin = false;

    m_minimizer_hist.clear();
    m_objective_hist.clear();
    m_norm_hist.clear();
  }

  void OptimizerMonitor::print( std::ostream& output, bool print_xmin ) const
  {
    // First check that there's something to print.
    if( m_norm_hist.empty() )
      {
        std::cerr << "Nothing to print from OptimizerMonitor!" << std::endl;
        queso_error();
      }

    this->print_header(output,print_xmin);

    unsigned int size = m_norm_hist.size();
    for(unsigned int i = 0; i < size; i++ )
      {
        this->print_iteration(i,output,print_xmin);
      }
  }

} // end namespace QUESO
