//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011 The PECOS Development Team
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
// 
// $Id$
//
//--------------------------------------------------------------------------

using namespace std;

#define _QUESO_Basic_API_DEF

#include <basic_interface.h>
#include <basic_classes.h>
#include <basic_int.h>
#include <grvy.h>

using namespace QUESO_Basic_API;

//-------------
// C Interface 
//-------------

extern "C" void QUESO_init(const char *inputfile)
{
  _QUESO_Basic = new QUESO_Basic_Class();
  _QUESO_Basic->Initialize(inputfile);
  return;
}

extern "C" void QUESO_statistical_inversion(double (*fp)(double *) )
{
  _QUESO_Basic->DefineParameterSpace ();
  _QUESO_Basic->Likelihood_Register  (fp);
  _QUESO_Basic->SolveInverseProblem  ();
  return;
}

extern "C" void QUESO_finalize()
{
  delete _QUESO_Basic;
  grvy_timer_finalize();

  printf("\n QUESO: Complete\n");

  grvy_timer_summarize();
  return;
}

//-------------------
// Fortran Interface 
//-------------------

char *f2c_char(char *,int);

extern "C" void queso_init_(char *inputfile,int _namelen)
{
  char *name = f2c_char(inputfile,_namelen);
  QUESO_init(name);

  delete[] name;
  return;
}

extern "C" void queso_statistical_inversion_(double (*fp)(double *) )
{
  QUESO_statistical_inversion(fp);
  return;
}

extern "C" void queso_finalize_()
{
  QUESO_finalize();
  return;
}


// f2c_char(): Convert evil Fortran character strings to C                 

char *f2c_char(char*input,int len)
{
  char* name = new char[len+1];

  strncpy(name,input,len);
  name[len]='\0';
  return(name);
}

