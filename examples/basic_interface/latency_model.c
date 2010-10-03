/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008,2009 The PECOS Development Team
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
 * Simple toy problem to illustrate the use of the Basic QUESO
 * interface to perform a statistical inverse problem.
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include <stdlib.h>
#include <stdio.h>
#include <queso_basic.h>
#include <mpi.h> 

char *queso_inputfile = "queso.inp";

double latency_likelihood(double *params);
void read_data(char *ifile,double *data);

static int num_data_pts = 3936;

int main(int argc, char *argv[])
{

  MPI_Init(&argc,&argv);

  printf("--> Initializing QUESO Environment...\n");

  /* Initialize QUESO environment and define input file */

  QUESO_init(queso_inputfile);	                   

 /* Register an application likelihood function with QUESO and run a
  * statistical inversion problem using MCMC sampling */

  QUESO_statistical_inversion(latency_likelihood);     

  /* Finalize the analysis and output statistics */

  QUESO_finalize();
  
  MPI_Finalize();
  return 0;

}

double latency_likelihood(double *params)
{
  
  int i;
  double likelihood;
  static double *data;
  static int     first_flag = 1;
  static double  variance = 1.*1.;  /* sigma^2 */

  if(first_flag)
    {
      data = (double *)calloc(num_data_pts,sizeof(double));

      if(data == NULL)
	{
	  printf("** Error: unable to alloc space for data array\n");
	  exit(1);
	}

      read_data("latency-all-nodes.txt",data);  
      first_flag = 0;
      printf("num_data_pts = %i\n",num_data_pts);
      printf("** End of first likelihood iteration\n");
      
    }

  likelihood = 0.0;

  for(i=0;i<num_data_pts;i++)
    likelihood += (params[0] - data[i])*(params[0] - data[i]);

  return( (likelihood)/(variance) );
}

void read_data(char *ifile,double *data)
{
  FILE *fp;
  //  char *ifile = "latency-all-nodes.txt";
  char tmp_string[20];
  float tmp_float;
  //  double *data;
  int i;



  fp = fopen(ifile,"r");

  if(fp == NULL)
    {
      printf("** Error: unable to open file %s",ifile);
      exit(1);
    }

   for(i=0;i<num_data_pts;i++)
    {
      fscanf(fp,"%8s %f",tmp_string,&tmp_float);
      data[i] = (double)tmp_float;
      printf("Read in %f for %s\n",data[i],tmp_string);
    }

  fclose(fp);

  printf("about to return from read_data\n");
  return;

}

