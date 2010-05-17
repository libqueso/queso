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
#include <queso.h>
#include <mpi.h> 
#include "grvy.h"

char *queso_inputfile = "queso.inp";

double latency_likelihood(double *params);
void read_data(char *ifile,double *data, int *hops);

//static int num_data_pts = 3936;
static int num_data_pts = 3892;
static int model = 2;

int main(int argc, char *argv[])
{

  MPI_Init(&argc,&argv);

  if(argc < 2)
    {
      printf("\nUsage: %s <input-file>\n\n",argv[0]);
      exit(1);
    }
  
  grvy_log_setlevel(GRVY_INFO);

  printf("--> Initializing QUESO Environment...\n");

  /* Initialize QUESO environment and define input file */

  QUESO_init(argv[1]);	                   

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
  static int num_params;
  double likelihood;
  double model_latency;
  static double *data;
  static    int *hops;
  static int     first_flag = 1;
  //  static double  variance = 1.*1.;  /* sigma^2 */
  //  static double  variance = 0.5;  /* sigma^2 */
  static double  variance = 0.01;  /* sigma^2 */

  if(first_flag)
    {
      data = (double *)calloc(num_data_pts,sizeof(double));
      hops =    (int *)calloc(num_data_pts,sizeof(int));

      if(data == NULL)
	{
	  printf("** Error: unable to alloc space for data array\n");
	  exit(1);
	}

      read_data("results.original",data,hops);  
      first_flag = 0;
      printf("num_data_pts = %i\n",num_data_pts);
      printf("** End of first likelihood iteration\n");
      
    }

  likelihood = 0.0;

  if(model == 1)
    {
      for(i=0;i<num_data_pts;i++)
	likelihood += (params[0] - data[i])*(params[0] - data[i]);
    }
  else if(model == 2)
    {
      for(i=0;i<num_data_pts;i++)
	{
	  model_latency = params[0] + params[1]*hops[i];
	  likelihood += (model_latency - data[i])*(model_latency - data[i]);
	  //	  printf("data = %f, model = %f, hops = %i\n",data[i],model_latency,hops[i]);
	}

      //      printf("alpha = %f, beta = %f\n",params[0],params[1]);
    }
  else
    {
      printf("unknown model\n");
      exit(1);
    }
	
      
  likelihood = likelihood/(1.*num_data_pts);

  return( (-likelihood)/(variance) );
}

void read_data(char *ifile,double *data, int *hops)
{
  FILE *fp;
  char tmp_string[20];
  float tmp_float;
  int   tmp_int;
  int i;

  fp = fopen(ifile,"r");

  if(fp == NULL)
    {
      printf("** Error: unable to open file %s",ifile);
      exit(1);
    }

   for(i=0;i<num_data_pts;i++)
    {
      fscanf(fp,"%8s %8s %f %i",tmp_string,tmp_string,&tmp_float,&tmp_int);
      data[i] = (double)tmp_float;
      hops[i] = (int)   tmp_int;
      printf("Read in %f (%i hops) for %s\n",data[i],hops[i],tmp_string);
    }

  fclose(fp);

  return;

}

