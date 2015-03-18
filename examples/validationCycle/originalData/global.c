// ******************* Global Arrhenius Kinetics ******************************
//                    (first order, one reaction)
//                          August 6, 2008
//
// Inputs:    (1) Temperature ramp rate:      beta (K/min)
//            (2) Pre-exponential:            A (s^-1)
// 	      (3) Activation energy:          E (kJ)
// 	      (4) Experimental data:          Me(i) vs. te(i) (K)
//
// Output:    (1) Sum of squared errors:      E2 = sum_i(Me(i)-Mt(i))


#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#define R 8.314472

int func(double t, const double M[], double f[], void *info)
{
  double* params = (double *)info;
  double A = params[0],E = params[1],beta = params[2];
  f[0] = -A*M[0]*exp(-E/(R*t))/beta;
  return GSL_SUCCESS;
}

#if 0
int jac(double t, const double M[], double *dfdM, double dfdt, void *info)
{
  double* params = (double *)info;
  double A = params[0],E = params[1],beta = params[2];
  dfdM = -A*exp(-E/(R*t))/beta;
  dfdt = -A*M[0]*E*exp(-E/(R*t))/(beta*R*t*t);
  return GSL_SUCCESS;
}
#endif

int main()
{

  double E2=0.;
  double A,E,beta,te[11],Me[11],Mt[11];
  int num_data;
  FILE *inp, *outp;

  inp = fopen("global.dat","r");
  fscanf(inp,"%lf %lf %lf",&A,&E,&beta);    /* read kinetic parameters */

  beta/=60.;	/* Convert heating rate to K/s */

  double params[]={A,E,beta};

  // read experimental data
  int i=0;
  int status;
  while (1){
    status = fscanf(inp,"%lf %lf",&te[i],&Me[i]);
    if (status == EOF) break;
    i++;
  }
  num_data = i;
  fclose(inp);

  // integration
  const gsl_odeiv_step_type * T = gsl_odeiv_step_gear1;

  gsl_odeiv_step *s = gsl_odeiv_step_alloc(T,1);
  gsl_odeiv_control *c = gsl_odeiv_control_y_new(1e-6,0.0);
  gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc(1);

  gsl_odeiv_system sys = {func, NULL, 1, (void *)params};

  double t = 0.1, t1 = 800.;
  double h = 1e-3;
  double M[1];
         M[0]=1.;

  outp = fopen("global.out","w");
  fprintf(outp,"Temp (K)    M\n------------------\n");

  i=0;
  double t_old=0., M_old[1];
  M_old[0]=1.;

  while (t < t1){
    int status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, t1, &h, M);

    if (status != GSL_SUCCESS) break;

    fprintf(outp,"%6.1lf %10.4lf\n",t,M[0]);

    if ( (t >= te[i]) && (t_old <= te[i]) ) {
      Mt[i] = (te[i]-t_old)*(M[0]-M_old[0])/(t-t_old) + M_old[0];
      E2+=(Me[i]-Mt[i])*(Me[i]-Mt[i]);
      // fprintf(outp,"%i %lf %lf %lf %lf\n",i,te[i],Me[i],Mt[i],E2);
      i++;
    }

    t_old=t;
    M_old[0]=M[0];

  }

  fprintf(outp,"For A = %g, E = %g, and beta = %.3lf\n",A,E,beta);
  fprintf(outp,"the sum of squared errors is %lf.",E2);

  gsl_odeiv_evolve_free(e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free(s);

  fclose(outp);

  return 0;

}
