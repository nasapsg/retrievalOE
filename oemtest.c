// ------------------------------------------------------------------------
// Program to test the optimal estimation module (oemfit)
// The main program implements also the loop to compute
// the optimal estimation solution.
// Liuzzi and Villanueva - NASA/GSFC - Sep/2018
// To compile: gcc oemtest.c oem.c -o oemtest
// ------------------------------------------------------------------------
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "oem.h" 

// Forward test model
void forward(double *p, double *x, double *y, long npts, int fitbase) 
{
  long i; double n0=pow(10.0,p[0]), ni;
  for (i=0;i<npts;i++) {
    ni = n0*exp(-p[1]*x[i]);
    ni = 1e18*ni*pow(x[i], 6.0);
    y[i] = log10(ni) + 1e-1*i*p[2];
  }  	  
}

// Main function
int main(int argc, char *argv[])
{
  // Simulate some test data  
  long npts=50, npar=3, i;
  double par0[]={2.3, 3000.0, 1.0}, ui;
  struct oem_data data; 
  data.npts = npts;
  data.x = malloc(npts*sizeof(double)); // Ordinate vector
  data.y = malloc(npts*sizeof(double)); // Flux vector
  data.e = malloc(npts*sizeof(double)); // Error vector
  data.se=NULL; data.sa=NULL;   // Force to compute covariances
  for (i=0;i<npts;i++) {
    ui = -4.0 + (2.0/5e1)*i;
    data.x[i] = pow(10.0, ui);
  }
  forward(par0, data.x, data.y, npts, 0);
  for (i=0;i<npts;i++) { 
    data.e[i]  = 0.1;
    data.y[i] += 2.0*((1.0*rand()/RAND_MAX)-0.5)*data.e[i];
  }

  // Prepare the parameters for the fit
  struct oem_par pars[npar];
  for (i=0; i<npar; i++) memset(&pars[i], 0, sizeof(struct oem_par));
  pars[0].value=4.0;     pars[0].minval=1.0;     pars[0].maxval=5.0;       pars[0].step=0.1*pars[0].value;
  pars[1].value=200.0;   pars[1].minval=100.0;   pars[1].maxval=6000.0;    pars[1].step=0.1*pars[1].value;
  pars[2].value=3.0;     pars[2].minval=0.1;     pars[2].maxval=10.0;      pars[2].step=0.1*pars[2].value;

  // Configuration parameters
  struct oem_config config;
  config.gamma = 1.0;  // Standard OE is 1
  config.xtol = 1e-3;  // Convergence tolerance
  config.maxiter = 50; // Restrict to 50 iterations
  config.fpcost = 0;   // The cost function only includes the signal chi-square

  // Reset results array
  struct oem_result results;  
  memset(&results, 0, sizeof(struct oem_result));
  results.npar = npar;

  // Perform the fit
  oem_fit(forward, pars, &data, &results, config);
  
  // Print the results
  printf("Chisq: %f (status %d)\n", results.chisq, results.status);
  for (i=0;i<npar;i++) printf("%e +/- %e\n", results.pfit[i], results.efit[i]);  
  return 1;  
}
