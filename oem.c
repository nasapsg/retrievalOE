// ----------------------------------------------------------------------------
// PSG optimal estimation module (OEM) - Liuzzi and Villanueva
// Latest updates Villanueva, NASA/GSFC, December 2024
// ----------------------------------------------------------------------------
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "oem.h"

// Main program that performs the fit
void oem_fit(oem_func forward, struct oem_par *pars, struct oem_data *data, struct oem_result *results, struct oem_config config)
{
  long npar = results->npar, i, j;
  long npts = data->npts, l;
  double delta, cost0, lcost, cost=0.0, dcost=0.0, lambda; int err;

  // Allocate operational arrays and vectors
  double *fper, *pd;
  fper     = malloc(npts*sizeof(double));   // Perturbed model
  pd       = malloc(npar*sizeof(double));   // Temporary parameter for the derivative
  oem_indx = malloc(npar*sizeof(long));     // Index of LU decomposition (used by LU inverse)
  oem_p    = malloc(npar*sizeof(double));   // p vector used by matrix inversion
  oem_b    = malloc(npar*sizeof(double));   // b vector used by matrix inversion
  oem_r    = malloc(npar*sizeof(double));   // r vector used by matrix inversion
  oem_vv   = malloc(npar*sizeof(double));   // vv vector used by LU matrix inversion
  oem_a    = mdim(npar, npar, 0);           // Temporary copy of the matrix to be inverted
  oem_diff = mdim(npts, 1, 0);              // Vector containing the model difference
  oem_diffT= mdim(1, npts, 0);              // Transpose of the model difference
  oem_npt2 = mdim(npts, npts, 0);           // Temporary npts square matrix
  oem_npa2 = mdim(npar, npar, 0);           // Temporary npar square matrix
  oem_onem = mdim(1, 1, 0);                 // Temporary 1x1 matrix
  oem_pdff = mdim(npar, 1, 0);              // Vector containing the parameter difference
  oem_pdffT= mdim(1, npar, 0);              // Transpose of the parameter difference
  oem_kT   = mdim(npar, npts, 0);           // Tranpose of Jacobians matrix
  oem_kT_sei = mdim(npar, npts, 0);         // Operational matrix
  oem_kT_sei_k = mdim(npar, npar, 0);       // Operational matrix
  oem_si   = mdim(npar, npar, 0);           // Inverse of the covariance matrix
  oem_s    = mdim(npar, npar, 0);           // Covariance matrix
  oem_par  = mdim(npar, 1, 0);              // Estimated new parameter vector

  // Allocate resources in the results structure (if necessary)
  if (results->prior==NULL) results->prior= malloc(npar*sizeof(double));
  if (results->pfit==NULL)  results->pfit = malloc(npar*sizeof(double));
  if (results->efit==NULL)  results->efit = malloc(npar*sizeof(double));
  if (results->fit==NULL)   results->fit  = malloc(npts*sizeof(double));
  if (results->ak==NULL)    results->ak   = mdim(npar, npar, 0);
  if (results->k==NULL)     results->k    = mdim(npts, npar, 0);

  // Prepare the parameters and their co-variances (if not provided)
  for (i=0;i<npar;i++) {
    if (pars[i].value<pars[i].minval) pars[i].value=pars[i].minval;
    if (pars[i].value>pars[i].maxval) pars[i].value=pars[i].maxval;
    if (pars[i].variance==0.0) pars[i].variance = (pars[i].maxval - pars[i].minval)/(2.0*5.0);
    if (pars[i].maxval==0.0 && pars[i].minval==0.0) { pars[i].maxval = pars[i].value+pars[i].variance*5.0;  pars[i].minval = pars[i].value-pars[i].variance*5.0; }
    if (pars[i].step==0.0) pars[i].step=1e-3*pars[i].maxval;
    results->prior[i] = pars[i].value;
    results->pfit[i] = pars[i].value;
  }

  // Define inverse of parameter covariance matrix
  oem_sai = mdim(npar, npar, 1);
  if (data->sa==NULL) {
    for (i=0;i<npar;i++) oem_sai[i][i] = 1.0/(pars[i].variance*pars[i].variance);
  } else {
    err = minverse(oem_sai, data->sa, npar, 1);
    if (err!=0) { results->status=-2; return; }
  }
  for (i=0;i<npar;i++) for (j=0;j<npar;j++) oem_sai[i][j] *= config.gamma + TINY; // Scale covariances by Gamma

  // Define inverse of spectral covariance matrix
  oem_sei = mdim(npts, npts, 1);
  if (data->se==NULL) {
    for (l=0;l<npts;l++) oem_sei[l][l] = 1.0/(data->e[l]*data->e[l]);
  } else {
    err = minverse(oem_sei, data->se, npts, 1);
    if (err!=0) { results->status=-2; return; }
  }

  // Set defaults config parameters
  if (config.xtol<=0) config.xtol=1e-3;
  if (config.maxiter<=0) config.maxiter=20;
  results->niter=0; results->nfev=0; results->status=0; results->dof=0.0; results->chisq=0.0;

  // Obtain an apriori model
  forward(results->pfit, data->x, results->fit, npts, 1); results->nfev++;

  // Compute initial cost function
  cost0 = oem_cost(data, results, results->pfit, config.fpcost); lcost = cost0;

  // Iterate and search for the solution
  while (results->status==0 && npar>0) {
    // Compute derivatives and Jacobians
    for (i=0;i<npar;i++) {
      for (j=0;j<npar;j++) pd[j] = results->pfit[j];
      if ((pd[i]+pars[i].step) > pars[i].maxval) delta = -pars[i].step; else delta = pars[i].step;
      pd[i] += pars[i].step;
      forward(pd, data->x, fper, npts, 1); results->nfev++;
      for (l=0;l<npts;l++) results->k[l][i] = (fper[l] - results->fit[l])/delta + TINY*data->e[l];
    }

    // Perform a Levenberg-Marquard exploratory search of the solution near the current state
    lambda = 0.0; oem_gscale = 1.0; dcost = 1;
    while (oem_gscale<=100.0 && dcost>0) {
      err = oem_inv(data, results);
      if (err!=0) { results->status =-2; break; }         // Check for error in the matrix inversion process
      if (config.mode==1) { results->status = 1; break; } // Only computing Jacobians

      // Advance solution and explore braking the search
      while (lambda<=100.0 && dcost>0) {
        results->npegged=0;
        for (l=0;l<npar; l++) {
          pd[l] = results->pfit[l] + (oem_par[l][0] - results->pfit[l])/(1.0 + lambda);
          if (pd[l]<pars[l].minval) { pd[l] = pars[l].minval; results->npegged++; pars[l].status=1; }
          else if (pd[l]>pars[l].maxval) { pd[l] = pars[l].maxval; results->npegged++; pars[l].status=1; }
          else { pars[l].status=0; }
        }
        forward(pd, data->x, results->fit, npts, 1); results->nfev++;  // Re-compute solution with the new parameters
        cost = oem_cost(data, results, pd, config.fpcost);             // Compute cost function
        dcost = (cost-lcost)/lcost;                                    // Relative change of the cost function

        // Explore a Levenberg-Marquardt braking approach
        if ((fabs(dcost) < config.xtol) && lambda==0.0) break;
        else if (dcost>0 && lambda==0.0) lambda=1.0;
        else if (dcost>0 && lambda>0.0) lambda *= 10.0;
      }

      // Status of convergence and explore increasing the weight of the prior or even resetting to the current state
      if ((fabs(dcost) < config.xtol) && lambda==0.0) {  // Solution converged
        break; 
      } else if (dcost>0 && lambda>100.0) { // Try increasing the tug to the a-priori and try again
        lambda = 0.0;
        oem_gscale *= 10.0;
      }
    }

    // Check for convergence
    results->niter++;
    if (results->status!=0) {
      break;
    } else if ((fabs(dcost) < config.xtol) && lambda<0.1) {
      for (l=0;l<npar; l++) results->pfit[l] = pd[l];
      results->status = 1; // Convergence in chi-square
    } else if (results->niter>=config.maxiter) {
      for (l=0;l<npar; l++) results->pfit[l] = pd[l];
      results->status = 2; // Reached the maximum number of iterations     
    } else if (dcost>0 && cost<cost0) {
      for (l=0;l<npar; l++) results->pfit[l] = pd[l];
      results->status = 3; // Divergence of the siluation but reached a better state
    } else if (dcost>0) {
      results->status = -1; // The retrieval diverged and could not fit the data
    } else if (results->status==0) {
      for (l=0;l<npar; l++) results->pfit[l] = pd[l];
      lcost = cost;
    }
  }

  // Generate the latest model
  if (config.mode==0 && results->status>0 && npar>0) {
    forward(results->pfit, data->x, results->fit, npts, 1); results->nfev++;  // Generate final model
    cost = oem_cost(data, results, results->fit, config.fpcost);                            // Compute cost function
  }

  // Calculate resulting variables
  mmult(results->ak, oem_s, oem_kT_sei_k, npar, npar, npar);                 // Calculate averaging kernels S KT*Se^-1*K
  results->dof=0.0; for (j=0;j<npar;j++) results->dof += results->ak[j][j];  // DOF is sum of diagonal terms of AK
  for (i=0;i<npar;i++) results->efit[i] = sqrt(oem_s[i][i])/(results->ak[i][i] + TINY);

  // Release the memory
  free(fper); free(pd);
  free(oem_indx); free(oem_p); free(oem_b); free(oem_r); free(oem_vv); 
  mfree(oem_a, npar); mfree(oem_diff, npts); mfree(oem_diffT, 1); mfree(oem_npt2, npts);
  mfree(oem_npa2, npar); mfree(oem_onem, 1); mfree(oem_pdff, npar); mfree(oem_pdffT, 1);
  mfree(oem_kT, npar); mfree(oem_kT_sei, npar); mfree(oem_kT_sei_k, npar);
  mfree(oem_si, npar); mfree(oem_s, npar); mfree(oem_par, npar);
}

// Inversion function
int oem_inv(struct oem_data *data, struct oem_result *results)
{
  int err;
  long npts = data->npts, l;
  long npar = results->npar, i, j;

  // Transpose Jacobians (K) [npts x npar] -> (KT) [npar x npts]
  mtranspose(oem_kT, results->k, npts, npar);

  // Calculate KT*Se^-1 [npar x npts]
  mmult(oem_kT_sei, oem_kT, oem_sei, npar, npts, npts);

  // Calculate KT*Se^-1*K  [npar x npar]
  mmult(oem_kT_sei_k, oem_kT_sei, results->k, npar, npts, npar);

  // Calculate S_a^-1 + KT*Se^-1*K [npar x npar]
  for (i=0;i<npar;i++) for (j=0;j<npar;j++) oem_si[i][j] = oem_gscale*oem_sai[i][j] + oem_kT_sei_k[i][j];

  // Calculate the covariance matrix (S_a^-1 + KT*Se^-1*K)^-1 [npar x npar]
  err = minverse(oem_s, oem_si, npar, 1);
  if (err!=0) return err;

  // Calculate linearized model variation K(x-xa) [npts x 1]
  for (l=0;l<npar; l++) oem_pdff[l][0] = results->pfit[l] - results->prior[l];
  mmult(oem_diff, results->k, oem_pdff, npts, npar, 1);

  // Add model difference, y-f + K(x-xa) [npts x 1]
  for (l=0;l<npts; l++) oem_diff[l][0] += data->y[l] - results->fit[l];

  // Calculate KT*Se^-1(y-f + K(x-xa)) [npar x 1]
  mmult(oem_pdff, oem_kT_sei, oem_diff, npar, npts, 1);

  // Calculate xi+1 = xi + (S_a^-1 + KT*Se^-1*K)^-1 KT*Se^-1(y-f + K(x-xa)) [npar x 1]
  mmult(oem_par, oem_s, oem_pdff, npar, npar, 1);
  for (l=0;l<npar; l++) oem_par[l][0] += results->prior[l];
  return 0;
}


// Program to compute the cost function
double oem_cost(struct oem_data *data, struct oem_result *results, double *pd, double fpcost)
{
  double cost;
  long npts=data->npts, npar=results->npar, l;

  // Calculate chi-square component of the cost
  for (l=0;l<npts; l++) oem_diff[l][0] = data->y[l] - results->fit[l];
  mtranspose(oem_diffT, oem_diff, npts, 1);
  mmult(oem_npt2, oem_diffT, oem_sei, 1, npts, npts);
  mmult(oem_onem, oem_npt2, oem_diff, 1, npts, 1);
  results->chisq = oem_onem[0][0];
  cost = results->chisq;

  // Calculate the cost based on deviation of the parameter space from the a-priori
  if (fpcost>0 && npar>0) {
    for (l=0;l<npar; l++) oem_pdff[l][0] = pd[l] - results->prior[l];
    mtranspose(oem_pdffT, oem_pdff, npar, 1);
    mmult(oem_npa2, oem_pdffT, oem_sai, 1, npar, npar);
    mmult(oem_onem, oem_npa2, oem_pdff, 1, npar, 1);
    cost += fpcost*oem_onem[0][0];
  }
  results->cost = cost;
  return cost;
}

// ----------------------------------------------------------------------------
// Matrix operations
// ----------------------------------------------------------------------------
// Allocate array (and set to 0 if asked)
double **mdim(long nx, long ny, int clear)
{
  double **m = malloc(nx*sizeof(double*)); long l;
  for (l=0;l<nx;l++) {
    m[l] = malloc(ny*sizeof(double));
    if (clear) memset(m[l], 0, ny*sizeof(double));
  }
 return m;
}

// De-allocate memory of array
void mfree(double **m, long nx)
{
  long l;
  for (l=0;l<nx;l++) free(m[l]);
  free(m);
}

// Matrix transpose
void mtranspose(double **knt, double **kn, long nx, long ny)
{
  long x,y;
  for (x=0;x<nx;x++) for (y=0;y<ny;y++) knt[y][x] = kn[x][y];
}

// Matrix multiplication C[nx:nz] = A[nx:ny] x B[ny:nz]
void mmult(double **c, double **a, double **b, long nx, long ny, long nz)
{
  long x,y,z; double sum=0.0; int mode=0;

  // Search for diagonality (simple test) to improve performance
  if (nx==ny && ny>1) if (a[0][1]==0.0) mode=1;
  if (ny==nz && ny>1) if (b[0][1]==0.0) { if (mode==1) mode=3; else mode=2; }

  // Perform the multiplication
  if (mode==0) {           // Standard-mode
    for (x=0;x<nx;x++) {
      for (z=0;z<nz;z++) {
        for (y=0;y<ny;y++) sum+= a[x][y]*b[y][z];
        c[x][z] = sum;
        sum = 0.0;
      }
    }
  } else if (mode==1) {   // A is diagonal+quadratic
    for (x=0;x<nx;x++) for (z=0;z<nz;z++) c[x][z] = a[x][x]*b[x][z];
  } else if (mode==2) {   // B is diagonal+quadratic
    for (x=0;x<nx;x++) for (z=0;z<nz;z++) c[x][z] = a[x][z]*b[z][z];
  } else {                // A and B are diagonal and quadratic
    for (x=0;x<nx;x++) c[x][x] = a[x][x]*b[x][x];
  }
}

// ludcmp: Given a matrix a[1..n][1..n], this routine replaces it
// by the LU decomposition of a rowwise permutation of itself
// This routine is used in combination with lubksb to solve linear equations or invert a matrix.
int ludcmp(double **a, long n, long *indx, double *d)
{
  long i,imax=0,j,k;
  double big,dum,sum,temp;

  *d=1.0;
  for (i=0;i<n;i++) {
    big=0.0;
    for (j=0;j<n;j++) if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big==0.0) return -1; // Singular matrix in routine ludcmp
    oem_vv[i]=1.0/big;
  }
  for (j=0;j<n;j++) {
    for (i=0;i<j;i++) {
      for (k=0,sum=a[i][j];k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<n;i++) {      
      for (k=0,sum=a[i][j];k<j;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=oem_vv[i]*fabs(sum)) >= big) {
        big=dum;
        imax=i;
      }
    }
    if (j != imax) {
      for (k=0;k<n;k++) {
        dum=a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum;
      }
      *d = -(*d);
      oem_vv[imax]=oem_vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    if (j != n-1) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<n;i++) a[i][j] *= dum;
    }
  }
  return 0;
}

// lubksb: solves the set of n linear equations A·X = B.
// This routine takes into account the possibility that b will begin with many zero elements,
// so it is efficient for use in matrix inversion.
void lubksb(double **a, long n, long *indx, double b[])
{
  long i,ii=0,ip,j; int fsum=0;
  double sum;

  for (i=0;i<n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (fsum) for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) { ii=i; fsum=1; }
    b[i]=sum;
  }
  for (i=n-1;i>=0;i--) {
    sum=b[i];
    for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}


// This routine constructs the Cholesky decomposition, A = L · LT of a matrix
// From numerical recipes
// https://www.grad.hr/nastava/gs/prg/NumericalRecipesinC.pdf
int choldc(double **a, int n, double *p)
{
  int i,j,k;
  double sum;
  for (i=0;i<n;i++) {
    for (j=i;j<n;j++) {
      for (sum=a[i][j],k=i-1;k>=0;k--) sum -= a[i][k]*a[j][k];
      if (i==j) {
        if (sum <= 0.0) return -1; // The input matrix for choldc is not positive definite
        p[i]=sqrt(sum); 
      } else {
        a[j][i]=sum/p[i];
      }
    }
  }
  return 0;
}

// Solves the set of n linear equations A · x = b via the Cholesky decomposed matrix
// where a is a positive-definite symmetric matrix
// From numerical recipes
void cholsl(double **a, int n, double *p, double *b, double *x)
{
  int i,k;
  double sum;
  for (i=0;i<n;i++) { // Solve L · y = b, storing y in x.
    for (sum=b[i],k=i-1;k>=0;k--) sum -= a[i][k]*x[k];
    x[i]=sum/p[i];
  }
  for (i=n-1;i>=0;i--) { // Solve LT · x = y.
    for (sum=x[i],k=i+1;k<n;k++) sum -= a[k][i]*x[k];
    x[i]=sum/p[i];
  }
}

// Matrix inversion
int minverse(double **a_inv, double **a, long n, int mode)
{
  long i, j, err;
  if (n<=1) mode=2;
  else if (a[0][1]==0.0 && a[1][0]==0.0) mode=2; // Simple test for diagonal matrix

  // Inversion via LU decomposition
  if (mode==0) {
    for(i=0;i<n;i++) for(j=0;j<n;j++) oem_a[i][j]=a[i][j];
    err = ludcmp(oem_a,n,oem_indx,oem_p);
    if (err!=0) return -1;
    for(j=0;j<n;j++) {
      for(i=0;i<n;i++) oem_b[i]=0.0;
      oem_b[j]=1.0;
      lubksb(oem_a,n,oem_indx,oem_b);
      for(i=0;i<n;i++) a_inv[i][j]=oem_b[i];
    }

  // Inversion via Cholesky decomposition
  } else if (mode==1) {
    for(i=0;i<n;i++) for(j=0;j<n;j++) oem_a[i][j]=a[i][j];
    err = choldc(oem_a,n,oem_p);
    if (err!=0) return -2;
    for(j=0;j<n;j++) {
      for(i=0;i<n;i++) oem_b[i]=0.0;
      oem_b[j]=1.0;
      cholsl(oem_a,n,oem_p,oem_b,oem_r);
      for(i=0;i<n;i++) a_inv[i][j]=oem_r[i];
    }

  // Simple inversion of the diagonal terms
  } else {
    for(i=0;i<n;i++) {
      if (n>1) memset(a_inv[i], 0, n*sizeof(double));
      a_inv[i][i]=1.0/a[i][i];
    }
  }
  return 0;
}
