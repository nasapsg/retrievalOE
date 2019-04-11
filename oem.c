// ----------------------------------------------------------------------------
// Optimal estimation module (OEM) - Routines
// Liuzzi and Villanueva - NASA/GSFC - Sep/2018
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
  double delta, lchisq; int status;

  // Allocate resources in the results structure (if necessary)
  if (results->pfit==NULL) results->pfit = malloc(npar*sizeof(double));
  if (results->efit==NULL) results->efit = malloc(npar*sizeof(double));
  if (results->fit==NULL)  results->fit  = malloc(npts*sizeof(double));
  if (results->kn==NULL)   results->kn   = mdim(npts, npar, 0);
  if (results->sn==NULL)   results->sn   = mdim(npar, npar, 0);
  if (results->dn==NULL)   results->dn   = mdim(npar, npts, 0);
  if (results->ak==NULL)   results->ak   = mdim(npar, npar, 0);

  // Prepare the parameters and their co-variances (if not provided)
  for (i=0;i<npar;i++) {
    if (pars[i].back==0.0)     pars[i].back = (pars[i].maxval + pars[i].minval)/2.0;
    if (pars[i].variance==0.0) pars[i].variance = (pars[i].maxval - pars[i].minval)/(2.0*5.0);
    if (pars[i].maxval==0.0 && pars[i].minval==0.0) { pars[i].maxval = pars[i].back+pars[i].variance*5.0;  pars[i].minval = pars[i].back-pars[i].variance*5.0; }
    if (pars[i].value<pars[i].minval) pars[i].value=pars[i].minval; 
    if (pars[i].value>pars[i].maxval) pars[i].value=pars[i].maxval; 
    results->pfit[i] = pars[i].value;
  }
  if (data->sa_inv==NULL) {
    data->sa_inv = mdim(npar, npar, 1);
    for (i=0;i<npar;i++) data->sa_inv[i][i] = 1.0/(pars[i].variance*pars[i].variance);
  }  
  
  // Calculate the data co-variance (if not provided)
  if (data->sy_inv==NULL) {
    data->sy_inv = mdim(npts, npts, 1);
    for (l=0;l<npts;l++) data->sy_inv[l][l] = 1.0/(data->e[l]*data->e[l]);
  }
  
  // Set defaults config parameters
  if (config.xtol<=0) config.xtol=1e-4;
  if (config.maxiter<=0) config.maxiter=50;
  results->niter=0; results->nfev=0; results->status=0; results->dof=0.0; results->chisq=0.0;
  
  // Prepare the variables for the Jacobians    
  double *ylo = malloc(npts*sizeof(double));   // Lower bound of the derivate
  double *yhi = malloc(npts*sizeof(double));   // Upper bound of the derivate
  double *pd  = malloc(npar*sizeof(double));   // Temporary parameter for the derivative
  
  // Obtain an apriori model
  forward(results->pfit, data->x, results->fit, npts); results->nfev++;
  
  // Compute chi-square
  oem_chisq(data, results);
  lchisq = results->chisq;
  if (npar==0) return;
  
  // Iterate and search for the solution  
  while (results->niter<config.maxiter && results->status==0) {
    // Compute derivatives and Jacobians
    results->npegged=0;
    for (i=0;i<npar;i++) {
      if (pars[i].status>1) { results->npegged++; continue; }
      for (j=0;j<npar;j++) pd[j] = results->pfit[j];
      if (pars[i].relstep>0) {
        delta = pars[i].back*pars[i].relstep;
      } else if (pars[i].step>0) {
        delta = pars[i].step;
      } else {
        delta = TINY;
      }
      if (pars[i].side==0) {
        pd[i] += delta;
        forward(pd, data->x, yhi, npts); results->nfev++;
        for (l=0;l<npts;l++) ylo[l]=results->fit[l];
      } else if (pars[i].side==1) {
        pd[i] -= delta;
        forward(pd, data->x, ylo, npts); results->nfev++;
        for (l=0;l<npts;l++) yhi[l]=results->fit[l];
      } else {
        pd[i] = results->pfit[i] - delta;
        forward(pd, data->x, ylo, npts); results->nfev++;
        pd[i] = results->pfit[i] + delta;
        forward(pd, data->x, yhi, npts); results->nfev++;
        delta *= 2.0;
      }
      for (l=0;l<npts;l++) results->kn[l][i] = (yhi[l] - ylo[l])/delta + TINY*data->e[l];
    }
                              
    // Perform the optimal estimation
    status = oem_inv(pars, data, results, config);
    results->dof=0.0; for (j=0;j<npar;j++) results->dof += results->ak[j][j];
    if (config.mode==1) return;
      
    // Re-compute solution with the new parameters
    forward(results->pfit, data->x, results->fit, npts); results->nfev++;
      
    // Compute chi-square
    oem_chisq(data, results);
      
    // Convergence criteria
    if (results->npegged >= npar) {
      results->status =-2; // All variables are pegged
    } else if (status>0) {
      results->status = 0; // Not there yet, some parameters need to be re-estimated
    } else if ((fabs(results->chisq - lchisq)/npts) < config.xtol) {
      results->status = 1; // Success - convergence in chi-square
    } else if (results->chisq > 2.0*lchisq) {
      results->status =-1; // Divergence in chi-square
    }
    lchisq = results->chisq;
    results->niter++;
  }
  if (results->niter>=config.maxiter) results->status=2;
  
  // Compute a-posteriori covariances of the parameters
  oem_err(data, results);
  
  // De-allocate temporary Jacobian variables
  free(yhi); free(ylo); free(pd);        
}


// Inversion function
int oem_inv(struct oem_par *pars, struct oem_data *data, struct oem_result *results, struct oem_config config)
{
  int status=0;
  long npts = data->npts, l;
  long npar = results->npar, i, j;

  // Transpose Jacobians (Kn) -> (KnT)
  double **knt = mdim(npar, npts, 0);
  mtranspose(knt, results->kn, npts, npar);

  // Calculate KnT*Sy^-1:
  double **knt_sy_inv = mdim(npar, npts, 1);
  mmult(knt_sy_inv, knt, data->sy_inv, npar, npts, npts);
  mfree(knt, npar, npts);

  // Calculate KnT*Sy^-1*Kn
  double **tmp = mdim(npar, npar, 1);
  mmult(tmp, knt_sy_inv, results->kn, npar, npts, npar);
  
  // Calculate KnT*Sy^-1*Kn + gamma*S_a^-1
  for (i=0;i<npar;i++) for (j=0;j<npar;j++) tmp[i][j] += (config.gamma+TINY)*data->sa_inv[i][j];

  // Calculate the inverse of the matrix  
  minverse(results->sn, tmp, npar);
  mfree(tmp, npar, npar);
  
  // Calculate the contribution functions
  mmult(results->dn, results->sn, knt_sy_inv, npar, npar, npts);
  mfree(knt_sy_inv, npar, npts);

  // Calculate the averaging kernels
  mmult(results->ak, results->dn, results->kn, npar, npts, npar);

  // Calculate the solution
  double *pnew = malloc(npar*sizeof(double));
  for (i=0;i<npar;i++) {
    pnew[i] = pars[i].back;
    for (l=0;l<npts;l++) pnew[i] += results->dn[i][l]*(data->y[l]-results->fit[l]);
    for (j=0;j<npar;j++) pnew[i] += results->ak[i][j]*(results->pfit[j]-pars[j].back);
    for (j=0;j<npar;j++) pnew[i] += results->sn[i][j]*config.gamma*data->sa_inv[j][j]*(results->pfit[j]-pars[j].back);
  }
  
  // Assign the value and check for limits
  for (i=0;i<npar;i++) {
    if (pars[i].status>1) continue;
    if (pnew[i]>pars[i].maxval || pnew[i]<pars[i].minval) {
      if (pnew[i]>pars[i].maxval) pnew[i] = pars[i].maxval;
      if (pnew[i]<pars[i].minval) pnew[i] = pars[i].minval;
      if (pars[i].status==0) {
        pars[i].status=1;
        status=1;
      } else {
        pars[i].status=2;
        status=1;
      }
    }
    results->pfit[i]=pnew[i];
  }
  free(pnew);
  return status;
}


// Program to compute the signal chi-square
void oem_chisq(struct oem_data *data, struct oem_result *results)
{
  long npts = data->npts, l;
  double **diff, **diffT, **tmp, **tmp2;
  
  // Calculate (y-fx)T * Sy_inv
  diff  = mdim(npts, 1, 0); for (l=0;l<npts; l++) diff[l][0] = data->y[l] - results->fit[l];
  diffT = mdim(1, npts, 0); mtranspose(diffT, diff, npts, 1);
  tmp   = mdim(1, npts, 0); mmult(tmp, diffT, data->sy_inv, 1, npts, npts);
  
  // Calculation of ((y-fx)T * Sy_inv) * (y-fx)
  tmp2  = mdim(1,    1, 0); mmult(tmp2, tmp, diff, 1, npts, 1);
  results->chisq = tmp2[0][0];
  mfree(diff, npts, 1); mfree(diffT, 1, npts); mfree(tmp, 1, npts); mfree(tmp2, 1, 1);
}


// Program to calculate the a-posteriori covariances of the parameters and the measurements
void oem_err(struct oem_data *data, struct oem_result *results)
{
  long npar = results->npar, i;
  long npts = data->npts;
  double **sy, **dn_sy, **dnt, **sxm;

  // Calculation of the parameters covariance matrix
  // Calculate Dn*Sy
  sy    = mdim(npts, npts, 1); minverse(sy, data->sy_inv, npts);
  dn_sy = mdim(npar, npts, 1); mmult(dn_sy, results->dn, sy, npar, npts, npts);
  mfree(sy, npts, npts);
  
  // Calculate Sxn_est = (Dn * Sy) * DnT
  dnt   = mdim(npts, npar, 0); mtranspose(dnt, results->dn, npar, npts);
  sxm   = mdim(npar, npar, 1); mmult(sxm, dn_sy, dnt, npar, npts, npar);
  for (i=0;i<npar;i++) results->efit[i] = sqrt(sxm[i][i])/results->ak[i][i];
  mfree(dn_sy, npar, npts); mfree(dnt, npts, npar); mfree(sxm, npar, npar);
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
void mfree(double **m, long nx, long ny) 
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

     
// ludcmp: Given a matrix a[1..n][1..n], this routine replaces it 
// by the LU decomposition of a rowwise permutation of itself
// This routine is used in combination with lubksb to solve linear equations or invert a matrix.
// Numerical recipes: https://www.astro.umd.edu/~ricotti/NEWWEB/teaching/ASTR415/NR3/legacy/nr2/C_211/index.htm
void ludcmp(double **a, long n, long *indx, double *d)
{
  long i,imax=0,j,k;
  double big,dum,sum,temp;
  double *vv = malloc(n*sizeof(double));
  
  *d=1.0;
  for (i=0;i<n;i++) {
    big=0.0;
    for (j=0;j<n;j++) if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big==0.0) { printf("Diagonal matrix\n"); return; }
    vv[i]=1.0/big;
  }
  for (j=0;j<n;j++) {
    for (i=0;i<j;i++) {
      sum=a[i][j];
      for (k=0;k<i;k++) { 
        sum -= a[i][k]*a[k][j];
      }
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<n;i++) {
      sum=a[i][j];
      for (k=0;k<j;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
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
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    if (j != n-1) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<n;i++) a[i][j] *= dum;
    }
  }
  free(vv);
}


// Matrix inversion
void minverse(double **a_inv, double **a, long n) 
{
  // Check for diagonality
  long mode=0, i, j; 
  if (n>1) { if (a[0][1]==0.0) mode=1; else mode=0; } else { mode=1; }
  
  // Perform the inversion
  if (mode==0) { 
    long   *indx = malloc(n*sizeof(long)); 
    memset(indx, 0, n*sizeof(long));
    double *col  = malloc(n*sizeof(double));
    double *d    = malloc(n*sizeof(double));
    ludcmp(a,n,indx,d); 
    for(j=0;j<n;j++) {
      for(i=0;i<n;i++) col[i]=0.0; 
      col[j]=1.0;
      lubksb(a,n,indx,col); 
      for(i=0;i<n;i++) a_inv[i][j]=col[i];
    }
    free(d); free(col); free(indx);
  } else {
    for(i=0;i<n;i++) a_inv[i][i]=1.0/a[i][i];  
  }
}

