// ----------------------------------------------------------------------------
// Optimal estimation module (OEM) - Header
// Liuzzi and Villanueva - NASA/GSFC - Sep/2018
// ----------------------------------------------------------------------------
#define TINY 1.0e-20

// Parameter structure 
struct oem_par {
  double value;        // A-priori value / starting value
  double back;         // Estimation background value  (if maxval/minval defined, then background is replaced with max+min/2)
  double variance;     // Variance of the parameter (if maxval/minval defined, then variance is replaced with max-min/2)
  double maxval;       // maximum limit boundary value
  double minval;       // lower limit boundary value
  double step;         // Step size for finite difference 
  double relstep;      // Relative step size for finite difference
  int status;          // Status of the fit for this variable (0:Normal, 1:Touched the limits, 2: Unconstrained)
  int side;            // Sidedness of finite difference derivative 
                       //     0 - one-sided forward derivative (f(x+h) - f(x)  )/h
                       //     1 - one-sided reverse derivative (f(x)   - f(x-h))/h
                       //     2 - two-sided derivative (f(x+h) - f(x-h))/(2*h) 
};

// Data structure
struct oem_data {
  long npts;           // Number of points
  double *x;           // Ordinate vector
  double *y;           // Signal vector
  double *e;           // Error vector
  double **sy_inv;     // Covariance vector for signal (optional, diagonal calculated based on the error) 
  double **sa_inv;     // Covariance vector for the parameters (optional, diagonal calculated based on parameter variance) 
};

// Configuration structure 
struct oem_config {
  double gamma;        // Levenberg-Marquardt estimation parameter
  double xtol;         // Relative parameter convergence criterium  Default: 1e-4
  int maxiter;         // Maximum number of iterations.             Default: 50
  int mode;            // Runtime mode (0:Fit, 1:Calculate Jacobians/Kernels/functions)  
};

// Result structure
struct oem_result {
  double chisq;        // Final chi^2 
  int npar;            // Number of parameters
  int npegged;         // Number of variables locked to the limit
  int niter;           // Number of iterations 
  int nfev;            // Number of function evaluations 
  int status;          // Fitting status code
  double dof;          // Degrees of freedom
  double *pfit;        // Optimum parameters (1D npar vector)
  double *efit;        // Parameter uncertainties (1-sigma), 1D npar vector 
  double *fit;         // Optimum model
  double **kn;         // Jacobians (npts x npar)
  double **sn;         // Estimated error covariance matrix (npar x npar)
  double **dn;         // Contribution functions (npar x npts)
  double **ak;         // Averaging kernels (npar x npar)
};

// Forward fitting function format 
typedef void (*oem_func)( 
  double *p,           // Parameters for the forward model
  double *x,           // Ordinate vector
  double *y,           // Resulting function
  long npts);          // Number of points

// Prototypes
void oem_fit(oem_func forward, struct oem_par *pars, struct oem_data *data, struct oem_result *results, struct oem_config config);
int  oem_inv(struct oem_par *pars, struct oem_data *data, struct oem_result *results, struct oem_config config);
void oem_chisq(struct oem_data *data, struct oem_result *results);
void oem_err(struct oem_data *data, struct oem_result *results);
double **mdim(long nx, long ny, int clear);
void mfree(double **m, long nx, long ny);
void mtranspose(double **kn, double **knt, long nx, long ny);
void mmult(double **a, double **b, double **c, long nx, long ny, long nz);
void madd(double **a, double **b, double **c, long nx, long ny, int flag);
void lubksb(double **a, long n, long *indx, double b[]);
void ludcmp(double **a, long n, long *indx, double *d);
void minverse(double **a_inv, double **a, long n); 
