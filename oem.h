// ----------------------------------------------------------------------------
// PSG Optimal estimation module (OEM) - Header
// ----------------------------------------------------------------------------
#define TINY 1.0e-20

// Parameter structure 
struct oem_par {
  double value;        // A-priori value
  double minval;       // lower limit boundary value
  double maxval;       // maximum limit boundary value
  double scale;        // Normalization factor applied to the values
  double zero;         // Correction of the zero value to ensure always positive values
  double variance;     // Variance of the parameter (if maxval/minval defined, then variance is replaced with max-min/20)
  double step;         // Step size for finite difference 
  int status;          // Status of the fit for this variable (0:Normal, 1:At one of the limits)
};

// Data structure
struct oem_data {
  long npts;           // Number of points
  double *x;           // Ordinate vector
  double *y;           // Signal vector
  double *e;           // Error vector
  double **se;         // Covariance vector for signal (optional, diagonal calculated based on the error) 
  double **sa;         // Covariance vector for the parameters (optional, diagonal calculated based on parameter variance) 
};

// Configuration structure 
struct oem_config {
  double fpcost;       // Factor to weight the parameter deviation in the cost function (0 is simply chi-square)
  double gamma;        // Regularization parameter (0:Classic LM/no-apriori, 1:classic Rodgers' formalism, 10:Heavily tailored to the a-priori)
  double xtol;         // Relative parameter convergence criterium
  int maxiter;         // Maximum number of iterations
  int mode;            // Runtime mode (0:Fit, 1:Calculate Jacobians/Kernels/functions)  
};

// Result structure
struct oem_result {
  double cost;         // Final cost function
  double chisq;        // Final chi^2 
  int npar;            // Number of parameters
  int npegged;         // Number of variables locked to the limit
  int niter;           // Number of iterations 
  int nfev;            // Number of function evaluations 
  int status;          // Fitting status code
  double dof;          // Degrees of freedom 
  double *prior;       // A-priori value of the parameters
  double *pfit;        // Optimum parameters (1D npar vector)
  double *efit;        // Parameter uncertainties (1-sigma), 1D npar vector (npar)
  double *fit;         // Optimum model (npts)
  double **ak;         // Averaging kernels (npar x npar)
  double **k;          // Jacobians (npts x npar)
};

// Forward fitting function format 
typedef void (*oem_func)( 
  double *p,           // Parameters for the forward model
  double *x,           // Ordinate vector
  double *y,           // Resulting function
  long npts,           // Number of points
  int fitbase);        // Flag indicating if base/gain should be fitted

// Prototypes
void oem_fit(oem_func forward, struct oem_par *pars, struct oem_data *data, struct oem_result *results, struct oem_config config);
int oem_inv(struct oem_data *data, struct oem_result *results);
double oem_cost(struct oem_data *data, struct oem_result *results, double *pd, double fprior);
double **mdim(long nx, long ny, int clear);
void mfree(double **m, long nx);
void mtranspose(double **kn, double **knt, long nx, long ny);
void mmult(double **a, double **b, double **c, long nx, long ny, long nz);
void madd(double **a, double **b, double **c, long nx, long ny, int flag);
int minverse(double **a_inv, double **a, long n, int mode); 
