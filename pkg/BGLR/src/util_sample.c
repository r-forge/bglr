#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rconfig.h>

#ifdef SUPPORT_OPENMP
  #include <omp.h>
  #define CSTACK_DEFNS 7
  #include "Rinterface.h"
#endif

/*
This function extracts a set of columns, the columns of the matrix are stacked 
in the vector *X (dim=nrows*ncols), columns=1,...,ncols
*/

SEXP extract_column(SEXP columns, SEXP rows, SEXP X)
{
   int i,j,r,jump;
   double *pX, *pxj;
   int *pcolumns;
   
   SEXP xj;
   
   r=INTEGER_VALUE(rows);
    
   PROTECT(X=AS_NUMERIC(X));
   pX=NUMERIC_POINTER(X);
   
   PROTECT(columns=AS_INTEGER(columns));
   pcolumns=INTEGER_POINTER(columns);
   
   PROTECT(xj = NEW_NUMERIC(r*LENGTH(columns)));
   
   pxj=NUMERIC_POINTER(xj);
   
   jump=-1;
   
   for(j=0; j<LENGTH(columns);j++)
   {
      for(i=0; i<r;i++) 
      {
	jump++;
	pxj[jump]=pX[i+(pcolumns[j]-1)*r];
      }
   } 
   
   UNPROTECT(3);
   
   return(xj);
}

/*
 * This is a generic function to sample betas in various models, including 
 * Bayesian LASSO, BayesA, Bayesian Ridge Regression, BayesCpi, etc.
 
 * For example, in the Bayesian LASSO, we wish to draw samples from the full 
 * conditional distribution of each of the elements in the vector bL. The full conditional 
 * distribution is normal with mean and variance equal to the solution (inverse of the coefficient of the left hand side)
 * of the following equation (See suplementary materials in de los Campos et al., 2009 for details),
   
    (1/varE x_j' x_j + 1/(varE tau_j^2)) bL_j = 1/varE x_j' e
 
    or equivalently, 
    
    mean= (1/varE x_j' e)/ (1/varE x_j' x_j + 1/(varE tau_j^2))
    variance= 1/ (1/varE x_j' x_j + 1/(varE tau_j^2))
    
    xj= the jth column of the incidence matrix
    
 *The notation in the routine is as follows:
 
 n: Number of rows in X
 pL: Number of columns in X
 XL: the matrix X stacked by columns
 XL2: vector with x_j' x_j, j=1,...,p
 bL: vector of regression coefficients
 e: vector with residuals, e=y-yHat, yHat= predicted values
 varBj: vector with variances, 
	For Bayesian LASSO, varBj=tau_j^2 * varE, j=1,...,p
	For Ridge regression, varBj=varB, j=1,...,p, varB is the variance common to all betas.
	For BayesA, varBj=varB_j, j=1,...,p
	For BayesCpi, varBj=varB, j=1,...,p, varB is the variance common to all betas
	
 varE: residual variance
 minAbsBeta: in some cases values of betas near to zero can lead to numerical problems in BL, 
             so, instead of using this tiny number we assingn them minAbsBeta
 
 */


SEXP sample_beta(SEXP n, SEXP pL, SEXP XL, SEXP xL2, SEXP bL, SEXP e, SEXP varBj, SEXP varE, SEXP minAbsBeta, SEXP ncores)
{
    double *xj, *pXL, *pxL2, *pbL, *pe, *pvarBj;
    double rhs,c,sigma2e, smallBeta;
    int j,i, rows, cols;
    int useCores, haveCores;
    SEXP list;
   
	
    #ifdef SUPPORT_OPENMP
	  R_CStackLimit=(uintptr_t)-1;
          useCores=INTEGER_VALUE(ncores);
          haveCores=omp_get_num_procs();
          if(useCores<=0 || useCores>haveCores) useCores=haveCores;
          omp_set_num_threads(useCores);
    #endif
	
    GetRNGstate();
	
    rows=INTEGER_VALUE(n);
    cols=INTEGER_VALUE(pL);
    sigma2e=NUMERIC_VALUE(varE);
    smallBeta=NUMERIC_VALUE(minAbsBeta);
	
    PROTECT(XL=AS_NUMERIC(XL));
    pXL=NUMERIC_POINTER(XL);

    PROTECT(xL2=AS_NUMERIC(xL2));
    pxL2=NUMERIC_POINTER(xL2);

    PROTECT(bL=AS_NUMERIC(bL));
    pbL=NUMERIC_POINTER(bL);

    PROTECT(e=AS_NUMERIC(e));
    pe=NUMERIC_POINTER(e);

    PROTECT(varBj=AS_NUMERIC(varBj));
    pvarBj=NUMERIC_POINTER(varBj);

    xj=(double *) R_alloc(rows,sizeof(double));

    for(j=0; j<cols;j++)
    {
	  rhs=0;
	  #pragma omp parallel 
	  {
	    #pragma omp for reduction(+:rhs) schedule(static)
	    for(i=0; i<rows; i++)
	    {
	      xj[i]=pXL[i+j*rows];
	      pe[i] = pe[i] + pbL[j]*xj[i];
	      rhs+=xj[i]*pe[i];
	    }
	  }
	  rhs=rhs/sigma2e;
  	  c=pxL2[j]/sigma2e + 1.0/pvarBj[j];
	  pbL[j]=rhs/c + sqrt(1.0/c)*norm_rand();
	  
	  #pragma omp for schedule(static)
	  for(i=0; i<rows; i++)
	  {
	    pe[i] = pe[i] - pbL[j]*xj[i];
	  }
          if(fabs(pbL[j])<smallBeta)
          {
             pbL[j]=smallBeta;
          }
        }
        
        // Creating a list with 2 vector elements:
	PROTECT(list = allocVector(VECSXP, 2));
	// attaching bL vector to list:
	SET_VECTOR_ELT(list, 0, bL);
	// attaching e vector to list:
	SET_VECTOR_ELT(list, 1, e);

  	PutRNGstate();

  	UNPROTECT(6);

  	return(list);
}

/*
 * This routine will be used to do updating of the coefficients in Bayes C-pi
 * OpenMP version
 */

SEXP d_e(SEXP p, SEXP n, SEXP X, SEXP d, SEXP b, SEXP error, SEXP varE, SEXP probInside,SEXP ncores)
{
  int i,j,rows,cols;
  double sigma2e, probIn, sum1,sum2, logOdds,tmp;
  //double logProbIn, logProbOut;
  double *pX, *perror, *pb, *eIn, *eOut;
  int *pd;
  int useCores, haveCores;
  SEXP list;
  
  cols=INTEGER_VALUE(p);
  rows=INTEGER_VALUE(n);
  sigma2e=NUMERIC_VALUE(varE);
  probIn=NUMERIC_VALUE(probInside);
    
  PROTECT(X=AS_NUMERIC(X));
  pX=NUMERIC_POINTER(X);

  PROTECT(d=AS_INTEGER(d));
  pd=INTEGER_POINTER(d);

  PROTECT(b=AS_NUMERIC(b));
  pb=NUMERIC_POINTER(b);
  
  PROTECT(error=AS_NUMERIC(error));
  perror=NUMERIC_POINTER(error);
  
  eOut=(double *) R_alloc(rows,sizeof(double));
  eIn=(double *) R_alloc(rows,sizeof(double));
 
  #ifdef SUPPORT_OPENMP
          R_CStackLimit=(uintptr_t)-1;
          useCores=INTEGER_VALUE(ncores);
          haveCores=omp_get_num_procs();
          if(useCores<=0 || useCores>haveCores) useCores=haveCores;
          omp_set_num_threads(useCores);
  #endif
 
  GetRNGstate();
  for(j=0; j<cols; j++)
  { 
     //Rprintf("%d\n",pd[j]);
     if(pd[j]==1)
     {
       #pragma omp for schedule(static)
       for(i=0; i<rows;i++)
       {
         eIn[i]=perror[i];
	 eOut[i]=perror[i] + pX[i+j*rows] * pb[j];
       }
     }else{
        #pragma omp for schedule(static)
        for(i=0; i<rows;i++)
	{
	  eOut[i]=perror[i];
	  eIn[i]=perror[i] - pX[i+j*rows] * pb[j];
	}
     }
     sum1=0;
     sum2=0;
     
     #pragma omp parallel 
     {
        #pragma omp for reduction(+:sum1,sum2) schedule(static)
	for(i=0; i<rows;i++)
	{
	  //sum1+=dnorm(eIn[i],0,sqrt(sigma2e),1); 
	  //sum2+=dnorm(eOut[i],0,sqrt(sigma2e),1);
            sum1+=eIn[i]*eIn[i];
            sum2+=eOut[i]*eOut[i];
	}
     }
     
     //logProbIn=log(probIn)+sum1;
     //logProbOut =log(1 - probIn)+sum2;
     //logOdds = logProbIn - logProbOut;

     logOdds=log(probIn/(1-probIn))+0.5/sigma2e*(sum2-sum1);
     tmp=exp(logOdds)/(1 + exp(logOdds));

     if(unif_rand()<tmp)
     {
         pd[j]=1;
     }else
     {
         pd[j]=0;
     }
     
     //Update error
     #pragma omp for schedule(static)
     for(i=0; i<rows;i++)
     {
         perror[i]=eOut[i] - pX[i+j*rows] * pb[j] * pd[j];
     }
     
  }
  
  // Creating a list with 2 vector elements
  PROTECT(list = allocVector(VECSXP, 2));

  // attaching b vector to list
  SET_VECTOR_ELT(list, 0,d);

  
  // attaching error vector to list:
  SET_VECTOR_ELT(list, 1, error);


  PutRNGstate();

  UNPROTECT(5);

  return(list);  
}


/*
Routines rgauss, rinvGauss and rinvGaussR taken from 
Package: SuppDists
Version: 1.1-8
Date: 2009/12/09
Title: Supplementary distributions
Author: Bob Wheeler <bwheelerg@gmail.com>
*/

void rgauss(double* normArray, int n, double mean, double sd)
{
        int i;
        GetRNGstate();
        for (i=0;i<n;i++) normArray[i]=rnorm(mean,sd);
        PutRNGstate();
}


/*
random inverse Gaussian values
Follows Mitchael,J.R., Schucany, W.R. and Haas, R.W. (1976). Generating
random roots from variates using transformations with multiple roots.
American Statistician. 30-2. 88-91.
*/

void rinvGauss(double* normArray,int n,double mu,double lambda)
{ 
        double b=0.5*mu/lambda;
        double a=mu*b;
        double c=4.0*mu*lambda;
        double d=mu*mu;

        rgauss(normArray,n,0,1);
        GetRNGstate();
        for (int i=0;i<n;i++) {
                if (mu<=0 || lambda<=0) {
                        normArray[i]=NA_REAL;
                }
                else {
                        double u=unif_rand();
                        double v=normArray[i]*normArray[i];     // Chi-square with 1 df
                        double x=mu+a*v-b*sqrt(c*v+d*v*v);      // Smallest root
                        normArray[i]=(u<(mu/(mu+x)))?x:d/x;     // Pick x with prob mu/(mu+x), else d/x;
                        if (normArray[i]<0.0){
                                v=x;
                        }
                }
        }
        PutRNGstate();
}

// Random function for R
void rinvGaussR(
        double *nup,
        double *lambdap,
        int *Np,
        int *Mp,        // length of nup and lambdap
        double *valuep
)
{
        int N=*Np;
        int M=*Mp;
        int D;
        int j;
        int k;
        int loc;
        int cloc;
        double *tArray;

        if (M==1)
                rinvGauss(valuep,N,*nup,*lambdap);
        else { // Allow for random values for each element of nu and lambda
                D=(N/M)+((N%M)?1:0);
                tArray=(double *)S_alloc((long)D,sizeof(double));
                loc=0;
                for (j=0;j<M;j++) {
                        rinvGauss(tArray,D,nup[j],lambdap[j]);
                        for (k=0;k<D;k++) {
                                cloc=loc+k*M;
                                if (cloc<N)
                                        valuep[cloc]=tArray[k];
                                else break;
                        }
                        loc++;
                }
        }
}

