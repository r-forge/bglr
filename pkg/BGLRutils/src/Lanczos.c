/*  BGRLutils/src/Lanczos.c
 *
 *  mgcv_trisymeig and extreme_eigenvalues adapted from the routines
 *  ggcv_trisymeig and Rlanczos routines in the mgcv R package, Version 1.7-22, http://cran.r-project.org/web/packages/mgcv/index.html
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 or 3 of the License
 *  (at your option).
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  A copy of the GNU General Public License is available at
 *  http://www.r-project.org/Licenses/
 */


#include "Lanczos.h"

void mgcv_trisymeig(double *d,double *g,double *v,int *n,int getvec,int descending) 
/* Find eigen-values and vectors of n by n symmetric tridiagonal matrix 
   with leading diagonal d and sub/super diagonals g. 
   eigenvalues returned in d, and eigenvectors in columns of v, if
   getvec!=0. If *descending!=0 then eigenvalues returned in descending order,
   otherwise ascending. eigen-vector order corresponds.  

   Routine is divide and conquer followed by inverse iteration. 

   dstevd could be used instead, with just a name change.
   dstevx may be faster, but needs argument changes.
*/ 
		    
{ char compz;
  double *work,work1,x,*dum1,*dum2;
  int ldz=0,info,lwork=-1,liwork=-1,*iwork,iwork1,i,j;

  if (getvec) { compz='I';ldz = *n;} else { compz='N';ldz=0;}

  /* workspace query first .... */
  F77_NAME(dstedc)(&compz,n,
		   d, g, /* lead and su-diag */
		   v, /* eigenvectors on exit */  
                   &ldz, /* dimension of v */
		   &work1, &lwork,
		   &iwork1, &liwork, &info);

   lwork=(int)floor(work1);if (work1-lwork>0.5) lwork++;
   work=(double *)calloc((size_t)lwork,sizeof(double));
   liwork = iwork1;
   iwork= (int *)calloc((size_t)liwork,sizeof(int));

   /* and the actual call... */
   F77_NAME(dstedc)(&compz,n,
		   d, g, /* lead and su-diag */
		   v, /* eigenvectors on exit */  
                   &ldz, /* dimension of v */
		   work, &lwork,
		   iwork, &liwork, &info);

   if (descending) { /* need to reverse eigenvalues/vectors */
     for (i=0;i<*n/2;i++) { /* reverse the eigenvalues */
       x = d[i]; d[i] = d[*n-i-1];d[*n-i-1] = x;
       dum1 = v + *n * i;dum2 = v + *n * (*n-i-1); /* pointers to heads of cols to exchange */
       for (j=0;j<*n;j++,dum1++,dum2++) { /* work down columns */
         x = *dum1;*dum1 = *dum2;*dum2 = x;
       }
     }
   }

   free(work);
   free(iwork);
   *n=info; /* zero is success */
}

void extreme_eigenvalues(double *A,double *U,double *D,int *n, int *m, int *lm,double *tol) {
/* Faster version of lanczos_spd for calling from R.
   A is n by n symmetric matrix. Let k = m + max(0,lm).
   U is n by k and D is a k-vector.
   m is the number of upper eigenvalues required and lm the number of lower.
   If lm<0 then the m largest magnitude eigenvalues (and their eigenvectors)
   are returned 

   Matrices are stored in R (and LAPACK) format (1 column after another).

   ISSUE: 1. Currently all eigenvectors of Tj are found, although only the next unconverged one
             is really needed. Might be better to be more selective using dstein from LAPACK. 
          2. Basing whole thing on dstevx might be faster
          3. Is random start vector really best? Actually Demmel (1997) suggests using a random vector, 
             to avoid any chance of orthogonality with an eigenvector!
          4. Could use selective orthogonalization, but cost of full orth is only 2nj, while n^2 of method is
             unavoidable, so probably not worth it.  
*/
  //Rprintf("Initializing...\n");
  int biggest=0,f_check,i,k,kk,ok,l,j,vlength=0,ni,pi,converged,incx=1;
  double **q,*v=NULL,bt,xx,yy,*a,*b,*d,*g,*z,*err,*p0,*p1,*zp,*qp,normTj,eps_stop,max_err,alpha=1.0,beta=0.0;
  unsigned long jran=1,ia=106,ic=1283,im=6075; /* simple RNG constants */
  const char uplo='U';
  eps_stop = *tol; 

  if (*lm<0) { biggest=1;*lm=0;} /* get m largest magnitude eigen-values */
  f_check = (*m + *lm)/2; /* how often to get eigen_decomp */
  if (f_check<10) f_check =10;
  kk = (int) floor(*n/10); if (kk<1) kk=1;  
  if (kk<f_check) f_check = kk;

  q=(double **)calloc((size_t)(*n+1),sizeof(double *));

  /* "randomly" initialize first q vector */
  q[0]=(double *)calloc((size_t)*n,sizeof(double));
  b=q[0];bt=0.0;
  for (i=0;i < *n;i++)   /* somewhat randomized start vector!  */
  { 
    jran=(jran*ia+ic) % im;   /* quick and dirty generator to avoid too regular a starting vector */
    xx=(double) jran / (double) im -0.5; 
    b[i]=xx;
    xx=-xx;
    bt+=b[i]*b[i];
  } 
  bt=sqrt(bt); 
  for (i=0;i < *n;i++) {
    b[i]/=bt;
    //Rprintf("b[i]=%f\n",b[i]);
  }

  /* initialise vectors a and b - a is leading diagonal of T, b is sub/super diagonal */
  a=(double *)calloc((size_t) *n,sizeof(double));
  b=(double *)calloc((size_t) *n,sizeof(double));
  g=(double *)calloc((size_t) *n,sizeof(double));
  d=(double *)calloc((size_t) *n,sizeof(double)); 
  z=(double *)calloc((size_t) *n,sizeof(double));
  err=(double *)calloc((size_t) *n,sizeof(double));
  for (i=0;i< *n;i++) err[i]=1e300;
  /* The main loop. Will break out on convergence. */
  for (j=0;j< *n;j++) 
  { /* form z=Aq[j]=A'q[j], the O(n^2) step ...  */
    /*for (Ap=A,zp=z,p0=zp+*n;zp<p0;zp++) 
      for (*zp=0.0,qp=q[j],p1=qp+*n;qp<p1;qp++,Ap++) *zp += *Ap * *qp;*/
    /*  BLAS versions y := alpha*A*x + beta*y, */
    F77_NAME(dsymv)(&uplo,n,&alpha,
		A,n,
		q[j],&incx,
		&beta,z,&incx);
    /* Now form a[j] = q[j]'z.... */
    for (xx=0.0,qp=q[j],p0=qp+*n,zp=z;qp<p0;qp++,zp++) xx += *qp * *zp;
    a[j] = xx;
 
    /* Update z..... */
    if (!j)
    { /* z <- z - a[j]*q[j] */
      for (zp=z,p0=zp+*n,qp=q[j];zp<p0;qp++,zp++) *zp -= xx * *qp;
    } else
    { /* z <- z - a[j]*q[j] - b[j-1]*q[j-1] */
      yy=b[j-1];      
      for (zp=z,p0=zp + *n,qp=q[j],p1=q[j-1];zp<p0;zp++,qp++,p1++) *zp -= xx * *qp + yy * *p1;    
  
      /* Now stabilize by full re-orthogonalization.... */
      
      for (i=0;i<=j;i++) 
      { /* form xx= z'q[i] */
        /*for (xx=0.0,qp=q[i],p0=qp + *n,zp=z;qp<p0;zp++,qp++) xx += *zp * *qp;*/
        xx = -F77_NAME(ddot)(n,z,&incx,q[i],&incx); /* BLAS version */
        /* z <- z - xx*q[i] */
        /*for (qp=q[i],zp=z;qp<p0;qp++,zp++) *zp -= xx * *qp;*/
        F77_NAME(daxpy)(n,&xx,q[i],&incx,z,&incx); /* BLAS version */
      } 
      
      /* exact repeat... */
      for (i=0;i<=j;i++) 
      { /* form xx= z'q[i] */
        /* for (xx=0.0,qp=q[i],p0=qp + *n,zp=z;qp<p0;zp++,qp++) xx += *zp * *qp; */
        xx = -F77_NAME(ddot)(n,z,&incx,q[i],&incx); /* BLAS version */
        /* z <- z - xx*q[i] */
        /* for (qp=q[i],zp=z;qp<p0;qp++,zp++) *zp -= xx * *qp; */
        F77_NAME(daxpy)(n,&xx,q[i],&incx,z,&incx); /* BLAS version */
      } 
      /* ... stabilized!! */
    } /* z update complete */
    

    /* calculate b[j]=||z||.... */
    for (xx=0.0,zp=z,p0=zp+*n;zp<p0;zp++) xx += *zp * *zp;b[j]=sqrt(xx); 
  
    /*if (b[j]==0.0&&j< *n-1) ErrorMessage(_("Lanczos failed"),1);*/ /* Actually this isn't really a failure => rank(A)<l+lm */
    /* get q[j+1]      */
    if (j < *n-1)
    { q[j+1]=(double *)calloc((size_t) *n,sizeof(double));
      for (xx=b[j],qp=q[j+1],p0=qp + *n,zp=z;qp<p0;qp++,zp++) *qp = *zp/xx; /* forming q[j+1]=z/b[j] */
    }

    /* Now get the spectral decomposition of T_j.  */

    if (((j>= *m + *lm)&&(j%f_check==0))||(j == *n-1))   /* no  point doing this too early or too often */
    { for (i=0;i<j+1;i++) d[i]=a[i]; /* copy leading diagonal of T_j */
      for (i=0;i<j;i++) g[i]=b[i]; /* copy sub/super diagonal of T_j */   
      /* set up storage for eigen vectors */
      if (vlength) free(v); /* free up first */
      vlength=j+1; 
      v = (double *)calloc((size_t)vlength*vlength,sizeof(double));
 
      /* obtain eigen values/vectors of T_j in O(j^2) flops */
    
      kk = j + 1;
      mgcv_trisymeig(d,g,v,&kk,1,1);
      /* ... eigenvectors stored one after another in v, d[i] are eigenvalues */

      /* Evaluate ||Tj|| .... */
      normTj=fabs(d[0]);if (fabs(d[j])>normTj) normTj=fabs(d[j]);

      for (k=0;k<j+1;k++) /* calculate error in each eigenvalue d[i] */
      { err[k]=b[j]*v[k * vlength + j]; /* bound on kth e.v. is b[j]* (jth element of kth eigenvector) */
        err[k]=fabs(err[k]);
      }
      /* and check for termination ..... */
      if (j >= *m + *lm)
      { max_err=normTj*eps_stop;
        if (biggest) { /* getting m largest magnitude eigen values */
	  /* only one convergence test is sane here:
             1. Find the *m largest magnitude elements of d. (*lm is 0)
             2. When all these have converged, we are done.
          */   
          pi=ni=0;converged=1;
          while (pi+ni < *m) if (fabs(d[pi])>= fabs(d[j-ni])) { /* include d[pi] in largest set */
              if (err[pi]>max_err) {converged=0;break;} else pi++;
	    } else { /* include d[j-ni] in largest set */
              if (err[ni]>max_err) {converged=0;break;} else ni++;
            }
   
          if (converged) {
            *m = pi;
            *lm = ni;
            j++;break;
          }
        } else /* number of largest and smallest supplied */
        { ok=1;
          for (i=0;i < *m;i++) if (err[i]>max_err) ok=0;
          for (i=j;i > j - *lm;i--) if (err[i]>max_err) ok=0;
          if (ok) 
          { j++;break;}
        }
      }
    }
  }
  /* At this stage, complete construction of the eigen vectors etc. */
  
  /* Do final polishing of Ritz vectors and load va and V..... */
  /*  for (k=0;k < *m;k++) // create any necessary new Ritz vectors 
  { va->V[k]=d[k];
    for (i=0;i<n;i++) 
    { V->M[i][k]=0.0; for (l=0;l<j;l++) V->M[i][k]+=q[l][i]*v[k][l];}
  }*/

  /* assumption that U is zero on entry! */

  for (k=0;k < *m;k++) /* create any necessary new Ritz vectors */
  { D[k]=d[k];
    for (l=0;l<j;l++)
    for (xx=v[l + k * vlength],p0=U + k * *n,p1 = p0 + *n,qp=q[l];p0<p1;p0++,qp++) *p0 += *qp * xx;
  }

  for (k= *m;k < *lm + *m;k++) /* create any necessary new Ritz vectors */
  { kk=j-(*lm + *m - k); /* index for d and v */
    D[k]=d[kk];
    for (l=0;l<j;l++)
    for (xx=v[l + kk * vlength],p0=U + k * *n,p1 = p0 + *n,qp=q[l];p0<p1;p0++,qp++) *p0 += *qp * xx;
  }
 
  /* clean up..... */
  free(a);
  free(b);
  free(g);
  free(d);
  free(z);
  free(err);
  if (vlength) free(v);
  for (i=0;i< *n+1;i++) if (q[i]) free(q[i]);free(q);
    
  //*n = j; /* number of iterations taken */
  
} /* end of Rlanczos */


SEXP selected_eigenvalues(SEXP R_A, SEXP R_n, SEXP R_vl)
{
    double vl=NUMERIC_VALUE(R_vl);
    int n=INTEGER_VALUE(R_n);
    double *A, *lambda, *Z, *U, *D, *work;
    int *iwork, *ifail; 
    
    A = (double *) R_alloc(n * (size_t) n, sizeof(double));
    Memcpy(A, REAL(R_A), (size_t) n * n);
    
    PROTECT(R_A);
    
    int nbiggest=1;
    int nsmallest=-1;
    
    U=(double *) R_alloc((size_t) n, sizeof(double));
    D=(double *) R_alloc(1,sizeof(double));
    double abstol=1e-2;

    extreme_eigenvalues(A,U,D,&n,&nbiggest,&nsmallest,&abstol);

    SEXP R_lambda = PROTECT(allocVector(REALSXP, n));
    lambda = REAL(R_lambda);
        
    SEXP R_Z=PROTECT(allocVector(REALSXP,n * (size_t) n));
    Z=REAL(R_Z); 
    for(int i=0; i<(n*n);i++) Z[i]=0.0;

    
    /* Compute selected eigenvalues */
    abstol=1e-4;
    double vu=1.1*D[0];
	char jobz='V', range='V', uplo='U';
	int one=1;
	int zero=0;
	int info; 
   	int il=-1;    //Ignored
   	int iu=-1;    //Ignored
   	int m;
   	int lwork=-1;
   	double wkopt;
   	iwork=(int *) R_alloc(5 * (size_t) n,sizeof(double));
   	ifail=(int *) R_alloc((size_t) n,sizeof(double));
   	
   	/*
     Query mode
   	*/                         
    F77_NAME(dsyevx)(&jobz,&range,&uplo, &n, A, &n, &vl, &vu, 
                     &il, &iu, &abstol, &m, lambda, Z, &n,                     
                     &wkopt, &lwork, iwork, ifail, &info );
                     
    lwork = (int)wkopt;
    work = (double*) R_alloc((size_t) lwork, sizeof(double));
                     
    /* 
     The real computation 
     */
    F77_NAME(dsyevx)(&jobz,&range,&uplo, &n, A, &n, &vl, &vu, 
                     &il, &iu, &abstol, &m, lambda, Z, &n,                     
                     work, &lwork, iwork, ifail, &info);
    
    SEXP list = PROTECT(allocVector(VECSXP, 2)); 
    
    SET_VECTOR_ELT(list, 0, R_lambda);
    SET_VECTOR_ELT(list, 1, R_Z); 
    
    UNPROTECT(4);   
    
    return(list);
}