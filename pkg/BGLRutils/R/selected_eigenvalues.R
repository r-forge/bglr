extreme_eigenvalues=function(A,tol=1E-6)
{
    B=as.vector(A)
    n=nrow(A)
    m=1
    lm=0
    u=rep(0,n*(m+lm))
    d=rep(0,m+lm)
    tol=1E-6
    
    values=.C("extreme_eigenvalues",as.double(B),as.double(u),as.double(d),as.integer(n),
              as.integer(m),as.integer(lm),as.double(tol))[[3]]
    values
}

selected_eigenvalues=function(x,min_eigenvalue=1e-10)
{
 	x =as.matrix(x)
    n =nrow(x)
    if (!n) stop("0 x 0 matrix")
    if (n != ncol(x)) stop("non-square matrix in 'eigen'")
    n = as.integer(n)
    if(is.na(n)) stop("invalid nrow(x)")
    		
    out=.Call("selected_eigenvalues",x,n,min_eigenvalue)
    out=list(values=out[[1]],
	     vectors=matrix(out[[2]],nrow=n,ncol=n,byrow=FALSE))
    out
}
