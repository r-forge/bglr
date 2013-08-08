
if(FALSE)
{
rm(list=ls())
setwd("/tmp/")
library(BGLR)
}

#This function fits a mixed model
#y=X*Beta + Z*u + e
#where y is the vector of phenotypes
#X is the matrix for FIXED EFFECTS
#Z matrix connection phenotypes and genotypes
#u ~ N(0, varU*A), A a genomic relationship matrix or pedigree
#e ~ N(0, varE*I)

#The model is fitted using the algorithm described in Zhou and Stephens, 2012.
#Genome-wide efficient mixed-model analysis for assiacition studies, Nature genetics, 44(7): 821-824

#Additional arguments
# * method estimation method for the variance components, can be maximum likelihood or restricted maximum likelihood
#          i.e. method=c("ML","RML")
#          lambda_ini=varU/varE
# * d,V: eigenvalues and eigenvectors for G = Z A Z' = U*D*U'


fitMixed=function(y, X=NULL, Z=NULL, A=NULL,d=NULL,U=NULL, BLUE=TRUE, BLUP=TRUE,method="ML", lambda_ini=NULL, tol=1e-4,maxiter=100)
{

	if(!is.vector(y)) stop("y must be a vector\n");
	
	n=length(y)
	
	if(is.null(X))
	{
		warning("X was not set, an intercept was included by default\n",immediate. = TRUE)		
		X=matrix(1,nrow=n,ncol=1)
	}else{
		if(!is.matrix(X)) stop("X must be a matrix\n")
	}
	
	if(is.null(U) & is.null(d))
	{
		if(is.null(A)) stop("A must be a positive semi definite matrix\n")
	
		if(is.null(Z))
		{
			warning("Z was set to Identity matrix\n",immediate. = TRUE)
			G=A
		}else{
			if(!is.matrix(Z)) stop("Z must be a matrix\n")	
			G=Z%*%A%*%t(Z)
		}
		
		out=eigen(G)
		d=out$values
		U=out$vectors
		
	}else{
		if(is.null(U)) stop("You are providing the eigenvectors, but not the eigenvalues");
		if(is.null(d)) stop("You are providing the eigenvalues, but not the eigenvectors");	
	}
	
	v_y=as.vector(crossprod(U,y))
	V=crossprod(U,X)
	
	#Initial value for lambda
	if(is.null(lambda_ini)) lambda_ini=2	

	diff=1
	i=0
	while(diff>tol & i<maxiter)
	{
		i=i+1
		#H=lambda_ini*G+I
		#Hinv=solve(H)
		#Astar=Hinv%*%X%*%solve(t(X)%*%Hinv%*%X)%*%t(X)%*%Hinv
		#P=Hinv-Astar
	
		#dL1=-0.5*sum(diag(Hinv%*%G)) + 0.5*n * t(y)%*%P%*%G%*%P%*%y/(t(y)%*%P%*%y)
		#Trace of Hinv*G, Optimized
		Tr_Hinv=sum(1/(lambda_ini*d+1))
		Tr_Hinv_G=(n-Tr_Hinv)/lambda_ini
    
		#y'*P*y
    	dbar=1/(lambda_ini*d+1)
    	tmp1=t(v_y*dbar)%*%V
    	tmp2=solve(t(V)%*%diag(1/(lambda_ini*d+1))%*%V)
    
    	#Check this, is suboptimal
    	tmp3=sum(v_y^2*dbar^2)
    	tmp4=t(v_y*dbar^2)%*%V%*%tmp2%*%t(tmp1)
    	tmp5=as.vector(t(v_y*dbar)%*%V%*%tmp2%*%t(V))*dbar
    	tmp6=sum(tmp5^2)
    
    	ty_P_y=sum(v_y^2/(lambda_ini*d+1))-tmp1%*%tmp2%*%t(tmp1)
    	ty_P_y=as.numeric(ty_P_y)
    	ty_P_P_y = as.numeric(tmp3-2*tmp4+tmp6)  
    
    	ty_P_G_P_y=(ty_P_y-ty_P_P_y)/lambda_ini
      
		dL1=-0.5*Tr_Hinv_G + 0.5*n * ty_P_G_P_y/(ty_P_y)
		dL1=as.numeric(dL1)
	
		#Check this, probably rounding errors
		  tmp7=sum(v_y^2*dbar^3)
	  	  tmp8=t(tmp5*dbar)%*%V
	      tmp9=as.numeric(tmp8%*%tmp2%*%t(tmp8))
	      tmp10=t(v_y*dbar^2)%*%V
	      tmp11=as.numeric(tmp10%*%tmp2%*%t(tmp10))
	      tmp12=as.numeric(t(v_y*dbar^3)%*%V%*%tmp2%*%t(tmp1))
	      tmp13=t(v_y*dbar^2)%*%V
	      tmp14=as.numeric(tmp8%*%tmp2%*%t(tmp13))
	
	      ty_P_P_P_y = tmp7 - tmp9 - 2*tmp11 - tmp12 + 3*tmp14 #Compare with t(y)%*%P%*%P%*%P%*%y 
		#Up to here
	
		ty_P_G_P_G_P_y = (ty_P_y + ty_P_P_P_y -2 *ty_P_P_y)/lambda_ini^2 
	
	
		#dL2=sum(diag(Hinv%*%G%*%Hinv%*%G))-0.5*n*(2*(t(y)%*%P%*%G%*%P%*%G%*%P%*%y)*(t(y)%*%P%*%y)-(t(y)%*%P%*%G%*%P%*%y)^2)/(t(y%*%P%*%y)^2)
		#Trace Hinv*G*Hinv*G, Optimized
		Tr_Hinv_Hinv=sum(1/(lambda_ini*d+1)^2)
		Tr_Hinv_G_Hinv_G=(n-2*Tr_Hinv+Tr_Hinv_Hinv)/lambda_ini^2
		dL2=Tr_Hinv_G_Hinv_G-0.5*n*(2*(ty_P_G_P_G_P_y)*(ty_P_y)-(ty_P_G_P_y)^2)/(ty_P_y)^2
		dL2=as.numeric(dL2)
	
		cat("iter=",i,"\n")
		#Check the sign
		#Fisher scoring algorithm, Lynch and Walsh, pag. 795.
		lambda_fin = lambda_ini + dL1/dL2
		cat("lambda=varU/varE=",lambda_fin,"\n")
		diff=abs(lambda_ini-lambda_fin)
		lambda_ini=lambda_fin
	}
	
	if(diff<tol & i<maxiter)
	{
		cat("Fisher scoring converged!\n")
		varE=ty_P_y/n
		varU=lambda_fin*varE
		alpha=NULL
		if(BLUE)
		{
			alpha=tmp2%*%t(tmp1)
		}
		dbar=1/(lambda_ini*d+1)
		uhat=NULL
		if(BLUP)
		{
			if(is.null(alpha)) alpha=tmp2%*%t(tmp1)
			#This is wrong Z' is still missing
			#we should have AZ' V^{-1} (y-X*alpha)
			uhat=sweep(U,2L,lambda_ini*d*dbar,FUN="*")%*%t(U)%*%(y-X%*%alpha)
		}
		ans=list(varE=varE, varU=varU,Beta=alpha,u=as.vector(uhat),message="Fisher scoring converged!",method="ML")
	}else{
		cat("Fisher scoring did not converge!\n")
		ans=list(varE=NULL, varU=NULL,Beta=NULL,u=NULL,message="Fisher scoring did not converge!",method="ML")
	}
	
	#return the goodies
	return(ans)
}

if(FALSE)
{
#Examples
#1) Wheat dataset
data(wheat)
A=wheat.A
y=wheat.Y[,1]
ETA=list(list(K=A,model="RKHS"))
set.seed(123)
fm=BGLR(y=y,ETA=ETA,nIter=5000)

 #fm_mixed=fitMixed(y,A=A)
out=eigen(A)
fm_mixed=fitMixed(y=y,A=A,d=out$values,U=out$vectors)

plot(fm$ETA[[1]]$u,fm_mixed$u)
cat("fm_mixed$varE=",fm_mixed$varE,"\n")
cat("fm$varE=",fm$varE,"\n")
cat("fm_mixed$varU=",fm_mixed$varU,"\n")
cat("fm$ETA[[1]]$varU=",fm$ETA[[1]]$varU,"\n")

#mouse dataset
data(mouse)
A=mouse.A
y=mouse.pheno$Obesity.BMI
ETA=list(list(K=A,model="RKHS"))
set.seed(123)

fm=BGLR(y=y,ETA=ETA,nIter=5000)

#fm_mixed=fitMixed(y,A=A)
out=eigen(A)
fm_mixed=fitMixed(y=y,A=A,d=out$values,U=out$vectors)

plot(fm$ETA[[1]]$u,fm_mixed$u)
cat("fm_mixed$varE=",fm_mixed$varE,"\n")
cat("fm$varE=",fm$varE,"\n")
cat("fm_mixed$varU=",fm_mixed$varU,"\n")
cat("fm$ETA[[1]]$varU=",fm$ETA[[1]]$varU,"\n")

}
