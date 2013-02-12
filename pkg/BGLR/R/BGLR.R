## Fixed Effects ##################################################################
setLT.Fixed=function(LT,n,j,y,weights,nLT,saveAt,rmExistingFiles)
{
    # Checked, Gustavo, Nov. 4, 2012
    LT$X=as.matrix(LT$X)
    LT$p=ncol(LT$X)
	
    if(any(is.na(LT$X)))
    { 
	stop(paste(" LP ",j," has NAs in X",sep=""))
    }

    if(nrow(LT$X)!=n)
    {
        stop(paste(" Number of rows of LP ",j,"  not equal to the number of phenotypes.",sep=""))
    }
    
    #This is very inefficient
    #weights
    #for (i in 1:n) 
    #{ 
    #  LT$X[i, ] = weights[i] * LT$X[i, ]  
    #}
    
    #LT$x2=rep(0,LT$p)
    #for(i in 1:LT$p)
    #{ 
    #  LT$x2[i]=sum(LT$X[,i]^2) 
    #}
    
    LT$X=sweep(LT$X,1L,weights,FUN="*")        #weights
    LT$x2=apply(LT$X,2L,function(x) sum(x^2))  #the sum of the square of each of the columns
	
    LT$b=rep(0,LT$p)
    LT$post_b=rep(0,LT$p)
    LT$post_b2=rep(0,LT$p)
    fname=paste(saveAt,"ETA_",j,"_b.dat",sep="")
    LT$NamefileOut=fname; 

    if(rmExistingFiles)
    { 
       unlink(fname) 
    }

    LT$fileOut=file(description=fname,open="w")
    LT$X=as.vector(LT$X)
    LT$varB=1e10
    return(LT)
}

## Gaussian Regression ############################################################
setLT.BRR=function(LT,y,n,j,weights,nLT,R2,saveAt,rmExistingFiles)
{
       
    LT$X=as.matrix(LT$X)
    LT$p=ncol(LT$X)
	
    if(any(is.na(LT$X)))
    { 
      stop(paste(" LP ",j," has NAs in X",sep=""))
    }
    
    if(nrow(LT$X)!=n)
    {
      stop(paste(" Number of rows of LP ",j,"  not equal to the number of phenotypes.",sep=""))
    }
 
    #This is very inefficient
    #weights
    #for (i in 1:n) 
    #{ 
    #   LT$X[i, ] = weights[i]*LT$X[i, ]  
    #}
    
    #LT$x2=rep(0,LT$p)
    #sumMeanXSq=0
    
    #for(i in 1:LT$p)
    #{ 
	#	LT$x2[i]=sum(LT$X[,i]^2) 
	#	sumMeanXSq=sumMeanXSq+mean(LT$X[,i])^2
    #}
    
    
    LT$X=sweep(LT$X,1L,weights,FUN="*")  #weights
    LT$x2=apply(LT$X,2L,function(x) sum(x^2))  #the sum of the square of each of the columns
    sumMeanXSq = sum((apply(LT$X,2L,mean))^2)
    

    if(is.null(LT$df0))
    {
	LT$df0=5
	cat(paste(" Degree of freedom of LP ",j,"  set to default value (",LT$df0,").\n",sep=""))
    }

    if(is.null(LT$R2))
    { 
        LT$R2=R2/nLT
    }

    if(is.null(LT$S0))
    {
        if(LT$df0<2) stop("df0>2 in BRR in order to set S0\n")

	LT$MSx=sum(LT$x2)/n-sumMeanXSq       
	LT$S0=((var(y,na.rm=TRUE)*LT$R2)/(LT$MSx))*(LT$df0-2)  
	cat(paste(" Scale parameter of LP ",j,"  set to default value (",LT$S0,") .\n",sep=""))
    }
	
    LT$b=rep(0,LT$p)
    LT$post_b=rep(0,LT$p)
    LT$post_b2=rep(0,LT$p)
    LT$varB=LT$S0/(LT$df0-2)
    LT$post_varB=0                 
    LT$post_varB2=0
    fname=paste(saveAt,"ETA_",j,"_varB.dat",sep=""); 
    
    if(rmExistingFiles)
    { 
       unlink(fname) 
    }

    LT$NamefileOut=fname
    LT$fileOut=file(description=fname,open="w")
    LT$X=as.vector(LT$X)

    return(LT)
}

## Bayesian LASSO ############################################################
setLT.BL=function(LT,y,n,j,weights,nLT,R2,saveAt,rmExistingFiles)
{
    LT$X=as.matrix(LT$X)
    LT$p=ncol(LT$X)
		
    if(any(is.na(LT$X)))
    {
       stop(paste("LP ",j," has NAs in X",sep=""))
    }

    if(nrow(LT$X)!=n)
    {
        stop(paste(" Number of rows of LP ",j,"  not equal to the number of phenotypes.",sep=""))
    }
 
    #This is very inefficient
    #weights
    #for (i in 1:n) 
    #{ 
    #   LT$X[i, ] = weights[i]*LT$X[i, ]  
    #}
	
    #LT$x2=rep(0,LT$p)
    #sumMeanXSq=0
    #for(i in 1:LT$p)
    #{ 
    #  LT$x2[i]=sum(LT$X[,i]^2)
    #  sumMeanXSq=sumMeanXSq+mean(LT$X[,i])^2
    #}
    
    
    LT$X=sweep(LT$X,1L,weights,FUN="*")  #weights
    LT$x2=apply(LT$X,2L,function(x) sum(x^2))  #the sum of the square of each of the columns
    sumMeanXSq = sum((apply(LT$X,2L,mean))^2)
	
    LT$MSx=sum(LT$x2)/n-sumMeanXSq

    # Prior
    if(is.null(LT$R2))
    { 
       LT$R2=R2/nLT
    }
	
    # Setting default value of lambda
    if(!is.null(LT$lambda))
    { 
    	if(LT$lambda<0)
        {
		stop(" lambda should be positive\n")
	}  
    }
    if(is.null(LT$lambda))
    {
		LT$lambda2=2*(1-R2)/(LT$R2)*LT$MSx
		LT$lambda=sqrt(LT$lambda2)
		cat(paste(" Initial value of lambda in LP ",j," was set to default value (",LT$lambda,")\n",sep=""))
    }else{
	if(LT$lambda<0) stop(" lambda should be positive\n");
        LT$lambda2=LT$lambda^2
    }
	
    # Checking lambda-type
    if(is.null(LT$type))
    {
		LT$type="gamma"
		cat(paste("  By default, the prior density of lambda^2 in the LP ",j,"  was set to gamma.\n",sep=""))
    }else{
		if(!LT$type%in%c("gamma","beta","fixed")) stop(" The prior for lambda^2 should be gamma, beta or a point of mass (i.e., fixed lambda).\n")
    }
    if(LT$type=="gamma")
    {
		if(is.null(LT$shape))
                {
			 LT$shape=1.1
			 cat(paste("  shape parameter in LP ",j," was missing and was set to ",LT$shape,"\n",sep=""))
		}
		
		if(is.null(LT$rate))
                {
			 LT$rate=(LT$shape-1)/LT$lambda2
			 cat(paste("  rate parameter in LP ",j," was missing and was set to ",LT$rate,"\n",sep=""))
		}	
    }
    
    if(LT$type=="beta")
    {
		if(is.null(LT$shape1))
		{
		    LT$shape1=1
		    cat(paste("  shape1 parameter in LP ",j," was missing and was set to ",LT$shape1,"\n",sep=""))
		}
		if(is.null(LT$shape2))
		{
		    LT$shape2=1
		    cat(paste("  shape2 parameter in LP ",j," was missing and was set to ",LT$shape1,"\n",sep=""))
		}
		if(is.null(LT$max))
		{
		    LT$max=10*LT$lambda
		    cat(paste("  max parameter in LP ",j," was missing and was set to ",LT$max,"\n",sep=""))
		}
    }

    LT$b=rep(0,LT$p)
    LT$post_b=rep(0,LT$p)
    LT$post_b2=rep(0,LT$p)
    
    tmp=((var(y,na.rm=TRUE)*R2/nLT)/(LT$MSx))
    LT$tau2=rep(tmp,LT$p)
    LT$post_tau2=0  
    LT$post_lambda=0
    
    fname=paste(saveAt,"ETA_",j,"_lambda.dat",sep="");
    
    if(rmExistingFiles)
    { 
       unlink(fname) 
    }
    
    LT$NamefileOut=fname
    LT$fileOut=file(description=fname,open="w")

    LT$X=as.vector(LT$X)
    return(LT)
}


#Reproducing kernel Hilbert spaces
setLT.RKHS=function(LT,y,n,j,weights,saveAt,R2,nLT,rmExistingFiles)  
{
    if(is.null(LT$V))
    {
	LT$K = as.matrix(LT$K)
        if(nrow(LT$K)!=ncol(LT$K)) stop(paste(" Kermel for linear term ",j, " is not a square matrix\n",sep=""))

	#?# cambie     T = diag(weights)   LT$K = T %*% LT$K %*% T por lo que sigue abajo diag(weights)lo hace muy ineficiente

	for(i in 1:nrow(LT$K))
        {
		#j can not be used as subindex because its value is overwritten
		for(m in i:ncol(LT$K))
                {    
				LT$K[i,m]=LT$K[i,m]*weights[i]*weights[m] ;
				LT$K[m,i]=LT$K[m,i]
		}
	}
        tmp =eigen(LT$K)
        LT$V =tmp$vectors
        LT$d =tmp$values
	rm(tmp)
    }else{
	if(any(weights!=1))
        { 
		cat(paste(" Warning, in LT ",j," Eigen decomposition was provided, and the model involves weights. Note: You should have weighted the kernel before computing eigen(K).\n",sep="")) 
        }
    }
    
    #Defaul value for tolD
    if (is.null(LT$tolD)) 
    {
       LT$tolD = 1e-10
       cat(paste("  Default value of minimum eigenvalue in LP ",j," was set to ",LT$tolD,"\n",sep=""))
    }
    
    #Removing elements whose eigenvalues < tolD
    tmp= LT$d > LT$tolD
    LT$levelsU = sum(tmp)
    LT$d = LT$d[tmp]
    LT$V = LT$V[, tmp]
    
    #Default degrees of freedom and scale parameter
    if (is.null(LT$df0)) 
    {
      LT$df0 = 5
      cat(paste("  default value of df0 in LP ",j," was missing and was set to ",LT$df0,"\n",sep=""))
    }
   
    if(is.null(LT$R2))
    { 
           LT$R2=R2/nLT
    }

    if (is.null(LT$S0)) 
    {
          if(LT$df0<2) stop("df0>2 in RKHS in order to set S0\n");

	  LT$S0=((var(y,na.rm=TRUE)*LT$R2)/(mean(LT$d)))*(LT$df0-2)
          cat(paste("  default value of S0 in LP ",j," was missing and was set to ",LT$S0,"\n",sep=""))
    }
    
    LT$u=rep(0,nrow(LT$K))
    
    LT$varU=LT$S0/(LT$df0-2)
       
    LT$uStar=rep(0, LT$levelsU)
    
    #Output files
    fname=paste(saveAt,"ETA_",j,"_varU.dat",sep="")
    LT$NamefileOut=fname; 

    if(rmExistingFiles)
    { 
       unlink(fname) 
    }

    LT$fileOut=file(description=fname,open="w")
    LT$post_varU=0
    LT$post_uStar = rep(0, LT$levelsU)
    LT$post_u = rep(0, nrow(LT$K))
    LT$post_u2 = rep(0,nrow(LT$K))
    
    #return object
    return(LT)
}

###Bayes B###########################################################################################################################################                    

#Pseudo BayesB 
#See Variable selection for regression models, 
#Lynn Kuo and Bani Mallic, 1998. 
setLT.BayesB=function(LT,y,n,j,weights,saveAt,R2,nLT,rmExistingFiles)
{ 

  #Be sure that your X is a matrix
  LT$X=as.matrix(LT$X) 
  
  LT$p=ncol(LT$X)

  #This is very inefficient
  #weights
  #for (i in 1:n) 
  #{ 
  #   LT$X[i, ] = weights[i]*LT$X[i, ]  
  #}
	
  #LT$x2=rep(0,LT$p)
  #sumMeanXSq=0
  #for(i in 1:LT$p)
  #{ 
  #    LT$x2[i]=sum(LT$X[,i]^2)
  #    sumMeanXSq=sumMeanXSq+mean(LT$X[,i])^2
  #}

  LT$X=sweep(LT$X,1L,weights,FUN="*")  #weights
  LT$x2=apply(LT$X,2L,function(x) sum(x^2))  #the sum of the square of each of the columns
  sumMeanXSq = sum((apply(LT$X,2L,mean))^2)
  LT$MSx=sum(LT$x2)/n-sumMeanXSq

  
  if(any(is.na(LT$X))){ stop(paste("LP ",j," has NAs in X",sep=""))}
  if(nrow(LT$X)!=n){stop(paste("   Number of rows of LP ",j,"  not equal to the number of phenotypes.",sep=""))}
  
  if(is.null(LT$R2))
  {
    LT$R2=R2/nLT
    cat(paste("  R2 in LP ",j," was missing and was set to ",LT$R2,"\n",sep=""))
  }
    
  if(is.null(LT$df0))
  {
    LT$df0= 5
    cat(paste("  DF in LP ",j," was missing and was set to ",LT$df0,"\n",sep=""))
  }

  if(is.null(LT$probIn))
  {
    LT$probIn=0.5
    cat(paste("  probIn in LP ",j," was missing and was set to ",LT$probIn,"\n",sep=""))
  } 

  if(is.null(LT$S0))
  {
     if(LT$df0<2) stop("df0>2 in BayesB in order to set S0\n");

     LT$S0=var(y, na.rm = TRUE)*LT$R2/(LT$MSx)*(LT$df0-2)/LT$probIn
     cat(paste(" Scale paramter in LP ",j," was missing and was set to ",LT$S0,"\n",sep=""))
  }
 
  LT$b=rep(0, LT$p)
  LT$d=rbinom(n = LT$p, size = 1, prob = LT$probIn)

  LT$varB = LT$varB=rep(LT$S0/(LT$df0-2),LT$p)
   
  LT$X=as.vector(LT$X)

  fname=paste(saveAt,"ETA_",j,"_varB.dat",sep="") 

  if(rmExistingFiles)
  { 
      unlink(fname) 
  }

  LT$fileOut=file(description=fname,open="w")
   
  LT$post_varB=0
  LT$post_varB2=0

  LT$post_b=rep(0,LT$p)
  LT$post_b2=rep(0,LT$p)
  
  #return object
  return(LT) 
}



#Bayes C (Habier et al., 2011)
setLT.BayesC=function(LT,y,n,j,weights,saveAt,R2,nLT,rmExistingFiles)
{ 

  #Be sure that your X is a matrix
  LT$X=as.matrix(LT$X) 
  
  LT$p=ncol(LT$X)
  
  #Drop these because it is very inefficient
  #for (i in 1:n) { LT$X[i, ] = weights[i]*LT$X[i, ]  }
	
  #LT$x2=rep(0,LT$p)
  #sumMeanXSq=0
  #for(i in 1:LT$p)
  #{ 
  #    LT$x2[i]=sum(LT$X[,i]^2)
  #    sumMeanXSq=sumMeanXSq+mean(LT$X[,i])^2
  #}

  LT$X=sweep(LT$X,1L,weights,FUN="*")  #weights
  LT$x2=apply(LT$X,2L,function(x) sum(x^2))  #the sum of the square of each of the columns
  sumMeanXSq = sum((apply(LT$X,2L,mean))^2)

  LT$MSx=sum(LT$x2)/n-sumMeanXSq

  
  if(any(is.na(LT$X))){ stop(paste("LP ",j," has NAs in X",sep=""))}
  if(nrow(LT$X)!=n){stop(paste("   Number of rows of LP ",j,"  not equal to the number of phenotypes.",sep=""))}
  
  if(is.null(LT$R2))
  {
    LT$R2=R2/nLT
    cat(paste("  R2 in LP ",j," was missing and was set to ",LT$R2,"\n",sep=""))
  }
    
  if(is.null(LT$df0))
  {
    LT$df0 = 5
    cat(paste("  DF in LP ",j," was missing and was set to ",LT$df0,"\n",sep=""))
  }

  if(is.null(LT$probIn))
  {
    LT$probIn=0.5
    cat(paste("  probIn in LP ",j," was missing and was set to ",LT$probIn,"\n",sep=""))
  }
  
  if(is.null(LT$counts))
  {
    LT$counts=10
    cat(paste("  Counts in LP ",j," was missing and was set to ",LT$counts,"\n",sep=""))
  }  
  
  LT$countsIn=LT$counts * LT$probIn
  LT$countsOut=LT$counts - LT$countsIn
  
  if(is.null(LT$S0)){

    if(LT$df0<2) stop("df0>2 in BayesC in order to set S0\n")

    LT$S0 = var(y, na.rm = TRUE)*LT$R2/(LT$MSx)*(LT$df0-2)/LT$probIn
    cat(paste("  Scale parameter in LP ",j," was missing and was set to ",LT$counts,"\n",sep=""))
  }
    
  LT$b=rep(0, LT$p)
  LT$d=rbinom(n = LT$p, size = 1, prob = LT$probIn)

  LT$varB = LT$S0
   
  LT$X=as.vector(LT$X)

  fname=paste(saveAt,"ETA_",j,"_parBayesC.dat",sep="") 

  if(rmExistingFiles)
  { 
    unlink(fname) 
  }

  LT$fileOut=file(description=fname,open="w")
   
  LT$post_varB=0
  LT$post_varB2=0
  LT$post_d=0
  LT$post_probIn=0
  LT$post_b=rep(0,LT$p)
  LT$post_b2=rep(0,LT$p)
  
  #return object
  return(LT)
}

#Bayes A, Mewissen et al. (2001).
#Prediction of Total Genetic Value Using Genome-Wide Dense Marker Maps
#Genetics 157: 1819-1829

setLT.BayesA=function(LT,y,n,j,weights,saveAt,R2,nLT,rmExistingFiles)
{  
  LT$X=as.matrix(LT$X)
  LT$p=ncol(LT$X)

  #Drop this because is very inefficient
  #for (i in 1:n) { LT$X[i, ] = weights[i]*LT$X[i, ]  }
	
  #LT$x2=rep(0,LT$p)
  #sumMeanXSq=0
  #for(i in 1:LT$p)
  #{ 
  #    LT$x2[i]=sum(LT$X[,i]^2)
  #    sumMeanXSq=sumMeanXSq+mean(LT$X[,i])^2
  #}
  
  LT$X=sweep(LT$X,1L,weights,FUN="*")  #weights
  LT$x2=apply(LT$X,2L,function(x) sum(x^2))  #the sum of the square of each of the columns
  sumMeanXSq = sum((apply(LT$X,2L,mean))^2)
  LT$MSx=sum(LT$x2)/n-sumMeanXSq

  if(is.null(LT$df0))
  {
     LT$df0 = 5  
     cat(paste("  DF in LP ",j," was missing and was set to ",LT$df0,".\n",sep=""))
  }
  if(is.null(LT$R2)){
  	LT$R2=R2/nLT
    cat(paste("  R2 in LP ",j," was missing and was set to ",LT$R2,"\n",sep=""))
  }
  if(is.null(LT$S0))
  {
     if(LT$df0<2) stop("df0>2 in BayesA in order to set S0\n")
     LT$S0 = var(y, na.rm = TRUE)*LT$R2/(LT$MSx)*(LT$df0-2)
     cat(paste(" Scale paramter in LP ",j," was missing and was set to ",LT$S0,"\n",sep=""))
  }
  # Improvement: Treat Scale as random, assign a gamma density 
  
  LT$b=rep(0,LT$p)
   
  LT$varB=rep(LT$S0/(LT$df0-2),LT$p)
  
  # No files for now, add one file when S0 is treated as random.
    
  LT$X=as.vector(LT$X)
  
  LT$post_varB=0
  LT$post_varB2=0
  
  LT$post_b=rep(0,LT$p)
  LT$post_b2=rep(0,LT$p)
  
  #return object
  return(LT)
}

##################################################################################################
welcome=function()
{
  cat("\n");
  cat("#--------------------------------------------------------------------#\n");
  cat("#        _\\\\|//_                                                     #\n");
  cat("#       (` o-o ')      BGLR v1.0                                     #\n");
  cat("#------ooO-(_)-Ooo---------------------------------------------------#\n");
  cat("#                      Bayesian Generalized Linear Regression        #\n");
  cat("#                      Gustavo de los Campos, gdeloscampos@gmail.com #\n");
  cat("#    .oooO     Oooo.   Paulino Perez, perpdgo@colpos.mx              #\n");
  cat("#    (   )     (   )   December, 2012                                #\n");
  cat("#_____\\ (_______) /_________________________________________________ #\n");
  cat("#      \\_)     (_/                                                   #\n");
  cat("#                                                                    #\n");
  cat("#------------------------------------------------------------------- #\n");
  cat("\n");
}
##################################################################################################

##################################################################################################
#The density of a scaled inverted chi-squered distribution
dScaledInvChisq=function (x, df, S)
{
    tmp = dchisq(S/x, df = df)/(x^2)
    return(tmp)
}


##################################################################################################
#The density function for lambda

dLambda=function (rate, shape, lambda) 
{
    tmp = dgamma(x = I(lambda^2), rate = rate, shape = shape) * 2 * lambda
    return(tmp)
}

##################################################################################################
#Metropolis sampler for lambda in the Bayesian LASSO

metropLambda=function (tau2, lambda, shape1 = 1.2, shape2 = 1.2, max = 200, ncp = 0)
{
    lambda2 = lambda^2
    l2_new = rgamma(rate = sum(tau2)/2, shape = length(tau2),
        n = 1)
    l_new = sqrt(l2_new)
    logP_old = sum(dexp(x = tau2, log = TRUE, rate = (lambda2/2))) +
        dbeta(x = lambda/max, log = TRUE, shape1 = shape1, shape2 = shape2) -
        dgamma(shape = sum(tau2)/2, rate = length(tau2), x = (2/lambda2),
            log = TRUE)
    logP_new = sum(dexp(x = tau2, log = TRUE, rate = (l2_new/2))) +
        dbeta(x = l_new/max, log = TRUE, shape1 = shape1, shape2 = shape2) -
        dgamma(shape = sum(tau2)/2, rate = length(tau2), x = (2/l2_new),
            log = TRUE)
    accept = (logP_new - logP_old) > log(runif(1))
    if (accept) {
        lambda = l_new
    }
    return(lambda)
}

##################################################################################################
#Startup function
.onAttach = function(library, pkg)
{
  Rv = R.Version()
  if(!exists("getRversion", baseenv()) || (getRversion() < "2.12"))
    stop("This package requires R 2.12.1 or later")
  assign(".BGLR.home", file.path(library, pkg),
         pos=match("package:BGLR", search()))
  BGLR.version = "1.0 (2012-12-21)"
  assign(".BGLR.version", BGLR.version, pos=match("package:BGLR", search()))
  if(interactive())
  {
    packageStartupMessage(paste("Package 'BGLR', ", BGLR.version, ". ",sep=""),appendLF=TRUE)
    packageStartupMessage("Type 'help(BGLR)' for summary information",appendLF=TRUE)
  }
  invisible()
}


##################################################################################################
#rtrun draws from a truncated univariate normal distribution using the inverse CDF algorithm
#Arguments:
#mu: mean
#sigma: standard deviation
#a: lower bound
#b: upper bound
#NOTES: 1) This routine was taken from bayesm package, December 18, 2012
#       2) The inputs are not checked
rtrun=function (mu, sigma, a, b) 
{
    FA = pnorm(((a - mu)/sigma))
    FB = pnorm(((b - mu)/sigma))
    return(mu + sigma * qnorm(runif(length(mu)) * (FB - FA) + FA))
}

#Extract the values of z such that y[i]=j
#z,y vectors, j integer
extract=function(z,y,j) subset(as.data.frame(z,y),subset=(y==j))

#Random number from inverse gaussian distribution
#Taken from SuppDists package, version 1.1-8, by Bob Wheeler <bwheelerg@gmail.com>
#arguments:
#n is sample size
#nu vector real and non-negative parameter -- the Wald distribution results when nu=1
#lambda vector real and non-negative parameter
#
#Density function
#f(x,nu,lambda)=sqrt[lambda/(2 pi x^3)]exp[-lambda(x-nu)^2/(2 x nu^2)]
rinvGauss=function (n, nu, lambda)
{
    n =if (length(n) > 1)
        length(n)
    else n
    N = max(length(nu), length(lambda))
    nu = rep(nu, length.out = N)
    lambda = rep(lambda, length.out = N)
    .C("rinvGaussR", as.double(nu), as.double(lambda), as.integer(n), as.integer(N), value = double(n))$value
}

#log-likelihood for ordinal data
#y: response vector
#predicted response vector, yHat=X%*%beta
#threshold
loglik_ordinal=function(y,yHat,threshold)
{
	sum=0
        n=length(y)
	for(i in 1:n)
        {
           sum=sum + log(pnorm(threshold[y[i] + 1]-yHat[i])-pnorm(threshold[y[i]]-yHat[i]))
        }
        return(sum)
}

##################################################################################################

#Arguments:
#y: data vector, NAs allowed
#response_type: can be "gaussian", "Bernoulli" or "ordinal",
#ETA: The linear predictor
#nIter: Number of MCMC iterations
#burnIn: burnIn
#thin: thin
#saveAt: string, where to save the information
#Se: Scale parameter for the prior for varE
#dfe: Degrees of freedom for the prior for varE
#weights: 
#R2
#ncores: number of cores used in computations. If ncores=1 then it will 
#        use OpenMP to perform the computations in UNIX like systems.
#Note: The function was designed to work with gaussian responses, some changes were made to deal binary and ordinal responses


#To add new method: 
#(a) create setLT, 
#(b) add it to the switch statement,
#(c) add code to update parameters in the Gibbs sampler, 
#(d) add code to save samples
#(e) add code to compute posterior means
#(f) Test:
#(f1) Test simple example without hyper-paramaeters, evaluate how  
#        default values were set
	#(f2)  Check posterior means and files
	#(f3)  Test simple example with some hyper-parameters give and 
	#         some set by default
#(f4) Run an example with a few missing values, compare to BLR 
#       example, check: (i) residual variance, (ii) plot of effects, (iii) plot 
#        of predictions in trn, (iv) plot of prediction in tst.


BGLR=function (y, response_type = "gaussian", a = NULL, b = NULL, 
    ETA = NULL, nIter = 1500, burnIn = 500, thin = 5, saveAt = "", 
    S0 = NULL, df0 = 5, R2 = 0.5, minAbsBeta = 1e-09, weights = NULL, 
    ncores = 1, verbose = TRUE, rmExistingFiles = TRUE) 
{
    welcome()

    if (!(response_type %in% c("gaussian", "Bernoulli", "ordinal")))  stop(" Only gaussian, Bernoulli and ordinal responses are allowed\n")

    if (saveAt == "") {
        saveAt = paste(getwd(), "/", sep = "")
    }

    a = as.vector(a)
    b = as.vector(b)
    y = as.vector(y)
    n = length(y)

    if (is.null(weights)) 
    {
        weights = rep(1, n)
    }

    sumW2 = sum(weights^2)
    nSums = 0

    whichNa = which(is.na(y))
    nNa = length(whichNa)

    Censored = FALSE

    if (response_type == "gaussian") 
    {
        if ((!is.null(a)) | (!is.null(b))) 
        {
            Censored = TRUE
            if ((length(a) != n) | (length(b) != n)) stop(" y, a and b must have the same dimension\n")
            if (any(weights != 1)) stop(" Weights are only implemented for Gausian uncensored responses\n")
        }
        mu = weighted.mean(x = y, w = weights, na.rm = TRUE)
    }
    post_mu = 0
    post_mu2 = 0

    fname = paste(saveAt, "mu.dat", sep = "")
    if (rmExistingFiles) 
    {
        unlink(fname)
    }
    else {
        cat(" Note: samples will be appended to existing files. \n")
    }

    fileOutMu = file(description = fname, open = "w")

    #Do a small change, for working with the data augmentation algorithm, Albert & Chib, 1993.
    if (response_type == "Bernoulli") {
        cat(" Prior for residual is not necessary, if you provided it, it will be ignored\n")
        if (any(weights != 1)) stop(" Weights are not supported \n")
       
        z = y
        
        phat = mean(z, na.rm = TRUE)
        mu = qnorm(sd = 1, mean = 0, p = phat)
        whichZero = which(z == 0)
        whichOne = which(z == 1)
        nzero = length(whichZero)
        none = length(whichOne)
        y[whichZero] = rtrun(mu = rep(mu, nzero), sigma = rep(1, nzero), a = rep(-Inf, nzero), b = rep(0, nzero))
        y[whichOne] = rtrun(mu = rep(mu, none), sigma = rep(1, none), a = rep(0, none), b = rep(Inf, none))
        if (nNa > 0) y[whichNa] = rnorm(n = nNa, mean = mu, sd = 1)
    }

    if (response_type == "ordinal") {
        cat(" Prior for residual is not necessary, if you provided it, it will be ignored\n")
        if (any(weights != 1)) stop(" Weights are not supported \n")
        if (nNa > 0) stop(" Missing values are not implemented yet for ordinary traits\n")
        
        z = y

        #initialize cut-points
        tmp = table(z)
        nclass = length(names(tmp))
        if (nclass <= 2) stop(paste(" Data vector y has only ", nclass, " differente values, it should have at least 3 different values\n"))
        threshold = c(-Inf, seq(-3.5, 3.5, length.out = nclass - 1), Inf)

        #initialize y based on the z-values
        y[z == 1] = -4
        y[z == nclass] = 4
        for (m in 2:(nclass - 1)) {
            y[z == m] = 0.5 * (threshold[m] + threshold[m + 1])
        }
	
	#mu, 
        #if nclass-1 of the thresholds can move freely, then the intercept must be removed from the model
        mu=0
        
        #posterior for thresholds
        post_threshold = 0
        post_threshold2 = 0
    }

    post_logLik = 0

    # yStar & yHat
    yStar = y * weights
    yHat = mu * weights
    
    if (nNa > 0) {
        yStar[whichNa] = yHat[whichNa]
    }

    post_yHat = rep(0, n)
    post_yHat2 = rep(0, n)

    # residual and residual variance
    e = (yStar - yHat)

    varE = var(e, na.rm = TRUE) * (1 - R2)
    sdE = sqrt(varE)

    if (is.null(S0)) {
        S0 = varE * (df0 - 2)
    }

    post_varE = 0
    post_varE2 = 0

    fname = paste(saveAt, "varE.dat", sep = "")

    if (rmExistingFiles) {
        unlink(fname)
    }

    fileOutVarE = file(description = fname, open = "w")

    nLT = ifelse(is.null(ETA), 0, length(ETA))

    if (nLT > 0) {
        for (i in 1:nLT) {
            if (!(ETA[[i]]$model %in% c("FIXED", "BRR", "BL", "BayesA", "BayesB","BayesC", "RKHS"))) 
            {
                stop(paste(" Error in ETA[[", i, "]]", " model ", ETA[[i]]$model, " not implemented (note: evaluation is case sensitive).", sep = ""))
            }
            ETA[[i]] = switch(ETA[[i]]$model, 
			      FIXED = setLT.Fixed(LT = ETA[[i]],  n = n, j = i, weights = weights, y = y, nLT = nLT, saveAt = saveAt, rmExistingFiles = rmExistingFiles), 
                              BRR = setLT.BRR(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles), 
                              BL = setLT.BL(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles), 
                              RKHS = setLT.RKHS(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles), 
                              BayesC = setLT.BayesC(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles), 
                              BayesA = setLT.BayesA(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles),
                              BayesB = setLT.BayesB(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles)
                              #BEN=setLT.BEN(LT=ETA[[i]],n=n,j=i,weights=weights,y=y,nLT=nLT,R2,saveAt=saveAt)
                              )
        }
    }
    # Gibbs sampler
    time = proc.time()[3]
    for (i in 1:nIter) {
        # intercept
        e = e + weights * mu
        rhs = sum(weights * e)/varE
        C = sumW2/varE
        sol = rhs/C
        mu = rnorm(n = 1, sd = sqrt(1/C)) + sol

        if (response_type == "ordinal") {
            mu = 0
        }

        e = e - weights * mu
        
        #deltaSS and deltadf for updating varE
        deltaSS = 0
        deltadf = 0

        if (nLT > 0) {
            for (j in 1:nLT) {
                ## Fixed effects ####################################################################
                if (ETA[[j]]$model == "FIXED") {
                  varBj = rep(ETA[[j]]$varB, ETA[[j]]$p)
                  ans = .Call("sample_beta", n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, 
                                             e, varBj, varE, minAbsBeta, ncores)
                  ETA[[j]]$b = ans[[1]]
                  e = ans[[2]]
                }#End of fixed effects

                ## Ridge Regression ##################################################################
                if (ETA[[j]]$model == "BRR") {
                  varBj = rep(ETA[[j]]$varB, ETA[[j]]$p)
                  ans = .Call("sample_beta", n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, 
                                             e, varBj, varE, minAbsBeta, ncores)
                  ETA[[j]]$b = ans[[1]]
                  e = ans[[2]]

                  DF = ETA[[j]]$df0 + ETA[[j]]$p
                  SS = sum(ETA[[j]]$b^2) + ETA[[j]]$S0
                  ETA[[j]]$varB = SS/rchisq(df = DF, n = 1)
                }# END BRR

                ## Bayesian LASSO ####################################################################
                if (ETA[[j]]$model == "BL") {
                  varBj = ETA[[j]]$tau2 * varE
                  ans = .Call("sample_beta", n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, 
                                             e, varBj, varE, minAbsBeta, ncores)
                  ETA[[j]]$b = ans[[1]]
                  e = ans[[2]]

                  nu = sqrt(varE) * ETA[[j]]$lambda/abs(ETA[[j]]$b)
                  tmp = NULL
                  try(tmp <- rinvGauss(n = ETA[[j]]$p, nu = nu, lambda = ETA[[j]]$lambda2))
                  if (!is.null(tmp)) {
                    if (!any(is.na(sqrt(tmp)))) {
                      ETA[[j]]$tau2 = 1/tmp
                    }
                    else {
                      warning(paste("tau2 was not updated in iteration",i, "due to numeric problems with beta\n",sep=" "),immediate. = TRUE)
                    }
                  }
                  else {
                    warning(paste("tau2 was not updated  in iteration",i,"due to numeric problems with beta\n",sep=" "),immediate. = TRUE)
                  }

                  #Update lambda 
                  if (ETA[[j]]$type == "gamma") {
                    rate = sum(ETA[[j]]$tau2)/2 + ETA[[j]]$rate
                    shape = ETA[[j]]$p + ETA[[j]]$shape
                    ETA[[j]]$lambda2 = rgamma(rate = rate, shape = shape, n = 1)
                    if (!is.na(ETA[[j]]$lambda2)) {
                      ETA[[j]]$lambda = sqrt(ETA[[j]]$lambda2)
                    }
                    else {
                      warning(paste("lambda was not updated in iteration",i, "due to numeric problems with beta\n",sep=" "),immediate. = TRUE)
                    }
                  }

                  if (ETA[[j]]$type == "beta") {
                    ETA[[j]]$lambda = metropLambda(tau2 = ETA[[j]]$tau2, 
                                                   lambda = ETA[[j]]$lambda, shape1 = ETA[[j]]$shape1, shape2 = ETA[[j]]$shape2, 
                                                   max = ETA[[j]]$max)
                    ETA[[j]]$lambda2 = ETA[[j]]$lambda^2
                  }

                  deltaSS = deltaSS + sum((ETA[[j]]$b/sqrt(ETA[[j]]$tau2))^2)
                  deltadf = deltadf + ETA[[j]]$p
                }#END BL

                ## RKHS ####################################################################
                if (ETA[[j]]$model == "RKHS") {
                  #error
                  e = e + ETA[[j]]$u
                  rhs = crossprod(ETA[[j]]$V, e)/varE
                  varU = ETA[[j]]$varU * ETA[[j]]$d
                  C = as.numeric(1/varU + 1/varE)
                  SD = 1/sqrt(C)
                  sol = rhs/C
                  tmp = rnorm(n = ETA[[j]]$levelsU, mean = sol, sd = SD)
                  ETA[[j]]$uStar = tmp
                  ETA[[j]]$u = as.vector(ETA[[j]]$V %*% tmp)
		  
                  #update error
                  e = e - ETA[[j]]$u
                   
                  #update the variance
                  tmp = ETA[[j]]$uStar/sqrt(ETA[[j]]$d)
                  SS = as.numeric(crossprod(tmp)) + ETA[[j]]$S0
                  DF = ETA[[j]]$levelsU + ETA[[j]]$df0
                  ETA[[j]]$varU = SS/rchisq(n = 1, df = DF)
                }#END RKHS

                ## Bayes C ################################################################################
                if (ETA[[j]]$model == "BayesC") {
                  #Update marker effects
                  mrkIn = ETA[[j]]$d == 1
                  pIn = sum(mrkIn)
                  if (pIn > 0) {
                    x = .Call("extract_column", which(mrkIn), n, ETA[[j]]$X)

                    ans = .Call("sample_beta", n, pIn, x, ETA[[j]]$x2[mrkIn], ETA[[j]]$b[mrkIn], 
                                               e, rep(ETA[[j]]$varB, pIn), varE, minAbsBeta, ncores)
                    ETA[[j]]$b[mrkIn] = ans[[1]]
                    e = ans[[2]]
                  }

                  if ((ETA[[j]]$p - pIn) > 0) {
                    ETA[[j]]$b[(!mrkIn)] = rnorm(n = (ETA[[j]]$p - pIn), sd = sqrt(ETA[[j]]$varB))
                  }

                  #Update indicator variables #?# discuss this ##
                  ans = .Call("d_e", ETA[[j]]$p, n, ETA[[j]]$X, ETA[[j]]$d, ETA[[j]]$b, e, varE, ETA[[j]]$probIn, ncores)

                  ETA[[j]]$d = ans[[1]]
                  e = ans[[2]]

                  #Update the variance of marker effects
                  SS = sum(ETA[[j]]$b^2) + ETA[[j]]$S0
                  DF = ETA[[j]]$df0 + ETA[[j]]$p
                  ETA[[j]]$varB = SS/rchisq(df = DF, n = 1)
                  mrkIn = sum(ETA[[j]]$d)
                  ETA[[j]]$probIn = rbeta(shape1 = (mrkIn + ETA[[j]]$countsIn + 1), 
                                          shape2 = (ETA[[j]]$p - mrkIn + ETA[[j]]$countsOut + 1), n = 1)
                }#End BayesC

                ## BayesA ##############################################################################
                if (ETA[[j]]$model == "BayesA") {
                  varBj = ETA[[j]]$varB
                  ans = .Call("sample_beta", n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, 
                                             e, varBj, varE, minAbsBeta, ncores)
                  ETA[[j]]$b = ans[[1]]
                  e = ans[[2]]
                   
                  #Update variances
                  SS = ETA[[j]]$S0 + ETA[[j]]$b^2
                  DF = ETA[[j]]$df0 + 1
                  ETA[[j]]$varB = SS/rchisq(n = ETA[[j]]$p, df = DF)
                }#End BayesA
                
                #BayesB
                if(ETA[[j]]$model=="BayesB")
                {
                      #Update marker effects
                      mrkIn=ETA[[j]]$d==1
                      pIn=sum(mrkIn)
                      if(pIn>0)
                      {
                        x=.Call("extract_column",which(mrkIn),n, ETA[[j]]$X)

                        ans = .Call("sample_beta", n, pIn, x, ETA[[j]]$x2[mrkIn], ETA[[j]]$b[mrkIn],
                                       e, ETA[[j]]$varB, varE, minAbsBeta,ncores)
                        ETA[[j]]$b[mrkIn] = ans[[1]]
                        e = ans[[2]]
                      }

                      if((ETA[[j]]$p-pIn)>0)
                      {
                          ETA[[j]]$b[(!mrkIn)]=rnorm(n=(ETA[[j]]$p-pIn),sd=sqrt(ETA[[j]]$varB))
                      }

                      #Update indicator variables
                      ans=.Call("d_e",ETA[[j]]$p,n,ETA[[j]]$X,ETA[[j]]$d,ETA[[j]]$b,e,varE,ETA[[j]]$probIn,ncores)
                      ETA[[j]]$d=ans[[1]]
                      e=ans[[2]]

                      #Update the variance component associated with the markers
                      SS = ETA[[j]]$b^2 + ETA[[j]]$S0
                      DF = ETA[[j]]$df0+1
                      ETA[[j]]$varB = SS/rchisq(df=DF, n = ETA[[j]]$p)
                }#End Bayes B

            }#Loop for
        }#nLT
        
        # yHat & missing values
        # Now only for gaussian responses
        yHat = yStar - e
        if (nNa > 0) {
            if (Censored) {
                yStar[whichNa] = rtrun(mu = yHat[whichNa], a = a[whichNa], b = b[whichNa], sigma = sdE)
            }
            else {
                yStar[whichNa] = yHat[whichNa] + rnorm(n = nNa, sd = sdE)
            }
            e[whichNa] = yStar[whichNa] - yHat[whichNa]
        }

        # residual variance for Gaussian
        # and Bernoulli families
        if (response_type == "gaussian") {
            SS = sum(e * e) + S0 + deltaSS
            DF = n + df0 + deltadf
            varE = SS/rchisq(n = 1, df = DF)
            sdE = sqrt(varE)
        }

        if (response_type == "Bernoulli") {
            varE = 1
            sdE = 1

            #Update yStar, this is the latent variable
            yStar[whichOne] = rtrun(mu = yHat[whichOne], sigma = rep(1, none), a = rep(0, none), b = rep(Inf, none))
            yStar[whichZero] = rtrun(mu = yHat[whichZero], sigma = rep(1, nzero), a = rep(-Inf, nzero), b = rep(0, nzero))

            #Update error
            e = yStar - yHat
        }

        if (response_type == "ordinal") {
            varE = 1
            sdE = 1
            
            #Update yStar, this is the latent variable
            for (m in 1:n) {
                yStar[m] = rtrun(mu = yHat[m], sigma = 1, a = threshold[z[m]], b = threshold[z[m] + 1])
            }

            #Update thresholds
            for (m in 2:nclass) {
                lo = max(max(extract(yStar, z, m - 1)), threshold[m - 1])
                hi = min(min(extract(yStar, z, m)), threshold[m + 1])
                threshold[m] = runif(1, lo, hi)
            }

            #Update error
            e = yStar - yHat
        }

        # Saving samples and computing running means
        if ((i%%thin == 0)) {
            if (nLT > 0) {
                for (j in 1:nLT) {

                  if (ETA[[j]]$model == "FIXED") {
                    write(ETA[[j]]$b, file = ETA[[j]]$fileOut, append = TRUE)
                  }

                  if (ETA[[j]]$model == "BRR") {
                    write(ETA[[j]]$varB, file = ETA[[j]]$fileOut, append = TRUE)
                  }

                  if (ETA[[j]]$model == "BL") {
                    write(ETA[[j]]$lambda, file = ETA[[j]]$fileOut, append = TRUE)
                  }

                  if (ETA[[j]]$model == "RKHS") {
                    write(ETA[[j]]$varU, file = ETA[[j]]$fileOut, append = TRUE)
                  }

                  if (ETA[[j]]$model == "BayesC") {
                    tmp = c(ETA[[j]]$probIn, ETA[[j]]$varB)
                    write(tmp, ncolumns = 2, file = ETA[[j]]$fileOut, append = TRUE)
                  }

                  if (ETA[[j]]$model == "BayesA") {
                    # Nothing here for now
                  }
                  
                  if(ETA[[j]]$model=="BayesB")
                  {
                        write(round(ETA[[j]]$varB,6),file=ETA[[j]]$fileOut,append=TRUE,ncolumns=length(ETA[[j]]$varB),sep=" ")
                  }
                }
            }

            #Output files
            write(x = mu, file = fileOutMu, append = TRUE)
            write(x = varE, file = fileOutVarE, append = TRUE)
            if (i > burnIn) {
                nSums = nSums + 1
                k = (nSums - 1)/(nSums)
                if (nLT > 0) {
                  for (j in 1:nLT) {

                    if (ETA[[j]]$model == "FIXED") {
                      ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
                    }

                    if (ETA[[j]]$model == "BRR") {
                      ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
                      ETA[[j]]$post_varB = ETA[[j]]$post_varB * k + (ETA[[j]]$varB)/nSums
                      ETA[[j]]$post_varB2 = ETA[[j]]$post_varB2 * k + (ETA[[j]]$varB^2)/nSums
                    }

                    if (ETA[[j]]$model == "BL") {
                      ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
                      ETA[[j]]$post_tau2 = ETA[[j]]$post_tau2 * k + (ETA[[j]]$tau2)/nSums
                      ETA[[j]]$post_lambda = ETA[[j]]$post_lambda * k + (ETA[[j]]$lambda)/nSums
                    }

                    if (ETA[[j]]$model == "RKHS") {
                      ETA[[j]]$post_varU = ETA[[j]]$post_varU * k + ETA[[j]]$varU/nSums
                      ETA[[j]]$post_uStar = ETA[[j]]$post_uStar * k + ETA[[j]]$uStar/nSums
                      ETA[[j]]$post_u = ETA[[j]]$post_u * k + ETA[[j]]$u/nSums
                      ETA[[j]]$post_u2 = ETA[[j]]$post_u2 * k + (ETA[[j]]$u^2)/nSums
                    }

                    if (ETA[[j]]$model == "BayesC") {
                      ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
                      ETA[[j]]$post_varB = ETA[[j]]$post_varB * k + (ETA[[j]]$varB)/nSums
                      ETA[[j]]$post_varB2 = ETA[[j]]$post_varB2 * k + (ETA[[j]]$varB^2)/nSums
                      ETA[[j]]$post_d = ETA[[j]]$post_d * k + (ETA[[j]]$d)/nSums
                      ETA[[j]]$post_probIn = ETA[[j]]$post_probIn * k + (ETA[[j]]$probIn)/nSums
                    }

                    if (ETA[[j]]$model == "BayesA") {
                      ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
                      ETA[[j]]$post_varB = ETA[[j]]$post_varB * k + (ETA[[j]]$varB)/nSums
                      ETA[[j]]$post_varB2 = ETA[[j]]$post_varB2 * k + (ETA[[j]]$varB^2)/nSums
                    }

                    if(ETA[[j]]$model=="BayesB")
                    {
                        ETA[[j]]$post_b=ETA[[j]]$post_b*k+ETA[[j]]$b/nSums
                        ETA[[j]]$post_b2=ETA[[j]]$post_b2*k+(ETA[[j]]$b^2)/nSums
                        ETA[[j]]$post_varB=ETA[[j]]$post_varB*k+(ETA[[j]]$varB)/nSums
                        ETA[[j]]$post_varB2=ETA[[j]]$post_varB2*k+(ETA[[j]]$varB2^2)/nSums
                    }
                  }
                }

                post_mu = post_mu * k + mu/nSums
                post_mu2 = post_mu2 * k + (mu^2)/nSums

                post_yHat = post_yHat * k + yHat/nSums
                post_yHat2 = post_yHat2 * k + (yHat^2)/nSums

                post_varE = post_varE * k + varE/nSums
                post_varE2 = post_varE2 * k + (varE^2)/nSums

                if (response_type == "ordinal") {
                  post_threshold = post_threshold * k + threshold/nSums
                  post_threshold2 = post_threshold2 * k + (threshold^2)/nSums
                 
                  #Check this
                  logLik=loglik_ordinal(z,yHat,threshold)
		  
                }

                if(response_type == "gaussian") {
                  tmpE = e/weights
                  tmpSD = sqrt(varE)/weights

                  if (nNa > 0) {
                    tmpE = tmpE[-whichNa]
                    tmpSD = tmpSD[-whichNa]
                  }
                  
	          logLik = sum(dnorm(tmpE, sd = tmpSD, log = TRUE))

                }#end gaussian

		if(response_type== "Bernoulli")
                {
		    #Be careful, the vector with the original response here is z
		    if (nNa > 0) {
                        pSuccess = pnorm(mean = 0, sd = 1, q = yHat[-whichNa])
                        tmp = z[-whichNa]
                        logLik = sum(log(ifelse(z[-whichNa] == 0, (1 - pSuccess), pSuccess)))
                  }
                  else {
                        pSuccess = pnorm(mean = 0, sd = 1, q = yHat)
                        tmp = z[-whichNa]
                        logLik = sum(log(ifelse(z == 0, (1 - pSuccess), pSuccess)))
                  }
                }

                post_logLik = post_logLik * k + logLik/nSums
            }
        }#end of saving samples and computing running means

        if (verbose) {
            cat("---------------------------------------\n")
            tmp = proc.time()[3]
            cat(c(paste(c("  Iter=", "Time/Iter=", "varE="), round(c(i, c(tmp - time), varE), 3), sep = "")), "\n")
            time = tmp
        }
    }#end of Gibbs sampler

    #Closing files
    close(fileOutVarE)
    close(fileOutMu)

    if (nLT > 0) {
        for (i in 1:nLT) {
            if (!is.null(ETA[[i]]$fileOut)) {
                close(ETA[[i]]$fileOut)
            }
            ETA[[i]]$fileOut = NULL
        }
    }
    
    #return goodies
    if (response_type == "Bernoulli" | response_type == "ordinal") {
        y = z
    }

    out = list(y = y, whichNa = whichNa, saveAt = saveAt, nIter = nIter, 
               burnIn = burnIn, thin = thin, minAbsBeta = minAbsBeta, 
               weights = weights, ncores = ncores, verbose = verbose, 
               response_type = response_type, df0 = df0, S0 = S0)

    out$yHat = post_yHat
    out$SD.yHat = sqrt(post_yHat2 - (post_yHat^2))
    out$mu = post_mu
    out$SD.mu = sqrt(post_mu2 - post_mu^2)
    out$varE = post_varE
    out$SD.varE = sqrt(post_varE2 - post_varE^2)
    
    #goodness of fit 
    out$fit = list()
    
    if(response_type=="gaussian")
    {
    	tmpE = (yStar - post_yHat)/weights
    	tmpSD = sqrt(post_varE)/weights
    
    	if (nNa > 0) {
        	tmpE = tmpE[-whichNa]
        	tmpSD = tmpSD[-whichNa]
    	}
    	out$fit$logLikAtPostMean = sum(dnorm(tmpE, sd = tmpSD, log = TRUE))

	if (Censored) {
            cdfA = pnorm(q = a[whichNa], sd = sqrt(post_varE), mean = post_yHat[whichNa])
            cdfB = pnorm(q = b[whichNa], sd = sqrt(post_varE), mean = post_yHat[whichNa])
            out$fit$logLikAtPostMean = out$fit$logLikAtPostMean + sum(log(cdfB - cdfA))
        }
    }
   
    if(response_type=="Bernoulli")
    {
	#Be careful, now the vector with response is y
	if(nNa>0)
        {
	   pSuccess = pnorm(mean = 0, sd = 1, q = post_yHat[-whichNa])
           out$fit$logLikAtPostMean = sum(log(ifelse(y[-whichNa] == 0, (1 - pSuccess), pSuccess)))
        }else
        {
           pSuccess = pnorm(mean = 0, sd = 1, q = post_yHat)
           out$fit$logLikAtPostMean = sum(log(ifelse(y == 0, (1 - pSuccess), pSuccess)))
        }
    }

    if(response_type=="ordinal")
    {
         out$fit$logLikAtPostMean = loglik_ordinal(y,post_yHat,post_threshold)
         out$threshold = post_threshold[-c(1, nclass + 1)]
         out$SD.threshold = sqrt(post_threshold2 - post_threshold^2)[-c(1, nclass + 1)] 
    }

    out$fit$postMeanLogLik = post_logLik
    out$fit$pD = -2 * (post_logLik - out$fit$logLikAtPostMean)
    out$fit$DIC = out$fit$pD - 2 * post_logLik

    # Renaming/removing objects in ETA
    if (nLT > 0) {
        for (i in 1:nLT) {

            if (ETA[[i]]$model != "RKHS") {
                ETA[[i]]$b = ETA[[i]]$post_b
                ETA[[i]]$SD.b = sqrt(ETA[[i]]$post_b2 - ETA[[i]]$post_b^2)
                tmp = which(names(ETA[[i]]) %in% c("post_b", "post_b2"))
                ETA[[i]] = ETA[[i]][-tmp]
            }
            
            if(ETA[[i]]$model=="RKHS")
            {
               ETA[[i]]$SD.u=sqrt(ETA[[i]]$post_u2 - ETA[[i]]$post_u^2)
            }

            if (ETA[[i]]$model %in% c("BRR", "BayesA", "BayesC")) {
                ETA[[i]]$varB = ETA[[i]]$post_varB
                ETA[[i]]$SD.varB = ETA[[i]]$post_varB2 - (ETA[[i]]$post_varB^2)
                tmp = which(names(ETA[[i]]) %in% c("post_varB", "post_varB2"))
                ETA[[i]] = ETA[[i]][-tmp]
            }
        }
        out$ETA = ETA
    }
    class(out) = "BGLR"
    return(out)
}


#This function will be a wrapper for BGLR
#the idea is to maintain the compatibility with the function BLR in 
#the package BLR that was released in 2010, updated in 2011 and 2012

#FIXME: thin2 parameter is missing in BGLR

BLR=function (y, XF = NULL, XR = NULL, XL = NULL, GF = list(ID = NULL, 
    A = NULL), prior = NULL, nIter = 1100, burnIn = 100, thin = 10, 
    thin2 = 1e+10, saveAt = "", minAbsBeta = 1e-09, weights = NULL, 
    ncores = 1) 
{

    ETA = NULL
    ETA = list()
    nLT = 0

    cat("This implementation is a simplified interface for the more general\n")
    cat("function BGLR, we keep it for backward compatibility with our package BLR\n")

    warning("thin2 parameter is not used any more and will be deleted in next releases\n",immediate. = TRUE);
   
    cat("Setting parameters for BGLR...\n")
    if (is.null(prior)) {
        cat("===============================================================\n")
        cat("No prior was provided, BGLR will be running with improper priors.\n")
        cat("===============================================================\n")
        prior = list(varE = list(S = NULL, df = 1), varBR = list(S = 0, 
            df = 0), varU = list(S = 0, df = 0), lambda = list(shape = 0, 
            rate = 0, type = "random", value = 50))
    }
    if (!is.null(XF)) {
        nLT = nLT + 1
        ETA[[nLT]] = list(X = XF, model = "FIXED")
    }
    if (!is.null(XR)) {
        nLT = nLT + 1
        ETA[[nLT]] = list(X = XR, model = "BRR", df0 = prior$varBR$df, 
            S0 = prior$varBR$S)
    }
    if (!is.null(XL)) {
        nLT = nLT + 1
        if (prior$lambda$type == "random") {
            if (is.null(prior$lambda$rate)) {
                cat("Setting prior for lambda^2 to beta\n")
                prior$lambda$type = "beta"
                prior$lambda$shape = NULL
                prior$lambda$rate = NULL
            }
            else {
                cat("Setting prior for lambda^2 to gamma\n")
                prior$lambda$type = "gamma"
                prior$lambda$max = NULL
                prior$lambda$shape1 = NULL
                prior$lambda$shape2 = NULL
            }
        }
        ETA[[nLT]] = list(X = XL, model = "BL", type = prior$lambda$type, 
                          rate = prior$lambda$rate, shape = prior$lambda$shape, 
                          max = prior$lambda$max, shape1 = prior$lambda$shape1, 
                          shape2 = prior$lambda$shape2, lambda = prior$lambda$value)
    }

    #FIXME: In original BLR IDS are used to buid A matrix Z,
    #and then run the model y=Zu+e, u~MN(0,varU*A), and then using the Cholesky factorization
    #it was possible to fit the model. The algorithm used here is different (Orthogonal variables)
    #and it may be that the IDS are not longer necessary

    if (!is.null(GF[[1]])) {
        nLT = nLT + 1
        ETA[[nLT]] = list(K = GF$A, model = "RKHS", df0 = prior$varU$df, 
            S0 = prior$varU$S)
        warning("IDs are not used any more and will be deleted in next releases...\n",immediate. = TRUE) 
    }

    cat("Finish setting parameters for BGLR\n")
    cat("Fitting model using BGLR...\n")
    out = BGLR(y = y, ETA = ETA, df0 = prior$varE$df, S0 = prior$varE$S, 
               nIter = nIter, burnIn = burnIn, thin = thin, saveAt = saveAt, 
               minAbsBeta = minAbsBeta, weights = weights, ncores = ncores)

    #Backward compatibility with BLR
    if (nLT > 0) {
        for (j in 1:nLT) {
            if (ETA[[j]]$model == "FIXED") {
                out$bF = out$ETA[[j]]$b
                out$SD.bF = out$ETA[[j]]$SD.b
            }
            if (ETA[[j]]$model == "BL") {
                out$bL = out$ETA[[j]]$b
                out$SD.bL = out$ETA[[j]]$SD.b
            }
            if (ETA[[j]]$model == "BRR") {
                out$bR = out$ETA[[j]]$b
                out$SD.bR = out$ETA[[j]]$SD.b
                out$varBR = out$ETA[[j]]$varB
                out$SD.bR = out$ETA[[j]]$SD.varB
            }
            if (ETA[[j]]$model == "RKHS") {
                out$u = out$ETA[[j]]$post_u
                out$SD.u =out$ETA[[j]]$SD.u
                out$varU = out$ETA[[j]]$post_varU
            }
        }
    }
    out$ETA = NULL
    class(out) = "BLR"
    return(out)
}


if(FALSE){

###BAYESIAN ELASTIC NET ###########################################################################################################################################   
#Bayesian Elastic Net LASSO
#Kyung et al., 2010
setLT.BEN=function(LT,n,j,y,weights,nLT,R2,saveAt)
{
    LT$X=as.matrix(LT$X)
    LT$p=ncol(LT$X)

    if(any(is.na(LT$X))){stop(paste("LP ",j," has NAs in X",sep=""))}
    if(nrow(LT$X)!=n){stop(paste("   Number of rows of LP ",j,"  not equal to the number of phenotypes.",sep=""))}

    for (i in 1:n) { LT$X[i, ] = weights[i]*LT$X[i, ]  }

    LT$x2=rep(0,LT$p)
    sumMeanXSq=0
    for(i in 1:LT$p)
    {
      LT$x2[i]=sum(LT$X[,i]^2)
      sumMeanXSq=sumMeanXSq+mean(LT$X[,i])^2
    }
    LT$MSx=sum(LT$x2)/n-sumMeanXSq

    # Prior
    if(is.null(LT$type))
    {
                LT$type="gamma"
                cat(paste("  By default, lambda1^2 and lambda2 in LP ",j,"  were set to gamma.\n",sep=""))
    }
   
    if(LT$type!="gamma") stop("Only gamma type priors are allowed for lambda1^2 and lambda2");

    #Lambda1
    if(is.null(LT$lambda1))
    {
        LT$lambda1=0.5
        cat(paste("  Initial value of lambda1 in LP ",j," was missing and was set to ",LT$lambda1,"\n",sep=""))
    }
    if(is.null(LT$shape1))
    {
        LT$shape1=10
        cat(paste("  shape1 parameter for lambda1^2 in LP ",j," was missing and was set to ",LT$shape1,"\n",sep=""))
    }
    if(is.null(LT$rate1))
    {
         LT$rate1=1
         cat(paste("  rate1 parameter for lambda1^2 in LP ",j," was missing and was set to ",LT$rate1,"\n",sep=""))
    } 

    #Lambda2
    if(is.null(LT$lambda2))
    {
        LT$lambda2=0.5
        cat(paste("  Initial value of lambda2 in LP ",j," was missing and was set to ",LT$lambda2,"\n",sep=""))
    }
    if(is.null(LT$shape2))
    {
        LT$shape2=10
        cat(paste("  shape2 parameter for lambda2 in LP ",j," was missing and was set to ",LT$shape2,"\n",sep=""))
    }
    if(is.null(LT$rate2))
    {
         LT$rate2=1
         cat(paste("  rate2 parameter for lambda2 in LP ",j," was missing and was set to ",LT$rate2,"\n",sep=""))
    }

    LT$b=rep(0,LT$p)
    LT$post_b=rep(0,LT$p)
    LT$post_b2=rep(0,LT$p)
    
    tmp=((var(y,na.rm=TRUE)*R2/nLT)/(LT$MSx))
    LT$tau2=rep(tmp,LT$p)
    LT$post_tau2=rep(0,LT$p)
    LT$post_lambda1=0
    LT$post_lambda2=0

    #Output_file
    fname=paste(saveAt,"ETA_",j,"_lambda1_lambda2.dat",sep="")
    LT$NamefileOut=fname
    LT$fileOut=file(description=fname,open="w")
   
    LT$X=as.vector(LT$X)
    return(LT)

}

##########################################
                    if(ETA[[j]]$model=="BEN")
                    {
                        #Update Betas
 			varBj = varE * (ETA[[j]]$tau2)/(1+ETA[[j]]$lambda2*ETA[[j]]$tau2)
                    	ans = .Call("sample_beta", n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2,  ETA[[j]]$b,
                                                e, varBj, varE, minAbsBeta,ncores)
                    	ETA[[j]]$b = ans[[1]]
                    	e=ans[[2]]
                        
                        #Update tau2
                        nu=sqrt(ETA[[j]]$lambda1^2*varE/ETA[[j]]$b^2)         
                        tmp=try(1/rinvGauss(ETA[[j]]$p,nu, ETA[[j]]$lambda1^2))
                	if(any(is.nan(tmp)) | any(is.infinite(tmp)))
                	{
                  		cat("Warning tau2 was not updated due to numeric problems with beta\n")
                	}else
                	{
                  		ETA[[j]]$tau2=tmp
                	} 

                        #Update lambda1 and lambda2
                        rate = sum(ETA[[j]]$tau2)/2 + ETA[[j]]$rate1
                        shape = ETA[[j]]$p + ETA[[j]]$shape1
                        ETA[[j]]$lambda1 = sqrt(rgamma(rate = rate, shape = shape, n = 1))

                        rate=sum(ETA[[j]]$b^2)/(2*varE) + ETA[[j]]$rate2
                        shape=ETA[[j]]$p/2+ETA[[j]]$shape2
                        ETA[[j]]$lambda2 = rgamma(1,shape=shape,rate=rate)

                        deltaSS=deltaSS+sum(ETA[[j]]$b^2*(1/ETA[[j]]$tau2+ETA[[j]]$lambda2))
                        deltadf=deltadf+ETA[[j]]$p

                    } #END BEN

##############################################################################################################################################	

}
