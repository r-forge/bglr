\name{fitMixed}
\alias{fitMixed}
\title{fitMixed}
\description{
  This function fits a mixed model. See details section. 
}

\usage{
  fitMixed(y, X=NULL, Z=NULL, A=NULL,
           d=NULL,U=NULL, BLUE=TRUE, 
           BLUP=TRUE,method="ML", 
           lambda_ini=NULL, 
           tol=1e-4,maxiter=100)
}

\arguments{
  \item{y}{(numeric, \eqn{n}{n}) the data-vector (NAs not allowed).}
  \item{X}{numeric incidence matrix for fixed effects, dimmension \eqn{n \times c}{n x c} (NAs not allowed).}
  \item{Z}{numeric incidence matrix for random effects, dimmension \eqn{n \times q}{n x q} (NAs not allowed).}
  \item{A}{numeric matrix for modelling covariance for y, for example in the animal model it can 
           be a numeric relationship matrix derived from a pedigree. Dimmension \eqn{n \times n}{n x n}.}
  \item{d}{eigenvalues for \eqn{\boldsymbol G}{G} matrix. \eqn{\boldsymbol G = \boldsymbol{ZAZ}'=\boldsymbol{UDU}'}{G = ZAZ'=UDU'}.}
  \item{U}{eigenvector for \eqn{\boldsymbol G}{G} matrix. \eqn{\boldsymbol G = \boldsymbol{ZAZ}'=\boldsymbol{UDU}'}{G = ZAZ'=UDU'}.}
  \item{BLUE}{logical, if TRUE the BLUPs are computed.}
  \item{BLUP}{logical, if TRUE the BLUPs are computed.}
  \item{method}{estimation method for the variance components, can be maximum likelihood or restricted maximum likelihood. The detault is "ML" and is the 
                only one implemented now.}
  \item{lambda_ini}{Initial value for \eqn{\sigma_u^2/\sigma_e^2}{varU/varE}.}
  \item{tol}{tiny real number for declaring convergence.}
  \item{maxiter}{integer number for specifying the number of iterations.}
}

\value{
   A list with the following components:

   \item{varE}{\eqn{\hat \sigma_e^2}{varE}.}
   \item{varU}{\eqn{\hat \sigma_u^2}{varU}.}
   \item{Beta}{\eqn{\hat \beta}{Beta}.}
   \item{u}{\eqn{\hat u}{u}.} 
   \item{message}{string that informs if the method converged or not.} 
   \item{method}{method used for computing variance components.}
}

\details{

This function fits a mixed model

\deqn{\boldsymbol y= \boldsymbol X \boldsymbol \beta +  \boldsymbol{Z u} +  \boldsymbol e }{y=X*Beta + Z*u + e}

where \eqn{\boldsymbol y}{y} is the vector of phenotypes, \eqn{\boldsymbol X}{X} is the matrix for FIXED EFFECTS,
\eqn{\boldsymbol Z}{Z} matrix connection phenotypes and genotypes, \eqn{\boldsymbol u \sim N(\boldsymbol 0, \sigma_u^2)}{u ~ N(0, varU*A)}, 
\eqn{\boldsymbol A}{A} a genomic relationship matrix or pedigree and \eqn{\boldsymbol e \sim N(\boldsymbol 0, \sigma_e^2)}{e~ N(0, varE*I)}.

The model is fitted using the algorithm described in Zhou and Stephens, 2012.

}

\references{
  Zhou, X. and Stephens, Matthew. 2012. Genome-wide efficient mixed-model analysis for association studies. \emph{Nature Genetics}. \bold{47}(7) 821-824.
}

\author{Paulino Perez, Gustavo de los Campos}

\examples{

\dontrun{

library(BGLR)
library(BGLRutils)
setwd(tempdir())

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

}
