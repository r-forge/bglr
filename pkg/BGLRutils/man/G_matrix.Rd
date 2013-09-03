\name{G_matrix}
\alias{G_matrix}
\title{G_matrix}
\description{
  Computes genomic relationship matrices allowing the user to specify several options, for 
  example centering the markers, weighting and scaling.
  
  The computations can be done in parallel. The number of cpu cores can be specified in the UNIX 
  commmand line with if you are using OpenBLAS or GOTOBlas2,

  $export OPENBLAS_NUM_THREADS=4

  or 

  $export GOTO_NUM_THREADS=4

  or

  export OMP_NUM_THREADS=4
}

\usage{
  G_matrix(n, p, n_submatrix, file, centers = NULL, weights = NULL, 
           center_internally = 1, standard_internally = 1)
}

\arguments{
  \item{n}{Number of rows.}
  \item{p}{Number of columns.}
  \item{n_submatrix}{Number of submatrices.}
  \item{file}{Input file coded as binary.}
  \item{centers}{vector for centering the columns of the input matrix.}
  \item{weights}{vector for weighting the columns of the input matrix.}
  \item{center_internally}{logical, if 1 then the routine will center the columns of the input matrix.}
  \item{standard_internally}{logical, if 1 then the routine will standardize each of the columns of the input matrix internally.}
}

\details{

This routine will calculate genomic relationship matrix (\eqn{\boldsymbol G}{G}) 
for a given matrix with markers information (\eqn{\boldsymbol X}{X}).

The Genomic relationship matrix is computed by dividing the matrix \eqn{\boldsymbol X}{X} 
into \eqn{k}{k} blocks (sub matrices):

\deqn{\boldsymbol G=\boldsymbol{XX}'=\boldsymbol X_1 \boldsymbol X_1' + \dots + \boldsymbol X_k \boldsymbol X_k'}{G = XX' = X_1 * X_1' + ... + X_k * X_k'}


The columns of \eqn{X}{X} can be centered by providing a vector with centers or 
the columns can be centered internally by subtracting the sample mean internally. The 
columns of \eqn{X}{X} can be weighted by a vector with weights and optionally standardized.

The element \eqn{i}{i}, column \eqn{j} of \eqn{X}{X} is centered and weighted as 
follows:

\deqn{Z_{ij}=(X_{ij}-c_j) \times w_j / d_j}{(X[i,j]-centers[j])*weights[j]/d[j]}

then

\deqn{\boldsymbol G=\boldsymbol{ZZ}'=\boldsymbol Z_1 \boldsymbol Z_1' + \dots + \boldsymbol Z_k \boldsymbol Z_k'}{Z = ZZ' = Z_1 * Z_1' + ... + Z_k * Z_k'}

where \eqn{c_j=0}{c[j]=0} if the vector with centers is not given or center_internally=0. 
The weights (\eqn{w_j=1}{w[j]=1}) if the vector with weights is not given.
\eqn{d_j}{d[j]} are the standard deviation for column \eqn{j}{j} and are computed internally if 
standard_internally is set to 1. If standard_internally is not set then \eqn{d_j=1}{d[j]=1}.

}


\author{Paulino Perez, Gustavo de los Campos}

\examples{

\dontrun{
library(BGLRutils)
genos_file = system.file("extdata/mouse/bin/genosInt.bin", package="BGLRutils")
A=G_matrix(n=1814,p=10346,file=genos_file,n_submatrix=2)
A=A/10346
}
}
