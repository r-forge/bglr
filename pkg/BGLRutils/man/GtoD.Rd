\name{GtoD}
\alias{GtoD}
\title{GtoD}
\description{
  Computes an Euclidean distance matrix using a genomic relationship matrix.
}

\usage{
  GtoD(G)
}

\details{

\deqn{D_{ij}=G_{ii}+G_{jj}-2G_{ij}}{D[i,j]=G[i,i]+G[j,j]-2G[i,j]}

}

\arguments{
  \item{G}{Genomic relationship matrix, dimmension \eqn{n \times n}{nxn}.}
}

\value{

  The matrix with Euclidean distances.

}

\author{Paulino Perez, Gustavo de los Campos}

\examples{

\dontrun{
  library(BGLRutils)
  genos_file = system.file("extdata/mouse/bin/genosInt.bin", package="BGLRutils")
  G=G_matrix(n=1814,p=10346,file=genos_file,n_submatrix=2,
             center_internally=0,standard_internally=0)

  D=GtoD(G)
}

}
