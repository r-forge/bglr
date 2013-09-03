\name{selected_eigenvalues}
\alias{selected_eigenvalues}
\concept{eigen vector}
\concept{eigen value}
\title{Spectral Decomposition of a Matrix}
\description{
  Computes eigenvalues and eigenvectors of numeric  matrices.
}
\usage{
selected_eigenvalues(x,min_eigenvalue = 1e-10)
}
\arguments{
  \item{x}{a numeric or complex matrix whose spectral decomposition is to be computed.  Logical matrices are coerced to numeric.}
  \item{min_eigenvalue}{smallest eigen_value to compute.}
}

\value{
  The spectral decomposition of \code{x} is returned as components of a
  list with components

  \item{values}{a vector containing the \eqn{p} eigenvalues of \code{x},
    sorted in \emph{decreasing} order, according to \code{Mod(values)}
    in the asymmetric case when they might be complex (even for real
    matrices).  For real asymmetric matrices the vector will be
    complex only if complex conjugate pairs of eigenvalues are detected.
  }
  \item{vectors}{a \eqn{p\times p}{p * p} matrix whose columns
    contain the eigenvectors of \code{x}.

    Recall that the eigenvectors are only defined up to a constant: even
    when the length is specified they are still only defined up to a
    scalar of modulus one (the sign for real matrices).
  }
}

\references{
  Anderson. E. and ten others (1999)
  \emph{LAPACK Users' Guide}.  Third Edition.  SIAM.\cr
  Available on-line at
  \url{http://www.netlib.org/lapack/lug/lapack_lug.html}.

  Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988)
  \emph{The New S Language}.
  Wadsworth & Brooks/Cole.
  Springer-Verlag Lecture Notes in Computer Science \bold{6}.

  Wilkinson, J. H. (1965) \emph{The Algebraic Eigenvalue Problem.}
  Clarendon Press, Oxford.
}

\seealso{
  \code{\link{svd}}, a generalization of \code{eigen}; \code{\link{qr}}, and
  \code{\link{chol}} for related decompositions.

  To compute the determinant of a matrix, the \code{\link{qr}}
  decomposition is much more efficient: \code{\link{det}}.
}

\examples{

\dontrun{

  library(BGLR)
  library(BGLRutils)
  data(wheat)
  out=selected_eigenvalues(wheat.A,min_eigenvalue=1e-5)

}

}
\keyword{algebra}
\keyword{array}
