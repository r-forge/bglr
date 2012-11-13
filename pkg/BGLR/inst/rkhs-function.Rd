\name{RKHS}
\alias{RKHS}
\title{RKHS}
\description{
  Implements Reproducing Kernel Hilbert Spaces.
}

\usage{
  RKHS(y, a = NULL, b = NULL, family = "gaussian", K, XF = NULL,
       df0 = 0, S0 = 0, nIter = 110, burnIn = 10, thin = 10, 
       saveAt = "", weights = NULL)
}

\arguments{
  \item{y}{(numeric, \eqn{n}{n}) the data-vector (NAs allowed).}
  \item{a,b}{upper and lower limit for non censored observations.}
  \item{family}{character string, the parametric family of the outcome, it can be "gaussian" or "bernoulli". The default value is "gaussian".}
  \item{K}{a list of lists, each component of the list will contain the kernel K, the prior degrees of freedom df0 and scale parameter S0 for 
          the variance componente associated to each random effect.}
  \item{XF}{(numeric, \eqn{n \times pF}{n x pF}) incidence matrix for \eqn{\boldsymbol b_F}{bF}, the default value is NULL.}
  \item{df0}{numeric, degrees of freedom for inverted chi squared distribution used as prior for variance component associted with the 
             error, i.e. \eqn{\sigma_{\boldsymbol \varepsilon}^2}{varE}, the default value is 0.}
  \item{S0}{numeric, scale parameter for inverted chi squared distribution used as prior for the variance component associated with the error, the default value is 0.}
  \item{nIter,burnIn, thin}{(integer) the number of iterations, burn-in and thinning.}
  \item{saveAt}{character string, this may include a path and a pre-fix that will be added to the name of the files that are saved as the program runs.}
  \item{weights}{(numeric, \eqn{n}{n}) a vector of weights, may be NULL.}
}

\details{
  Implements a linear regression model of the form:

  \deqn{\boldsymbol y=\boldsymbol{\mu}+\boldsymbol X_F \boldsymbol b_F + \boldsymbol u_1 + \boldsymbol u_2 +\cdots+\boldsymbol u_q +\boldsymbol{e}}{y=1mu+XF*bF + u1 + u2 +...+uq +e}

  where \eqn{\boldsymbol y}{y} is a vector of phenotypes, \eqn{\mu} is an intercept, \eqn{\boldsymbol X_F} is a matrix with fixed effects,
  \eqn{\boldsymbol u_1,...,\boldsymbol u_q} are random effects with the following prior \eqn{\boldsymbol u_j \sim  N(\boldsymbol 0, \sigma^2_j \boldsymbol K_j)}
  here \eqn{\boldsymbol K_j} is a positive definite matrix.

  The problem is solved by using a Gibbs sampler. This model allows to implement predigree regression, GBLUP and Reproducing Kernel Hilbert Spaces.

}

\value{
   A list with posterior means, posterior standard deviations, and the parameters used to fit the model:
   
   \item{$mu}{the posterior mean of the intercept.}
   \item{$varE}{the posterior mean of \eqn{\sigma_{\boldsymbol \varepsilon}^2}{varE}.}
   \item{$yHat}{the posterior mean of \eqn{\boldsymbol y=\boldsymbol{\mu}+\boldsymbol X_F \boldsymbol b_F + \boldsymbol u_1 + \boldsymbol u_2 +\cdots+\boldsymbol u_q +\boldsymbol{e}}{y=1mu+XF*bF + u1 + u2 +...+uq +e}}
   \item{$weights}{vector of weights.}
   \item{$whichNa}{a vector indicating which entries in \eqn{\boldsymbol y}{y} were missing.}
   \item{$fit}{a list containing the following elements:
               \itemize{
                  \item $postMeanLogLik: the posterior mean of the Log-Likelihood.
                  \item $pD: estimated effective number of parameters, Spiegelhalter \emph{et al.} (2002).
                  \item $DIC: the deviance information criterion, Spiegelhalter \emph{et al.} (2002). 
                  \item $bF: the posterior mean of \eqn{\boldsymbol b_F}{bF}.
                  \item $SD.bF: the corresponding posterior standard deviation.
               }
             }
   \item{$K}{a list of lists. Each element in the list will contain the following elements:
              \itemize{
                \item $u 
                \item $uStar
                \item $varU
                \item $cumMSa
                \item $df0
                \item $S0
                \item $tolD 
              }
            }
   \item{$nIter}{the number of iterations made in the Gibbs sampler.}
   \item{$burnIn}{the nuber of iteratios used as burn-in.}
   \item{$thin}{the thin used.}
   \item{$df0}{degrees of freedom for the inverse chi-squared distribution used as prior for the variance component associated with the error.} 
   \item{$S0}{scale parameter for inverted chi-squared distribution used as prior for variance componente associated with the the error.}

The posterior means returned by RKHS are calculated after burnIn is
passed and at a thin as specified by the user.

\bold{Save}. The routine will save samples of \eqn{\mu}{mu}, \eqn{\boldsymbol b_F}{b_F}, variance components.

}

\references{
de los Campos, G., D. Gianola, and G.J.M. Rosa. 2009b. Reproducing kernel Hilbert
spaces regression: a general framework for genetic evaluation. \emph{J. Anim. Sci.} 
87:1883-1887, doi:10.2527/jas.2008-1259.

de los Campos, G., D. Gianola, and G.J.M. Rosa, K.A. Weigel, and J. Crossa.
2010. Semi-parametric genomic-enabled prediction of genetic values using
reproducing kernel Hilbert spaces methods. \emph{Genet. Res.} \bold{92:} 295-308.

Spiegelhalter, D.J., N.G. Best, B.P. Carlin and A. van der Linde. 2002. Bayesian measures of model complexity and 
fit (with discussion). \emph{Journal of the Royal Statistical Society}, Series B (Statistical Methodology) \bold{64} (4): 583-639.
}

\author{
  Gustavo de los Campos, John Hickey, Paulino Perez
}
