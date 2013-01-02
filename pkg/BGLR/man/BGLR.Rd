\name{BGLR}
\alias{BGLR}
\title{Bayesian Generalized Linear Regression}
\description{
  The BGLR (`Bayesian Generalized Linear Regression') function
  was designed to fit parametric regression models using different
  types of shrinkage methods. Several of the models implemented in this function  
  were presented in de los Campos \emph{et al.} (2009, 2010).
}
\usage{
  BGLR(y, response_type = "gaussian", a=NULL, b=NULL,ETA = NULL, nIter = 1500,
       burnIn = 500, thin = 5, saveAt = "", S0 = NULL, 
       df0 =5, R2 = 0.5, minAbsBeta = 1e-09, weights = NULL,
       ncores = 1, verbose = TRUE, rmExistingFiles = TRUE)

}
\arguments{
   \item{y}{(numeric, \eqn{n}{n}) the data-vector (NAs allowed).}
   \item{response_type}{string, specify the distribution of the response variable, right now only gaussian, Bernoulli and ordinal responses are allowed.}
   \item{a}{vector for specifying lower bound for censored observations, default value NULL. See details.}
   \item{b}{vector for specifying upper bound for censored observations, default value NULL. See details.}
   \item{ETA}{A list of predictors and prior specifications for the regression coefficients. For example the prior for the 
             regression coefficients can be that used in Bayesian LASSO, Bayesian ridge regression, BayesA, BayesB, BayesC-pi, Elastic Net LASSO, etc. See details below.}
   \item{weights}{(numeric, \eqn{n}{n}) a vector of weights, may be NULL.}
   \item{nIter,burnIn, thin}{(integer) the number of iterations, burn-in and thinning.}
   \item{saveAt}{(string) this may include a path and a pre-fix that will be added to the name of the files that are saved as the program runs.}
   \item{ncores}{Number of CPUs used in computations, default value 1.}
   \item{S0}{The scale parameter for the scaled inverse-chi squared distribution for \eqn{\sigma_e^2}{sigma[e]^2}.}
   \item{df0}{The degrees of freedom for the scaled inverse-chi squared distribution for \eqn{\sigma_e^2}{sigma[e]^2}.}
   \item{R2}{...}
   \item{minAbsBeta}{The minimum absolute value of the components of \eqn{\boldsymbol \beta_L}{bL} to avoid numeric problems when sampling from \eqn{\boldsymbol \tau^2}{\tau^2}, default \eqn{1 \times 10^{-9}}{1e-9}.}
   \item{verbose}{logical, if TRUE prints iteration history, defalt TRUE.}
   \item{rmExistingFiles}{logical, if TRUE removes existing output files from previous runs, default value is TRUE.}
}

\details{

The program run a Gibbs sampler for the regression model given below.

\bold{Likelihood}. The equation for the data is:

%                             __ Hbeta                   __ Hu          
%y =   1 mu  +  X   beta   +  \          X    beta    +  \          u   + e 
%                F      F     /__ h = 1   Rh      Rh     /__ h = 1   h 

\deqn{
 \begin{array}{lr}
 \boldsymbol y= \boldsymbol 1 \mu + \boldsymbol X_F \boldsymbol \beta_F + \sum_{h=1}^{H\beta} \boldsymbol X_{Rh} \boldsymbol \beta_{Rh} + \sum_{h=1}^{Hu} \boldsymbol u_h + \boldsymbol \varepsilon & (1)
 \end{array}
}{y=1*mu + XF*bF + sum(XR[h]*bR[h],h) + sum(u[h],h) + e ...(1)}


where \eqn{\mu}{mu} is an effect common to all individuals, \eqn{\boldsymbol X_F=\{ x_{Fij} \}}{XF} 
represent covariates whose effects \eqn{\boldsymbol \beta_{F}=\{\beta_{Fj}\}}{bF} will be estimated 
shrinkage (the so-called `fixed effecs', e.g., age, sex), \eqn{\boldsymbol X_{Rh}=\{x_{Rhij}\}}{XR[h]} 
represent covariates whose effects \eqn{\boldsymbol \beta_{Rh}=\{\boldsymbol \beta_{Rhj}\}}{bR[h]} will be treated as `random effects'
and will be estimated using shrinkage estimation methods (non-flat priors in a Bayesian context) and 
\eqn{\boldsymbol u_h=\{u_{hi}\}}{u[h]} are random effects used to describe, for example, a regression on a 
pedigree or a RKHS regression on markers.

\bold{Prior}

The model specification is complete once we assign a prior distribution to the model unknowns. 
The intercept \eqn{\mu}{mu} and \eqn{\boldsymbol \beta_F}{bF} are assigned flat priors, while
\eqn{\boldsymbol \beta_{Rh}}{bR[h]}, \eqn{\boldsymbol u_h}{u[h]} and \eqn{\sigma_e^2}{sigma[e]^2} are assigned non 
flat priors, denoted as \eqn{p(\boldsymbol \beta_R)}{p(bR)}, \eqn{p(\boldsymbol u)}{p(u)} and \eqn{p(\sigma_e^2)}{p(sigma[e]^2)}, 
respectively. The structure of the priors is as follows:

\deqn{
  \begin{array}{lr}
   p(\mu,\boldsymbol \beta_F, \boldsymbol \beta_{R1},...,\boldsymbol \beta_{RH\beta},\boldsymbol u_1,...,\boldsymbol u_{Hu},\sigma_e^2) \propto % 
   \left\{\prod_{h=1}^{H\beta} p (\boldsymbol \beta_{Rh})\right\} \left\{\prod_{h=1}^{Hu} p(\boldsymbol u_h) \right\} \chi^{-2}(\sigma_e^2| df,S), & (2)
\end{array}
}{p(mu,bF,bR[1],...,bR[Hb],u[1],...,u[Hu],sigma[e]^2) propto prod(p(bR[h]),h) prod(u[h],h) Inv_Scaled_chisq(sigma[e]^2 | df,S) ...(2)}

where \eqn{\chi^{-2} (\sigma^2 | df,S)}{Inv_Scaled_chisq(sigma^2 | df,S)} is a scaled-inverse Chi-square 
density assigned to \eqn{\sigma^2}{sigma^2} with degree of freedom and 
scale parameter \eqn{df}{df} and \eqn{S}{S} respectively.

The prior distribution assigned to \eqn{p(\boldsymbol u_h | \boldsymbol \theta_{uh})}{p(u[h]|theta[uh])} is multivariate normal centered at
zero and with covariance \eqn{ \sigma_{uh}^2 \boldsymbol K_{uh}}{sigma[uh]^2 K[uh]} where \eqn{\boldsymbol K_{uh}}{K[uh]} is a positive definite-matrix
 and \eqn{\sigma_{uh}^2}{sigma[uh]^2} is an unknow variance paramter. The prior assigned to this parameter is 
a scaled inverse chi-squared so that 

\deqn{
  \begin{array}{lr}
    p(\boldsymbol u_h, \sigma_{uh}^2 ) = N(\boldsymbol u_h| \boldsymbol 0, \sigma_{uh}^2 \boldsymbol K_{uh}) \chi^{-2} (S_{uh}, df_{uh})& (3)
  \end{array}
}{p(u[h],sigma[uh]^2)=N(u[h]|0,sigma[uh]^2 K[uh]) Inv_Scaled_chisq(df[uh],S[uh])}

Following standard assumptions of Bayesian regression models, regression coefficients
are assigned IID priors; therefore: 
\eqn{p(\boldsymbol \beta_{Rh} | \boldsymbol \theta_{Rh})=\left\{ \prod_{j=1}^{p_{Rh}}  p(\beta_{Rhj} | \boldsymbol \theta_{Rh}) \right\} p(\boldsymbol \theta_{Rh})}{p(bR[h]|thetaR[h])=prod(bR[hj]|thetaR[h],j) p(thetaR[h])},
where \eqn{p(\beta_{Rhj} | \boldsymbol \theta_{Rh})}{p(bR[hj]|thetaR[h])}  can be a double exponential distribution, 
a normal distribution, etc., \eqn{\boldsymbol \theta_{Rh}}{thetaR[h]} is a vector of
unknown indexing the prior density assigned to marker effects and \eqn{p(\boldsymbol \theta_{Rh})}{p(thetaR[h])} is the prior assigned
to these unknowns. 

Collecting assumptions we have:

\deqn{
  \begin{array}{ll}
    p(\mu, \boldsymbol \beta_F,\boldsymbol \beta_{R1},....,\boldsymbol \beta_{RH\beta},\boldsymbol u_1,...,\boldsymbol u_{Hu}, \sigma^2) %% 
             & \propto \prod_{h=1}^{H\beta} \left\{\prod_{j=1}^{p_{Rh}}  p(\beta_{Rhj} | \boldsymbol \theta_{Rh}) \right\} p(\boldsymbol \theta_{Rh}) \\
             & \times \left\{ \prod_{h=1}^{Hu} N(\boldsymbol u_h | \boldsymbol 0, \sigma_{uh}^2 \boldsymbol K_{uh}) \chi^{-2} (\sigma_{uh}^2 | S_{uh},df_{uh}) \right\} \\
             & \times \chi^{-2} (\sigma^2 | df, S)
  \end{array}
}{p(mu,bF,bR[1],...,bR[Hb],u[1],...,u[Hu],sigma^2) propto prod(prod(p(bR[hj]|thetaR[h]),j)p(thetaR[h]),h) * prod(N(u[h]|0,sigma[uh]^2*K[uh]) * Inv_Scaled_chisq(sigma[uh]^2 | S[uh],df[uh]),h) * Inv_Scaled_chisq(sigma^2| df, S)} 

\bold{Special cases}

\emph{Bayesian Gaussian Regression (BGR)}

A common approach in Bayesian shrinkage estimation is to assign independent and identically distributed (IDD) conditional Gaussian priors with 
unknown variance. This can be impemented by setting 

\deqn{
       p(\boldsymbol \beta_{Rh} | \boldsymbol \theta_{Rh}) = \left\{ \prod_{j=1}^{p_{Rh}} N(\beta_{Rhj} | 0, \sigma_{\beta h}^2) \right\} \chi^{-2} (\sigma_{\beta h}^2 | df_{\beta h}, S_{\beta h}).
     }{p(beta[Rh]| theta[Rh])=prod(N(beta[Rhj]|0,sigma[beta h]^2),j) * Inv_Scaled_chisq(sigma[beta h]^2 | df[beta h], S[beta h])}

When \eqn{\sigma_{\beta h}^2}{sigma[beta h]^2} is known, using this prior yield estimates which are equivalent to those of a RR. In a BGR the excent of 
shrinkage is controlled by the variance (or noise-to-signal) ratio \eqn{\lambda_h=\sigma_e^2/\sigma_{\beta h}^2}{lambda[h]=sigma[e]^2/sigma[beta h]^2}. 
This quantity is the same for all regression coefficients included in \eqn{\boldsymbol \beta_{Rh}}{beta[Rh]}; this may not be appropiate if some markers are located
in regions harboring QTL while others are located in regions which are not associated to genetic variance. To overcome this problem, 
alternative shrinkage procedures such as those described below can be used.

\emph{Mixtures of scaled-normal densities}

This class of mixtures can be used as prior of marker effects to obtain a type of shrinkage different than that
of a BGR. Examples of this are the double-exponential (DE) and scaled-t densities, which are commonly used as prior of
maker effects in Whole Genomic Prediction (WPG). The results models are known as the Bayesian LASSO (BL) and BayesA respectively.
Relative to the Gaussian density used in BGR, the DE and the scaled-t densities have higer mass at zero and thicker tails, inducing a 
different type of shrinkage. The DE and scaled-t prior densities can be represented as mixtures of scaled normal-densities
of the form

\deqn{
       p(\beta_{Rhj} | H) = \int N(\beta_{Rhj} | 0, \sigma_{\beta hj}^2) p(\sigma_{\beta hj}^2 | H)  \partial \sigma_{\beta hj}^2
    }{p(beta[Rhj]| H)= Integrate(N(beta[Rhj]|0,sigma[beta hj]^2)) * p(sigma[beta hj]^2 | H) d sigma[beta hj]^2}

where \eqn{\sigma_{\beta hj}^2}{sigma[beta hj]^2} is a marker-specific variance parameter, 
\eqn{p(\sigma_{\beta hj}^2 | H)}{p(sigma[beta hj]^2|H)} is a prior density assigned to this variance parameter
and \eqn{H}{H} is a set of hyperparameters which may be specified a-priori or estimated from the 
data. When \eqn{p(\sigma_{\beta hj}^2 | H)}{p(sigma[beta hj]^2 | H)} is an exponential (scaled-inverse chi-square) density, 
the resulting marginal prior density of marker effects is a double-exponential (scaled-t). 

\emph{Pedigree-based regressions}

They represent a generalization of the concept of `family history' to complex genealogies.
These regressions have been used over more than 5 decades for prediction of genetic values in animal and
plant breeding applications. Pedigree regressions can be implemented by setting \eqn{\boldsymbol K_{hu}=\boldsymbol A}{K[hu]=A} where 
\eqn{\boldsymbol A=\left\{a(i,i')\right\}}{A=a(i,i')} is a matrix whose entries are twice the coefficient of kinship between individuals, 
which can be computed from a pedigree.

\emph{Reproducing Kernel Hilbert Regressions(RKHS)}

RKHS are used for semi-parametric regressions in applications as diverse as scatter-plot smoothing
(smooting spline), spacial statistics (Kriging), gene expression or WGP. Estimates from RKHS can 
be motivated as the solution to a penalized optimization problem of as posterior modes in certain class
of Bayesian models. A Bayesian formulation of RKHS can be implemented  by simply setting 
\eqn{\boldsymbol K_{uh}=\left\{ K_{uh} (i,i') \right\}}{K[uh]=K[uh](i,i')} to be a matrix whose entries contain 
the evaluations of a reproducing kernel at pairs of poings \eqn{(i,i')}{(i,i')}. In WGP
models the reproducing kernel, \eqn{K(i,i')=K(\boldsymbol z_i, \boldsymbol z_i')}{K(i,i')=K(z[i],z[i'])}, maps
from pairs of marker genotypes \eqn{(\boldsymbol z_i, \boldsymbol z_i')}{z[i],z[i']} onto  co-variance function.
For instance, using the Gaussian kernel, \eqn{K(i,i')=\exp\left( -\omega ||\boldsymbol z_i - \boldsymbol z_{i'}||^2\right)}{K(i,i')=exp(-w||z[i]-z[i']||)},
where \eqn{||\boldsymbol z_i - \boldsymbol z_{i'}||}{z[i]-z[i']}  is a Euclidean distance between the two 
vectors of marker genotypes and \eqn{\omega}{w} is a bandwith parameter.

\emph{Censored outcomes}

In BGLR censored outcomes are dealt with as a missing data problem. BGLR handles three types of 
censoring: left, right and interval censored. For an interval censored data-point the information 
available is \eqn{a_i< y_i < b_i}  where: \eqn{a_i}  and \eqn{b_i}  are known lower and upper bounds 
and \eqn{y_i} is the actual phenotype which for censored data points is un-observed. Rigth censoring 
occurs when \eqn{b_i} is also uknown, therefore, the only information available is  \eqn{a_i<y_i}. In 
a time-to-event setting this means that we know that time to event exceeded the time at 
censoring given by \eqn{a_i}. Left censoring occurs when \eqn{b_i} is unknown; therefore, 
the only information available is: \eqn{y_i<b_i}. In BGLR censored outcomes are then specified 
with three vectors, \eqn{\boldsymbol y}, \eqn{\boldsymbol a} and \eqn{\boldsymbol b}. The configuration of the triplet   for un-censored, right-censored, 
left-censored and interval censored are described in the table below. 

\tabular{lccc}{
                     \tab a             \tab y          \tab b  \cr
   Un-censored       \tab NULL          \tab \eqn{y_i}  \tab NULL \cr
   Right censored    \tab \eqn{a_i}     \tab NA         \tab \eqn{\infty} \cr
   Left censored     \tab \eqn{-\infty} \tab NA         \tab \eqn{b_i} \cr
   Interval censored \tab \eqn{a_i}     \tab NA         \tab \eqn{b_i}
 } 

The only modification introduced in the Gibbs sampler required for handling censored data points 
consist of sampling, at each iteration of the Gibbs sampler, the censored phenotypes form 
the corresponding fully-conditional densities which in BGLR are truncated normal densities.

\emph{Binary outcomes}

They can be modeled using the threshold model, or probit link. Here, probability of success is  \eqn{P(Y_i=1)=\Phi(\eta_i)}  
where \eqn{\Phi(\cdot)} is the standard normal cumulative distribution function 
(also known as normal probit link) and  \eqn{\eta_i} is a linear predictor which 
can include the type of fixed or random effects handled by BGLR. In order to run a regression for 
binary outcomes, the response must be coded with 0's (failure) and 1's (success), and 
the argument response_type should be set to "Bernoulli". More details about this model
can be found in Albert & Chib (1993).


\emph{Ordinal outcomes}

They can be modeled using also the threshold model. Here we model, 
\eqn{\pi_{ij}=P(Y_i \leq j)=\Phi(\eta_{ij})}, where \eqn{\eta_{ij}=\gamma_j- \boldsymbol x_i' \boldsymbol \beta},
where \eqn{\Phi(\cdot)} is the standard normal cumulative distribution function, \eqn{\gamma_j} is a threshold, the thresholds must satisfy, 
\eqn{-\infty = \gamma_0 < \gamma_1 < \cdots < \gamma_J=\infty}, \eqn{J} is the cardinality of 
\eqn{\boldsymbol{y}}{y}. In order to run a regression for ordinal outomes, the response must be
coded as \eqn{1,...,J}, and the data should be ordered accordingly, the argument response_type should 
be set to "ordinal". More details about this model can be found in Albert & Chib (1993).

}

\value{

A list with posterior means, posterior standard deviations, and the parameters used to fit the model:

}

\references{

Albert J,. S. Chib. 1993. Bayesian Analysis of Binary and Polychotomus Response Data. \emph{JASA}, \bold{88}: 669-679.

de los Campos G., H. Naya, D. Gianola, J. Crossa, A. Legarra, E. Manfredi, K. Weigel and J. Cotes. 2009.
Predicting Quantitative Traits with Regression Models for Dense Molecular Markers and Pedigree. \emph{Genetics} \bold{182}: 375-385.

de los Campos, G., D. Gianola, G. J. M., Rosa, K. A., Weigel, and J. Crossa. 2010.  Semi-parametric genomic-enabled prediction of genetic values using 
reproducing kernel Hilbert spaces methods. \emph{Genetics Research}, \bold{92}:295-308.

Park T. and G. Casella. 2008. The Bayesian LASSO. \emph{Journal of the American Statistical Association} \bold{103}: 681-686.

Spiegelhalter, D.J., N.G. Best, B.P. Carlin and A. van der Linde. 2002. Bayesian measures of model complexity and 
fit (with discussion). \emph{Journal of the Royal Statistical Society}, Series B (Statistical Methodology) \bold{64} (4): 583-639.
}

\author{
Gustavo de los Campos, Paulino Perez Rodriguez,
}
\examples{

\dontrun{
#Demos
library(BGLR)

#BayesA
demo(BA)

#BayesB
demo(BB)

#Bayesian LASSO
demo(BL)

#Bayesian Ridge Regression
demo(BRR)

#BayesCpi
demo(BayesCpi)

#RKHS
demo(RKHS)

#Binary traits
demo(Bernoulli)

#Ordinal traits
demo(ordinal)

#Censored traits
demo(censored)

}

}
\keyword{models}
