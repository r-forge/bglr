\name{rkhs}
\alias{rkhs}
\title{rkhs}
\description{ 
  Implements Reproducing Kernel Hilbert Spaces.
  The behaviour of the module is controlled by the par file.
  This module uses the RKHS function.
}

\arguments{

  \item{Email}{character string, E-mail of the analyst, the system will send an e-mail once the analyses finish.}
  \item{genFile}{character string, path and file name for the genotype file in gf2 format.}
  \item{phenFile}{character string, path and file name for the phenotype file in ph1 format.}
  \item{traitName}{string vector, columns in the phenFile that gives the phenotype (or adjusted phenotype), to be used for training.}
  \item{covNames}{string vector, columns in the phenFile that contain covariates  to be used as fixed effects.}
  \item{facNames}{string vector, columns in the phenFile that contain factors to be used as fixed effects.}
  \item{naiveImpute}{logical, if TRUE missing genoytypes are imputed with the mean.}
  \item{ctrMrk}{logical, if TRUE genotypes are centered.}
  \item{stdMrk}{logical, if TRUE genotypes are standarized to a unit variance.}
  \item{lowestAlleleFreq}{those marker whose minor Allele frequency is less than this value will be removed from the analysis. Default value 1/200.}
  \item{maxFreqNA}{those markers whose frecuency is biger than this value will be removed from the amalysis, default value 0.2.}
  \item{removeGenotypesNotConnectedToPhenotypes}{logical, if TRUE genotypes not Connected to Phenotypes will be removed. The default value is TRUE.}
  \item{nIter}{Number of Iterations for the Gibbs sampler.}
  \item{burIn}{Number of samples used as burn-in.}

  There are three options to provide the kernels, 

  \itemize{
    \item Option 1: The kernels are calculated by using a matrix of molecular markers, in this case
          the user should provide a file with the genotypic information matrix in pf2 format. The user also has to select
          the number of kernels and the kernel type ('gaussian' or 'linear').
    \item Option 2: Provide the kernels in lt files.
    \item Option 3: Provide eigenvalue descomposition of the kernels.
  }

}

\author{
  Gustavo de los Campos, John Hickey, Paulino Perez
}

\examples{

EMail<-'gdeloscampos@gmail.com'

# PROVIDE EITHER genFile OR KFiles or EVDFiles #############

removeGenotypesNotConnectedToPhenotypes<-FALSE

# option 1: genotype file.
genFile<-'/gensys/datasets/wheat599/gentoypes.gf2'

# if TRUE, missing genoytypes are imputed with the mean
naiveImpute<-TRUE                 
ctrMrk<-TRUE
stdMrk<-TRUE
lowestAlleleFreq<-1/200
maxFreqNA<-.2

nKernels<-3

# linear=> K=XX'  gaussian=>K=exp(-h*D2)
kernelType<-'linear' 

# option 2: kernels (LT-file).
# character vector, one path per kernel

KFiles<-NULL 
   
# option 3: eigenvalue decomposition of kernels (*.rda)
EVDFiles<-NULL

# END of Kernels ############################################

# phenotype file

# path and filename where phenotypes are (format, see pf1)
phenFile<-'/gensys/datasets/wheat599/phenotypes.pf1'       
 
# columns in the phenFile that gives the phenotype 
#(or adjusted phneotype) to be used for trait.

traitName<-'phen_1'              

# fixed effects
covNames<-NULL
facNames<-NULL

# Parameters to the Gibbs sampler
nIter<-1500
burnIn<-500
# Prior

h2<-0.5
 
#Scale parameters are computed in mainCode.R

df0<-1e8

saveK<-TRUE

}


