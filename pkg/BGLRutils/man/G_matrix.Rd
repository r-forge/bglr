\name{G_matrix}
\alias{G_matrix}
\title{G_matrix}
\description{
  Computes genomic relationship matrices allowing the user to specify several options, for 
  example centering the markers, weighting and scaling.
  
  The computations are done in parallel. The number of cpu cores can be specified in the UNIX 
  commmand line with

  $export omp_num_threads=8
}

\usage{
  G_matrix(dataBase,fileOut,centers,
           weights,scale,gensys.loc)
}

\arguments{
  \item{dataBase}{Location of the dataBase.}
  \item{fileOut}{The name of the output file.}
  \item{center}{Vector for centering the markers, if NULL it will not center the markers. The default value is NULL.}
  \item{weights}{Vector for weighting the markers, if NULL it will not weight the markers. The defaul value is NULL.}
  \item{scale}{Scalar, it will divide the G matrix by using this number, default value 1.}
  \item{gensys.loc}{location of the gensys software. The system will use the environment variable gensys_home to obtain this 
                   information if it is not provided by the user.}
}

\author{Paulino Perez, Gustavo de los Campos, John Hickey}

\examples{

#G_matrix(dataBase='~/Desktop/genDatabases/trainingwheatGBS')


}
