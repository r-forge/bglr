\name{D_matrix}
\alias{D_matrix}
\title{D_matrix}
\description{
  Computes Euclidean distance matrices allowing the user to specify several options, for 
  example weighting and scaling.
  
  The computations are done in parallel. The number of cpu cores can be specified in the UNIX 
  commmand line with

  $export omp_num_threads=8
}

\usage{
  D_matrix(dataBase,fileOut,
           weights,scale,gensys.loc)
}

\arguments{
  \item{dataBase}{Location of the dataBase.}
  \item{fileOut}{The name of the output file.}
  \item{weights}{Vector for weighting the markers, if NULL it will not weight the markers. The defaul value is NULL.}
  \item{scale}{Scalar, it will divide the G matrix by using this number, default value 1.}
  \item{gensys.loc}{location of the gensys software. The system will use the environment variable gensys_home to obtain this 
                   information if it is not provided by the user.}
}

\author{Paulino Perez, Gustavo de los Campos, John Hickey}

\examples{

#D_matrix(dataBase='~/Desktop/genDatabases/trainingwheatGBS')

}
