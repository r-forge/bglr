\name{recodeGenotypes}
\alias{recodeGenotypes}
\title{recodeGenotypes}
\description{
  This function recodes a binary file with allele codes to 
  0,1,2 and NA format.

 \tabular{ll}{
     \bold{Value} \tab \bold{Interpretation} \cr
        A  \tab         A \cr
        C  \tab         C \cr
        G  \tab         G \cr
        T  \tab         T \cr
        R  \tab       A o G \cr
        Y  \tab        C o T \cr
        N  \tab        uknown base \cr
        +  \tab          + \cr
        O  \tab        - o + \cr
        -  \tab        deletion \cr
        H  \tab        reserved \cr
        V  \tab        reserved
  }
}

\usage{
  recodeGenotypes(dataBase,naCode=9)
}


\arguments{
  \item{dataBase}{Location of the dataBase.}
  \item{naCode}{missing value indicator in the dataBase, the default value is 9.}
}

\author{Paulino Perez, Gustavo de los Campos, John Hickey}
