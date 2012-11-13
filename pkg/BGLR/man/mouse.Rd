\name{mouse}
\alias{mouse}
\title{mouse dataset}
\description{

The mice data comes from an experiment carried out to detect and locate QTLs for complex traits in a 
mice population (Valdar et al. 2006a; 2006b). This data has already been analyzed for comparing genome-assisted 
genetic evaluation methods (Legarra et al. 2008). The data file consists of 1882 individuals 
(168 full-sib families), each genotyped for 10,946 polymorphic markers. The trait here 
here is body mass index (BMI), pre-corrected by body weight, season, month and day. Mice were 
housed in 359 cages; on average, each litter was allocated into 2.84 cages.

}

\usage{
  data(mouse)
}

\format{

 Matrix mouse.A contains the pedigree. The vector mouse.cageID contains the cage id,
 the vector mouse.ID contains an in for each mouse. The matrix mouse.X contains the markes information and 
 mouse.y corresponds to the body mass index.  

}

\source{
   \url{http://gscan.well.ox.ac.uk} 
 }

\references{

Legarra A., Robert-Granie, E. Manfredi, and J. M. Elsen, 2008 Performance of genomic selection in mice. Genetics 180:611-618.

Valdar, W., L. C. Solberg, D. Gauguier, S. Burnett, P. Klenerman et al., 2006a Genome-wide genetic association of complex traits in heterogeneous stock mice. Nat. Genet. 38:879-887.

Valdar, W., L. C. Solberg, D. Gauguier, W. O. Cookson, J. N. P. Rawlis et al., 2006b Genetic and environmental effects on complex traits in mice. Genetics, 174:959-984.

}


\keyword{datasets}
