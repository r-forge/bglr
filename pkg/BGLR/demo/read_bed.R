library(BGLR)
bed = system.file("extdata/sample.bed", package="BGLR")
n=120
p=20
out=rep(0,n*p)

out=.C("read_bed",as.character(bed),as.integer(n),as.integer(p),as.integer(out))[[4]]

#Recode snp to 0,1,2 format using allele 1
# 0 --> 0
# 1 --> 1
# 2 --> NA
# 3 --> 2

out[out==2]=NA
out[out==3]=2

X=matrix(out,nrow=p,ncol=n,byrow=TRUE)
X=t(X)

