recodeGenotypes<-function(dataBase,naCode=9){

  MAP<-read.table(paste(dataBase,'binary/map.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)  
  
  
  alleleOne<-matrix(unlist(strsplit(MAP$alleles,split='/')),byrow=TRUE,ncol=2)[,1]
  p<-length(alleleOne)
  IDs<-scan(paste(dataBase,'/binary/IDs.txt',sep=''),what=character(),quiet=TRUE)
  n<-length(IDs)
  fileOut<-paste(dataBase,'/binary/genosInt.bin',sep='')
  unlink(fileOut)
  fileOut<-file(fileOut,open='wb')
  fileIn<-paste(dataBase,'/binary/genosText.bin',sep='')
 
  naCounts<-rep(0,p)
  alleleCounts<-rep(0,p)
  
  for(i in 1:n){
        x<-readBinGenText(filename=fileIn,whichGenotype=i,p=p)
        x<-matrix(byrow=TRUE,ncol=2,data=unlist(strsplit(x,split='')))
        isNa<-(x[,1]=='N'|x[,1]=='R'|x[,1]=='Y'|
               x[,2]=='N'|x[,2]=='R'|x[,2]=='Y'
            )
        counts<-(x[,1]==alleleOne)+(x[,2]==alleleOne)
        counts<-ifelse(x[,1]=='H'|x[,2]=='H',1,counts)
        geno<-paste(ifelse(isNa,naCode,counts),collapse='')
        writeChar(object=geno,con=fileOut,eos=NULL)
		
     	tmp<-as.integer(strsplit(geno,split='')[[1]])
  	    naCounts<-naCounts+as.integer(tmp==9)

        alleleCounts<-alleleCounts+ifelse(tmp==9,0,tmp)
		
        cat(i,"/",n,"\n")
  }
   freqNAs<-naCounts/n
   freqAlleleOne<-alleleCounts/2/(n-naCounts)
   
   MAP$freqNAs<-freqNAs
   MAP$freqAleleOne<-freqAlleleOne
   write.table(MAP,file=paste(dataBase,'binary/map.txt',sep=''),col.names=TRUE,row.names=FALSE)
   
}

## allele codes
#  Value    interpretation
#  A           A
#  C           C
#  G           G
#  T           T
#  R          A o G
#  Y          C o T
#  N         i uknown base
#  +            +
#  O          - o +
#  -          deletion
#  H          reserved
#  V          reserved