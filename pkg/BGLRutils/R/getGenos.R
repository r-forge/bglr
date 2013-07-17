getGenos<-function(genDatabase,IDs=NULL,markerNames=NULL,minMAF=0,maxFreqNA=1,na.code=9,addColnames=TRUE,addRownames=TRUE)
{
 
   MAP<-read.table(file=paste(genDatabase,'/binary/map.txt',sep=''),header=TRUE)
   p<-ifelse(is.null(markerNames),nrow(MAP),length(whichMarkers))

   if(is.null(markerNames)){
      keepMarker<-rep(TRUE,p)
   }else{
      keepMarker<-MAP$rs%in%markerNames
   }
   tmp<-(MAP$freqNAs>=(1-maxFreqNA))
   keepMarker<-keepMarker&tmp
   MAF=ifelse(MAP$freqAleleOne>0.5,1-MAP$freqAleleOne,MAP$freqAleleOne)
   tmp<-MAF>=minMAF
   keepMarker<-keepMarker&tmp
   whichMarkers<-which(keepMarker)
   
   if(length(whichMarkers)==0){ stop('No marker satisfied the required MAF and frequ of NA conditions') }
  
   filein<-paste(genDatabase,'/binary/genosInt.bin',sep='')
   IDsDB<-scan(paste(genDatabase,'/binary/IDs.txt',sep=''),what=character(),quiet=TRUE)
   if(is.null(IDs)){IDs<-IDsDB }
   
   
   tmp<-IDs%in%IDsDB
   
   if(sum(tmp)==0){
      stop('None of the IDs exist in the database') 
   }else{
      if(sum(tmp)<length(IDs)){ 
	          print(paste('Warning, only  ', sum(tmp), ' out of the ', length(IDs), '  IDs appear in the database',sep='')) 
	  }
   }
   
   IDs<-IDs[which(tmp)]
   n<-length(IDs)
   X<-matrix(nrow=n,ncol=length(whichMarkers),NA)
   if(addRownames){ rownames(X)<-IDs }
   if(addColnames){ colnames(X)<-MAP$rs[whichMarkers] }
 
   cat("\n################################################## 100%% \n");
   for(i in 1:n)
   {
      #tmp<-i%in%seq(from=1,to=n,by=100)
      #if(tmp)
      #{
      #  cat('Loading genotype ',i,' (',round(100*i/n),'%)',sep='');cat('\n')
      #}

      tmp<-which(IDsDB==IDs[i])
      tmp<-readBinGenInt( filename=filein,p=p,whichGenotype=tmp)
      tmp<-tmp[whichMarkers]
      tmp<-ifelse(tmp==na.code,NA,tmp)
      X[i,]<-tmp

      progress1<-as.integer(i/n*100);
      progress2<-as.integer((i+1)/n*100)
      if(progress1%%2==0)
      {
           if(progress2%%2!=0) cat("#");
      }
   }

   cat("\n")

   return(X)
}	
