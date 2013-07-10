readLT1<-function(filename,IDs=NULL,IDsG=NULL,verbose=FALSE){

  inFile<-file(description=filename,open='r')

  if(is.null(IDs))
  {
        IDs<-scan(file=inFile,what=character(),nlines=1,quiet=TRUE)
        n<-length(IDs)
        if(n<=1) stop("Since you did not provide IDs I was trying to \n use first line as ID, but it has only one field\n");
  }

  n<-length(IDs)
  G<-matrix(nrow=n,ncol=n,NA)
  colnames(G)<-IDs
  rownames(G)<-IDs

  for(i in 1:n){
    tmp<-scan(inFile,what=numeric(),nlines=1,quiet=TRUE)
    G[i,1:i]<-tmp
    G[1:i,i]<-tmp
  }
  close(inFile)
  if(!is.null(IDsG)){
    tmp<-IDs%in%IDsG
    G<-G[tmp,tmp]
  }
  return(G)
}

readBinLT1=function(filename)
{
  inFile=file(description=filename,open='rb')
  n=readBin(inFile,what=integer(),n=1)
  G=matrix(nrow=n,ncol=n,NA)
  for(i in 1:n){
    tmp=readBin(inFile,what=double(),n=i)
    G[i,1:i]=tmp
    G[1:i,i]=tmp
  }
  call_length=readBin(inFile,what=integer(),n=1)
  call=readChar(inFile,nchars=call_length)
  close(inFile)

  return(list(G=G,call=call));
}

readBinEigen=function(filename)
{
  inFile=file(description=filename,open='rb')
  n=readBin(inFile,what=integer(),n=1)
  m=readBin(inFile,what=integer(),n=1)
  vectors=matrix(nrow=n,ncol=m,NA)

  #Read eigenvectors
  for(i in 1:n)
  {
     tmp=readBin(inFile,what=double(),n=m)
     vectors[i,]=tmp
  }
  
  #Read eigenvalues
  values=readBin(inFile,what=double(),n=m)
  
  call_length=readBin(inFile,what=integer(),n=1)
  call=readChar(inFile,nchars=call_length)
  close(inFile)

  return(list(values=values,vectors=vectors,call=call));

}
