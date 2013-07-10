readBinGenText<-function(filename,whichGenotype,p){
   nChar<-p*2
   filein<-file(filename,open='rb')
   seek(con=filein,where=(whichGenotype-1)*(nChar))
   x<-readChar(con=filein,nchar=nChar)
   close(filein)
   return(x)
}

readBinGenInt<-function(filename,whichGenotype,p){
   nChar<-p
   filein<-file(filename,open='rb')
   seek(con=filein,where=(whichGenotype-1)*(nChar))
   x<-readChar(con=filein,nchar=nChar)
   x<-as.integer(unlist(strsplit(x,split='')))

   close(filein)
   return(x)
}

