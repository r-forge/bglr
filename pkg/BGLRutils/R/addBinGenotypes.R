addBinGenotypes<-function(fileIn,dataBase,nrow,ncol,skip)
{
  #ID file
  IDs<-scan(fileIn,what=character(),nlines=1,quiet=TRUE)[-(1:skip)]
  tmp<-'IDs.txt'%in%list.files(paste(dataBase,'/binary',sep=''))
  if(tmp)
  {
    oldIDs<-scan(file=paste(dataBase,'/binary/IDs.txt',sep=''),what=character(),quiet=TRUE)
    if(any(IDs%in%oldIDs))
    {
        stop('Some of the IDs in fileIn are already listed in the database')
    }else{ 
        rm(oldIDs) 
	gc()
    }
  }
  
  IDFile<-paste(dataBase,'/binary/IDs.txt',sep='')
  write(x=IDs,file=IDFile,append=TRUE)

  #MAP
  hasMap<-'map.txt'%in%list.files(paste(dataBase,'/binary',sep=''))
  
  if(hasMap)
  {
    mapFile<-paste(dataBase,'binary/map2.txt',sep='')
    fields<- paste(strsplit(paste(paste(1:(skip-1),',',collapse=''),skip,collapse=''),split=' ')[[1]],collapse='')
    cmd<-paste(' cut ',fileIn,' -f ',fields, ' -s \t > ', mapFile,sep='')
    print(cmd)
    #system(cmd)
    
    map1<-read.table(file=paste(dataBase,'binary/map.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)
    map2<-read.table(file=paste(dataBase,'binary/map2.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)
    unlink(paste(dataBase,'binary/map2.txt',sep=''))
    if(all.equal(map1,map2)==TRUE)
    {
	#Do nothing
    }else{
	stop('Existing map does not match the map of the new data-file')
    }
  }

  mapFile<-paste(dataBase,'binary/map.txt',sep='')
  map1<-read.table(file=fileIn,
                   colClasses = c(rep("character",skip),rep("NULL",ncol-skip)))  
  write.table(file=mapFile,map1,row.names=FALSE,col.names=FALSE)

  #transpose in binary format
  fileOut=paste(normalizePath(dataBase),"/binary/genosText.bin",sep="")

  out=.C("read_write_gbs",
       as.character(normalizePath(fileIn)),
       as.integer(nrow),
       as.integer(ncol),
       as.integer(skip),
       as.character(fileOut))

  cat("File ",fileOut,"created succesfully\n");
}
