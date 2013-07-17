#Euclidean distance matrix
#D[i,j]=G[i,i]+G[j,j]-2G[i,j];

GtoD=function(G)
{
   n=nrow(G)
   D=matrix(0,nrow=n,ncol=n)
   for(i in 1:n)
   {
        cat(i,"/",n,"\n")
        for(j in i:n)
        {
           D[i,j]=G[i,i]+G[j,j]-2*G[i,j]
           D[j,i]=D[i,j]
        }
   }
  return(D)
}

