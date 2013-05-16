#This function reads a submatrix stored in a binary file
#The file contains a matrix of dimensions n x p, stored in row major order.
#rows=n
#columns=p
#from_column integer in {1,...,p}
#to_column integer in {1,...,p}

read_sub_matrix=function(n,p,from_column,to_column,file,centers,weights,center_internally=0,standard_internally=0)
{
        x=rep(0,n*(to_column-from_column+1))
	out=.C("read_sub_matrix",
    		as.double(x),
    		as.integer(n),
    		as.integer(p),
    		as.integer(from_column),
    		as.integer(to_column),
    		as.character(file),
    		as.double(centers),
    		as.double(weights),
    		as.integer(center_internally),
    		as.integer(standard_internally))
                out=out[[1]];
                dim(out)=c(n,to_column-from_column+1);
		return(out);
}

G_matrix=function(n,p,n_submatrix,file,centers,weights,center_internally=0,standard_internally=0,NPROWS=0, NPCOLS=0, MB=16, RFLAG=1, SPAWN=1)
{
   to_column=0;
   delta=as.integer(p/n_submatrix);
   
   G=matrix(0,nrow=n,ncol=n)

   for(k in 1:n_submatrix)
   {
      from_column=to_column+1;
      to_column=delta*k;
      if(k==n_submatrix) to_column=p;
      cat("Submatrix: ",k,"\n");
      cat("from_column: ",from_column,"\n");
      cat("to_column: ",to_column,"\n");
      Xtmp=read_sub_matrix(n,p,from_column,to_column,file,centers,weights,center_internally,standard_internally)
      cat("Computing... ")
      G=G+sla.tcrossprod(Xtmp,NPROWS,NPCOLS,MB,RFLAG,SPAWN)
      cat("Done\n")
   }
   return(G)
}

if(FALSE)
{       
        rm(list=ls())
        library(RScaLAPACK)

        #Example 1
	n=5416
	p=2826
	n_submatrix=3
	file="/Users/paulino/Desktop/genosInt3K.bin"
	centers=rep(0,p)
	weights=rep(1,p)
	G=G_matrix(n,p,n_submatrix,file,centers,weights,center_internally=0,standard_internally=0,NPROWS=2,NPCOLS=2,MB=16,RFLAG=1,SPAWN=1)
        load("~/Desktop/X3K.rda")
        for(i in 1:n)
        {
          for(j in 1:p)
          {
              if(X[i,j]==9) X[i,j]=0
          }
        }
        G1=tcrossprod(X)
        d=as.vector(G)-as.vector(G1)
        summary(d)
        
        # Example 2
	n=9229
	p=500568
        n_submatrix=40
	centers=rep(0,p)
	weights=rep(1,p)
	file="/Users/paulino/Desktop/genosInt.bin"
        
        G=G_matrix(n,p,n_submatrix,file,centers,weights,center_internally=0,standard_internally=0,NPROWS=2,NPCOLS=2,MB=16,RFLAG=1,SPAWN=1)
}

D_matrix=function(n,p,file)
{
}
