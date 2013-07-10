#This function reads a submatrix stored in a binary files
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
                as.character(normalizePath(file)),
                as.double(centers),
                as.double(weights),
                as.integer(center_internally),
                as.integer(standard_internally))
                out=out[[1]];
                dim(out)=c(n,to_column-from_column+1);
                return(out);
}

G_matrix=function(n,p,n_submatrix,file,centers=NULL,weights=NULL,center_internally=1,standard_internally=1)
{
   to_column=0;
   delta=as.integer(p/n_submatrix);

   if(is.null(centers))
   {
     centers=rep(0,p)
   }

   if(is.null(weights))
   {
     weights=rep(1,p)
   }

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
      G=G+tcrossprod(Xtmp)
      cat("Done\n")
   }
   return(G)
}


