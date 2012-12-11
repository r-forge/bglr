
#This function will read a bed file (binary file for genotypes in plink format)
read_bed=function(bed_file,bim_file,fam_file,na.strings=c("0","-9"),verbose=FALSE)
{
	#Extended map file (this gives the number of snps)
	bim=read.table(bim_file, comment.char="", as.is=TRUE, na.strings=na.strings)
	snps=as.character(bim[,2])

	#First 6 columns of ped file (this gives the number of individuals)
	fam=read.table(fam_file, comment.char="", as.is=TRUE, na.strings=na.strings)

	n=nrow(fam)
	p=length(snps)

	out=rep(0,n*p)

        if(verbose)
        {
           verbose=1
        }else 
        {   
           verbose=0
        }

	out=.C("read_bed",as.character(bed_file),as.integer(n),as.integer(p),as.integer(out),as.integer(verbose))[[4]]
        return(list(n=n,p=p,x=out))
}
