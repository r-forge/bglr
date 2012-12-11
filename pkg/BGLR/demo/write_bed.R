library(BGLR)

bed_file = system.file("extdata/sample.bed", package="BGLR")

#Extended map file (this gives the number of snps)
bim_file = system.file("extdata/sample.bim", package="BGLR")

#First 6 columns of ped file (this gives the number of individuals)
fam_file = system.file("extdata/sample.fam", package="BGLR")

out=read_bed(bed_file=bed_file,bim_file=bim_file,fam_file=fam_file,verbose=TRUE)

.C("write_bed",as.character("/tmp/test.bed"), as.integer(out$n), as.integer(out$p), as.integer(out$x))
