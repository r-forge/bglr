#include <stdio.h>
#include <R.h>

/*
  This function reads an arbitrary long line from a text file
  dynamically allocating memory for the buffer.
  Arguments: 
  fp: pointer to file
  nchar: number of characters that the function was able to read. If this number is negative that means that
         EOF (end of file). If nchar=0 then is a white line, if nchar>0 then this corresponds to a non empty line.
*/

#define BUF_SIZE 1000

char * read_string(FILE *fp, int *nchar)
{
      	char *buffer=NULL;
        char c;
        int count=0;
        int buf_size=BUF_SIZE;
        int flag=1;

        buffer = (char *) malloc(buf_size);
                  
        if(buffer!=NULL)
        {
            while (!feof(fp) && flag )
            {
                c = fgetc(fp);
                if(count==buf_size)
                {
                     buf_size *=2;
                     char *tmp = (char *) realloc(buffer, buf_size);
                     if(!buffer)
                     {
                        free(buffer);
                        error("cannot allocate buffer in read_string");
                     }else{
			buffer = tmp;
                     }
                }
                
                if(c=='\n')
                {
                   buffer[count]='\0';
                   flag=0;
                }else{
                   if(c!=EOF) 
                   {
                     buffer[count]=c;           //FIXME: This is weird if I do not take this into account a strange symbol appears as an extra line
                     count++;
                   }  
                }
            }
        }else
        {
          	error("Unable to allocate memory for buffer in read_string\n");
        }
        *nchar=count-1;
        return(buffer); 
}

/*
Experimental rotines for Plink support
http://pngu.mgh.harvard.edu/~purcell/plink/
*/

/*

Read ped file

The PED file is a white-space (space or tab) delimited file: the first six columns are mandatory:

     Family ID
     Individual ID
     Paternal ID
     Maternal ID
     Sex (1=male; 2=female; other=unknown)
     Phenotype

The IDs are alphanumeric: the combination of family and individual ID should uniquely identify a person. 
A PED file must have 1 and only 1 phenotype in the sixth column. The phenotype can be either a quantitative 
trait or an affection status column: PLINK will automatically detect which type (i.e. based on whether 
a value other than 0, 1, 2 or the missing genotype code is observed). 
*/

/*
WARNING: the ped_file should be absolute
         use normalizeaPath in R to avoid problems 
*/

void read_ped(char **ped_file)
{
	FILE *input;
        char *Line;
        int count=0;
        int nchar;

  	input=fopen(ped_file[0],"r");
      
        if(input!=NULL)
        {
                while(!feof(input))
                {
                        count++;
			Line=read_string(input,&nchar);
                        if(nchar>=0) Rprintf("%d \t %d \t %s\n",count,nchar,Line);
		}
                fclose(input);
        }else{
                error("It was not possible to open %s",ped_file[0]);
        }
}


/*
 Bed format,
  
 http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml

 The first 3 bytes have a special meaning. The first two are fixed, a 
 'magic number' that enables PLINK to confirm that a BED file is really 
 a BED file. That is, BED files should always start 01101100 00011011. 
 The third byte indicates whether the BED file is in SNP-major or 
 individual-major mode: a value of 00000001 indicates SNP-major (i.e. list 
 all individuals for first SNP, all individuals for second SNP, etc) 
 whereas a value of 00000000 indicates individual-major 
 (i.e. list all SNPs for the first individual, list all SNPs for 
 the second individual, etc)

 01101100=0x6c
 00011011=0x1b
 00000001=0x01

 use the xxd command to view the binary file
 
 Example:

 $xxd sample.bed
 
 0000000: 6c1b 01fa ff57 bfab ffff effe bffe ffff  l....W..........
 
 $xxd -b sample.bed
 
 0000000: 01101100 00011011 00000001 11111010 11111111 01010111  l....W


For the genotype data, each byte encodes up to four genotypes (2 bits per genoytpe). The coding is
     00  Homozygote "1"/"1"
     01  Heterozygote
     11  Homozygote "2"/"2"
     10  Missing genotype

The only slightly confusing wrinkle is that each byte is effectively read backwards. That is, if we label each of the 8 position as A to H, we would label backwards:

     01101100
     HGFEDCBA

and so the first four genotypes are read as follows:

     01101100
     HGFEDCBA

           AB   00  -- homozygote (first)
         CD     11  -- other homozygote (second)
       EF       01  -- heterozygote (third)
     GH         10  -- missing genotype (fourth)

Finally, when we reach the end of a SNP (or if in individual-mode, the end of an individual) we skip to the start of a new byte (i.e. skip any remaining bits in that byte).
It is important to remember that the files test.bim and test.fam will already have been read in, so PLINK knows how many SNPs and individuals to expect.

The routine will return a vector of dimension n*p, with the snps stacked. The vector contains integer codes:

Int code	Genotype	
0		00
1		01
2		10
3		11

Recode snp to 0,1,2 Format using allele "1" as reference

0 --> 0
1 --> 1
2 --> NA
3 --> 2

*/

void read_bed(char **bed_file, int *n, int *p, int *out, int *verbose)
{
 	FILE *input;
        unsigned char magic[3];   // First three bytes in the file
        int nbyte,i,j,k,l,bytes_read;
	char *buffer=NULL;
        unsigned char c;
        const unsigned char recode[4] = {'0', '2', '1', '3'};
  	const unsigned char mask = '\x03'; //This corresponds to 00000011 
        unsigned char code;

        input=fopen(bed_file[0],"rb");

        if(input!=NULL)
        {
           
           if(fread(magic,1,3,input)!=3)
           {
             	error("Unable to read the first 3 bytes in %s ", bed_file[0]);
           }else{
             	if(magic[0]!='\x6C' || magic[1]!='\x1B') error("%s file is not a valid .bed file (%X, %X), magic number error\n", bed_file[0],magic[0], magic[1]);
             	if(magic[2]!='\x01' && magic[2]!='\x00') error("only snp and individual major order are supported\n");
           }
           
           if(magic[2]=='\x01')
           {
              	if(*verbose) Rprintf("Start reading in snp major order...\n");
		nbyte=1+(*n-1)/4;
                if(*verbose) Rprintf("Number of bytes/snp = %d \n",nbyte);
                if(*verbose) Rprintf("Hex dump by snp \n");

                buffer = (char *) malloc(nbyte);                 
                
                if(buffer!=NULL)
                {
                        //Loop over SNPs
                        for(j=0; j<*p;j++)
                        {
                                bytes_read=fread(buffer,1,nbyte,input);
                                if(bytes_read==nbyte) 
                                {
                                        //Loop over individuals
                                        l=-1;
                                        if(*verbose) Rprintf("%d\t: ",j+1);
					for(i=0; i<nbyte;i++)
                                        {
                                          c=buffer[i];
                                          if(*verbose) Rprintf("%x ",c);
                                          for(k=0; k<4;k++)
                                          {
                                                 l++;
                                                 code = c & mask;
                                                 c=c>>2;  //Right shift (two bits)
                                                 //Some pices of information are meaningless if the number of individuals IS NOT a multiple of 4
                                                 //at the end of the snp
                                                 if(l<(*n))    
                                                 {
                                                    //Rprintf("%c",recode[code]);
                                                    out[l+j*(*n)]=recode[code]-'0';
                                                 }
                                          }
                                          if(*verbose){
                                          Rprintf(" ");
                                          if((i+1)%16==0) 
                                          {
                                          	Rprintf("\n");
                                                Rprintf("\t: ");
                                          }}
                                        }
                                        if(*verbose) Rprintf("\n");
                                }else{
                                        error("Unexpected number of bytes read from %s, expecting: %d, read: %d",bed_file[0],nbyte,bytes_read);
                                }
                        }
			free(buffer); 
  			fclose(input);
                }else{
                        error("Unable to allocate memory for buffer in read_bed\n");
                }
           }
           if(magic[2]=='\x00')
           {
		error("Individual major order not implemented yet"); 
           }  
       }else{
	     	error("It was not possible to open %s", bed_file[0]);
      }
}

/*

This rountine write a BED file in snp major order

The routine recibe a vector of dimension n*p, with the snps stacked. 
The vector contains integer codes:

Int code        Genotype        
0               00
1               01
2               10
3               11

*/

void write_bed(char **bed_file, int *n, int *p, int *out)
{
	FILE *output;
        const unsigned char recode[4] = {0x00, 0x02, 0x01, 0x03};
        unsigned char byte,mask; 
        int i,j,l;

        output=fopen(bed_file[0],"wb");
        
        if(fopen!=NULL)
        {		
                /*
                Write the magic number
                */
                fputc(0x6C, output); 
                fputc(0x1B, output);
                 
                /* 
                 Order, 
                 right now only SNP major order is supported 
                */
                fputc(0x01, output);

   		//Loop over SNPs
                for(j=0; j<*p;j++)
                {
			//Loop over individuals, 
                        byte=0x00;
                        l=-1;
                        for(i=0; i<*n; i++)
                        { 
                                l++;
                                mask=recode[out[i+j*(*n)]] << 2 * l;
                                byte = byte | mask ;
                                if((i+1)%4==0)
                                {
                                   fputc(byte,output);
                                   byte=0x00;
                                   l=-1;
                                }
                        }
                        //This includes n<4 or n that is not multiple of 4
                        if((*n)%4!=0) fputc(byte,output);
		}       
                fclose(output);
                
        }else{
		error("It was not possible to open %s", bed_file[0]);
       }
}

