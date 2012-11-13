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
