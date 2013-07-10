#include<R.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

void read_matrix(char **X, int rows, int columns, int skip_columns, char *Input_file)
{
    char *Buffer;
    char *token=NULL;
    size_t nbytes=1000; 
    FILE *ptr;
    int i,j,k; 

    Buffer= (char *) malloc (nbytes + 1);
   
    if(Buffer==NULL)
    {
        error("Unable to allocate space for buffer\n");
    }
    
    ptr=fopen(Input_file,"r");
    if(ptr!=NULL)
    {
    	i=0;
	Rprintf("Loading GBS data...");
       	while(getline(&Buffer, &nbytes, ptr)>0)
    	{
	        //printf("Line=%d\n",i);
                /*Skip the first line(s)*/
	    	if(i>0)
		{
		  		token=strtok(Buffer,"\t");
		  		j=0;
		  		k=0;
		  		while(token!=NULL)
		  		{
		        		if(j>(skip_columns-1))
					{
			  			//printf("k=%d\n",k);
			  			X[i-1][2*k]=token[0];
			  			X[i-1][2*k+1]=token[0];
			  			k++;
			  			//printf("%c\n",token[0]);
					}
					token = strtok(NULL,"\t");
                			j++;
		  		}
		}
      		i++;
    	}

    	Rprintf("Done\n");
        if(i!=rows)
        {
                 error("The file has MORE/LESS lines that those indicated by the user\n");
        }
        fclose(ptr);
    }else{
        error("Unable to open input file with incidence matrix\n");

    }	
}

void read_write_gbs(char **input_file, int *n, int *p, int *skip_columns, char **output_file)
{
	const int skip_rows=1;
	int i;
	int j;
	char **X;
	char space=' ';
		
	FILE *out;                      //Original file in binary format
		
	/*
	   Reserve space
	*/

    	X = (char**)malloc((*n-skip_rows) * sizeof(char*));
    	for (j = 0 ; j< (*n-skip_rows); j++)
	{
    		X[j] = (char*)malloc(2*((*p)-(*skip_columns))* sizeof(char));
	}
		
	read_matrix(X,*n,*p, *skip_columns, input_file[0]);
	
	/*
	 * Transpose of the file
	 */
	 
	out=fopen(output_file[0],"a+b");
	
	if(out!=NULL)
	{
        Rprintf("Writing original data in transposed binary format...");
	    for(j=0; j<(*p-*skip_columns); j++)
	    {
	       for(i=0; i<(*n-skip_rows); i++)
	       {
	          fwrite(&X[i][2*j],sizeof(char),1,out);
	          fwrite(&X[i][2*j+1],sizeof(char),1,out);
	       }
	       //fwrite(&space, sizeof(char),1,out);
	    }
	  fclose(out);
          Rprintf("Done\n");
	}else{
	  error("Unable to open file for writting in binary format...\n");
	}
	
	/*
	  Free space
	*/
	
	for (j = 0; j < (*n-skip_rows);  j++)
	  free(X[j]);                 // STEP 1: DELETE THE COLUMNS
	 free(X);                      // STEP 2: DELETE THE ROWS
}

/*
This function reads a submatrix stored in a binary file
The file contains a matrix of dimensions n x p, stored in row major order.
rows=n
columns=p
from_column integer in {1,...,p}
to_column integer in {1,...,p}
*/

/*

Discuss this:
How to deal with missing values?
Right now we are keeping out of the calculations cells with missing values
Mean and standard deviations are computed removing the missing values
Alternative: Impute missing values using naive imputation

*/
void read_sub_matrix(double *X, int *rows, int *columns, int *from_column, int *to_column,char **Input_file, double *centers, double *weights, int *center_internally, int *standard_internally)
{
	int seek;
	int i, j;
        FILE *ptr;
        char *xi;
        double x;
        double sum_x;
        double sum_x2;
        double *d;
        int n;
        off_t file_length;
        off_t offset;
        int progress1, progress2;

	Rprintf("Loading incidence matrix...\n");
        Rprintf("################################################## 100%% \n");

	if(*from_column <1 || *from_column>*columns)
        {
		error("Function read_submatrix: from_column out of range, current value is %d, allowed values 1,...,%d\n",*from_column,*columns);
        }
        
	if(*to_column>*columns || *to_column<=*from_column)
        {
		error("Function read_submatrix: to_column out of range, current value is %d, allowed values %d,...,%d\n",*to_column,*from_column+1,*columns);
	}

        xi = (char *) malloc((*to_column-*from_column+1)*sizeof(char));
       
 	if(xi==NULL)
      	{
        	error("Function read_sub_matrix: Unable to allocate memory for xi\n");
      	}

        d=(double *) malloc((*to_column-*from_column+1)*sizeof(double));

        if(d==NULL)
        {
		error("Function read_submatrix: Unable to allocate memory for d\n");
        }

        for(j=0; j<(*to_column - *from_column+1);j++) d[j]=1.0;

        ptr=fopen(Input_file[0],"rb");
        if(ptr!=NULL)
    	{
                fseeko(ptr, 0, SEEK_END);
                file_length = ftello(ptr);

        	for(i=0; i<*rows;i++)
        	{
                        offset=(off_t)i*(*columns) + *from_column-1;
                        if(offset>file_length) 
                        {
                            error("Function read_submatrix: fseeko64 error\n You are trying to jump beyond the end of the file, check n and p and verify that %s is not corrupted \n",Input_file[0]);
                            
                        }

                        seek=fseeko(ptr,offset,SEEK_SET);
                        if(seek==0)
                        {
          			fread(xi,sizeof(char),(*to_column - *from_column+1),ptr);
          			for(j=0; j<(*to_column - *from_column + 1);j++)
          			{
                                        x=xi[j]-'0';
                                        X[i+(j*(*rows))]=x;     //row + column*NUMROWS Column major order 
          			}
			}else{
				error("Function read_sub_matrix: fseeko64 error\n");
		       }

                       progress1=(int)((double)(i+1)/(*rows)*100);
                       progress2=(int)((double)(i+2)/(*rows)*100);

                       if(progress1%2==0)
                       {
                             if(progress2%2!=0) Rprintf("#");
                       }
                       
        	}
                
                /*
                Computing mean and standard deviations if necessary
                */  

                if(*center_internally || *standard_internally)
                {
			Rprintf("\nComputing  mean and standard deviations...\n");
                        Rprintf("################################################## 100%% \n");
                
                	for(j=0; j<(*to_column - *from_column+1);j++)
                	{
                    		sum_x=0;
                                sum_x2=0;
                                n=0;

                    		for(i=0; i<*rows;i++)
                    		{
                       			x=X[i+(j*(*rows))];
                       			if((int)x!=9)
                                        {
                                           if(*center_internally) sum_x+=x;
                                           if(*standard_internally) sum_x2+=x*x;
                                           n++;
                                        }
                       		}

                                if(*center_internally) centers[*from_column-1+j]=sum_x/n;
                                if(*standard_internally) d[j]=sqrt((sum_x2-sum_x*sum_x/(n))/(n-1)); 

                                progress1=(int)((double)(j+1)/(*to_column - *from_column+1)*100);
                                progress2=(int)((double)(j+2)/(*to_column - *from_column+1)*100);

                       		if(progress1%2==0)
                       		{
                             		if(progress2%2!=0) Rprintf("#");
                       		}
                       	 }
                }
                
                
                if(*standard_internally)
                {
                   for(j=0; j<(*to_column - *from_column+1);j++)
                   {
                         if(d[j]==0) Rprintf("\nWarning, marker %d is monomorphic and will be out of the computations",*from_column+j);
                   }
                }

                Rprintf("\nWeighting, centering and standardizing...\n");
                Rprintf("################################################## 100%% \n");

                for(i=0;i<*rows;i++)
		{
			for(j=0;j<(*to_column - *from_column+1);j++)
			{
			   x=X[i+(j*(*rows))];
			   if((int)x!=9 && d[j]!=0)
			   {
				x=(x-centers[*from_column-1+j])*weights[*from_column-1+j]/d[j];
                           }else{
				x=0;
                           }
                           X[i+(j*(*rows))]=x;
			}

		        progress1=(int)((double)(i+1)/(*rows)*100);
                        progress2=(int)((double)(i+2)/(*rows)*100);

                        if(progress1%2==0)
                        {
                              if(progress2%2!=0) Rprintf("#");
                        }
		}

        	fclose(ptr);
    	}else{
        	error("Function read_sub_matrix: Unable to open input file with incidence matrix\n");
    	}

        Rprintf("\n");

		free(xi);
}
