# file BGLR/methods.R
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

#BGLR: A Statistical Package for Whole-Genome Regression & Prediction
#Authors: Gustavo de los Campos & Paulino Perez Rodriguez
#Birmingaham, Alabama, 2013, 2014

print.BGLR=function(x,...)
{
  	if(!inherits(x, "BGLR")) stop("This function only works for objects of class `BGLR'\n");
}

summary.BGLR=function(object,...)
{
	tmp<-paste('--------------------> Summary of data & model <--------------------')
   	cat(tmp,'\n\n')

   	tmp<-paste(' Number of phenotypes=', sum(!is.na(object$y)))
        cat(tmp,'\n')

	cat(' Min (TRN)= ', min(object$y,na.rm=TRUE),'\n')
	cat(' Max (TRN)= ', max(object$y,na.rm=TRUE),'\n')
	cat(' Variance of phenotypes (TRN)=', round(var(object$y,na.rm=TRUE),4),'\n') 
 	cat(' Variance of weighted phenotypes (TRN)=',round(var(object$y*object$weights,na.rm=TRUE),4),'\n')
        cat(' Residual variance=',round(object$varE,4),'\n')

        n<-length(object$y)

        if(any(is.na(object$y)))
	{     
     		tst<-which(is.na(object$y))

     		cat(' N-TRN=',n-length(tst), ' N-TST=',length(tst),'\n')
     
     		cat(' Correlation TRN:\n')

     		cat('   - Standard= ',round(cor(object$y[-tst],object$yHat[-tst]),4),'\n')
  
     		cat('   - Weighted= ',round(cor((object$weights*object$y)[-tst],(object$weights*object$yHat)[-tst]),4),'\n')
   }else{
     cat(' N-TRN=',n,'  N-TST=0', '\n\n')

     cat(' Correlation TRN:\n')
     
     cat('   - Standard= ',round(cor(object$y,object$yHat),4),'\n')

     cat('   - Weighted= ',round(cor(object$weights*object$y,object$weights*object$yHat),4),'\n')

   }
  
   cat('\n')
   cat(' -- Information about the linear term -- \n')
   cat('\n')
   cat(' Model includes intercept by default\n')
   
   for(k in 1:length(object$ETA))
   {
        if(object$ETA[[k]]$model=="FIXED")
	{
		cat(" Coefficientes in ETA[",k,"] are asigned a flat prior\n")	
	}else{
		if(object$ETA[[k]]$model=="RKHS")
		{
			cat(" ETA[",k,"] assumed as a random effect normally distributed with \n covariance (or its eigendecoposition) provided by user \n");
		}else{
			cat(" Coefficientes in ETA[",k,"] modeled as in ", object$ETA[[k]]$model,"\n")
		}
	}
   } 

   cat('\n------------------------------------------------------------------\n');	
}

residuals.BGLR=function(object,...)
{
	object$y-object$yHat
}

predict.BGLR=function(object,newdata,...)
{
	object$yHat
}

effects.BGLR=function(object,...)
{
	object$ETA	
}

