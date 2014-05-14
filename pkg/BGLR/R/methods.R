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
	object
}

residuals.BGLR=function(object,...)
{

}

predict.BGLR=function(object,newdata,...)
{

}

effects.BGLR=function(object,...)
{

}

