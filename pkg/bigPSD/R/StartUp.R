#Functions executed when loading and unloading the package

.onAttach <- function(library, pkg)
{
  	Rv = R.Version()
  	if(!exists("getRversion", baseenv()) || (getRversion() < "2.13.0"))
    		stop("This package requires R 2.13.0 or later")
  	assign(".bigPSD.home", file.path(library, pkg),
         	pos=match("package:bigPSD", search()))
  	bigPSD.version = "0.6.1 (2013-05-11)"
  	assign(".bigPSD.version", bigPSD.version, pos=match("package:bigPSD", search()))
  	if(interactive())
  	{
    		packageStartupMessage(paste("Package 'bigPSD', ", bigPSD.version, ". ",sep=""),appendLF=TRUE)
    		packageStartupMessage("Type 'help.start()' for more information",appendLF=TRUE)
  	}
  	invisible()
}

.onUnload <- function(libpath)
{
	library.dynam.unload("bigPSD", libpath)
}

