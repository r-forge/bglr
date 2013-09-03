#startup function
.onAttach = function(library, pkg)
{
  Rv = R.Version()
  if(!exists("getRversion", baseenv()) || (getRversion() < "2.12.1"))
    stop("This package requires R 2.12.1 or later")
  assign(".BGLRutils.home", file.path(library, pkg),
         pos=match("package:BGLRutils", search()))
  BGLRutils.version = "1.0 (2013-09-02)"
  assign(".BGLRutils.version", BGLRutils.version, pos=match("package:BGLRutils", search()))
  if(interactive())
  {
    packageStartupMessage(paste("Package 'BGLRutils', ", BGLRutils.version, ". ",sep=""),appendLF=TRUE)
    packageStartupMessage("Type 'help.start()' for help",appendLF=TRUE)
  }
  invisible()
}

