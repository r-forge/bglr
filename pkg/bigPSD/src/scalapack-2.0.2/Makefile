############################################################################
#
#  Program:         ScaLAPACK
#
#  Module:          Makefile
#
#  Purpose:         Top-level Makefile
#
#  Creation date:   March 20, 1995
#
#  Modified:        February 15, 2000
#
#  Send bug reports, comments or suggestions to scalapack@cs.utk.edu
#
############################################################################

include SLmake.inc

PRECISIONS = double

############################################################################
#
#  The library can be set up to include routines for any combination of the
#  four PRECISIONS.  First, modify the ARCH, ARCHFLAGS, RANLIB, F77, CC,
#  F77FLAGS, CCFLAGS, F77LOADER, CCLOADER, F77LOADFLAGS, CCLOADFLAGS and
#  CDEFS definitions in SLmake.inc to match your library archiver, compiler
#  and the options to be used.
#
#  The command
#       make
#  without any arguments creates the library of precisions defined by the
#  environment variable PRECISIONS as well as the corresponding testing
#  executables,
#       make lib
#  creates only the library,
#       make exe
#  creates only the testing executables.
#       make example
#  creates only the example
#
#  The name of the library is defined in the file called SLmake.inc and
#  is created at this directory level.
#
#  To remove the object files after the library and testing executables
#  are created, enter
#       make clean
#
############################################################################

all: blacslib toolslib pblaslib redistlib scalapacklib

clean: cleanlib

blacslib:
	( cd BLACS; $(MAKE) lib )

pblaslib:
	( cd PBLAS/SRC; $(MAKE) $(PRECISIONS) )

redistlib:
	( cd REDIST/SRC; $(MAKE) integer $(PRECISIONS) )

scalapacklib:
	( cd SRC; $(MAKE) $(PRECISIONS) )

toolslib:
	( cd TOOLS; $(MAKE) $(PRECISIONS) )

cleanlib:
	( cd BLACS; $(MAKE) clean )
	( cd PBLAS/SRC; $(MAKE) clean )
	( cd SRC; $(MAKE) clean )
	( cd TOOLS; $(MAKE) clean )
	( cd REDIST/SRC; $(MAKE) clean )
	( rm -f $(SCALAPACKLIB) )

