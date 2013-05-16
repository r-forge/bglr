#
# RSCALAPACK/src/Makefile
#
CDEFS=-DAdd_
LIBS=-L/opt/openmpi-1.6.4/lib -lmpi  -lutil -lpthread scalapack-2.0.2/libscalapack.a  $(LAPACK_LIBS) $(BLAS_LIBS)  
MPI_INCLUDE_PATH=-I/opt/openmpi-1.6.4/include
CFLAGS= $(MPI_INCLUDE_PATH) -DOMPI -fPIC -O3 -pipe -g 
INTRINSICLIBS=-lgfortran -ldl -lm
exec_dir=../exec

PA_OBJECTS=ParallelAgent.o PAdistData.o 

CR_OBJECTS=CRscalapack.o callpdgesv.o callpdgeqrf.o callpdgesvd.o callpdgemm.o\
		callpdpotrf.o callpdpotri.o callpdsyevd.o CRcollectData.o \
		CRdistData.o callpdgetcp.o

OTHER_OBJECTS= Kernel_utils.o Plink_utils.o 

all: libscalapack.a bigPSD.so CRDriver

#Here we generate a SLmake.inc file
#This file is used when compiling ScaLAPACK
#Refer to ScaLAPACK documentation for further information

libscalapack.a: 
	rm -f scalapack-2.0.2/SLmake.inc
	@echo "Generating SLmake.inc"
	@echo "CDEFS = $(CDEFS)" > scalapack-2.0.2/SLmake.inc
	@echo "CC = $(CC)" >> scalapack-2.0.2/SLmake.inc
	@echo "FC = $(FC)" >> scalapack-2.0.2/SLmake.inc
	@echo "NOOPT = -O0 $(MPI_INCUDE_PATH)" >> scalapack-2.0.2/SLmake.inc
	@echo "FCFLAGS = -O3 $(MPI_INCLUDE_PATH)" >> scalapack-2.0.2/SLmake.inc
	@echo "CCFLAGS = -O3 $(MPI_INCLUDE_PATH)" >> scalapack-2.0.2/SLmake.inc
	@echo "FCLOADER = $(FC)" >> scalapack-2.0.2/SLmake.inc
	@echo "CCLOADER = $(CC)" >> scalapack-2.0.2/SLmake.inc
	@echo "FCLOADFLAGS = -O3 $(MPI_INCLUDE_PATH)" >> scalapack-2.0.2/SLmake.inc
	@echo "CCLOADFLAGS = -O3 $(MPI_INCLUDE_PATH)" >> scalapack-2.0.2/SLmake.inc
	@echo "ARCH = ar" >>  scalapack-2.0.2/SLmake.inc
	@echo "ARCHFLAGS= cr" >> scalapack-2.0.2/SLmake.inc
	@echo "RANLIB = ranlib" >> scalapack-2.0.2/SLmake.inc
	@echo "SCALAPACKLIB = libscalapack.a" >> scalapack-2.0.2/SLmake.inc
	@echo "Done"
	( cd scalapack-2.0.2; $(MAKE) )
        
bigPSD.so: $(PA_OBJECTS) $(OTHER_OBJECTS)
	$(SHLIB_LINK) $(PA_OBJECTS) $(OTHER_OBJECTS) $(CFLAGS) $(LIBS) -o $@


CRDriver: CRDriver.o $(CR_OBJECTS)
	$(CC) CRDriver.o $(CR_OBJECTS) $(LIBS) $(CFLAGS) $(INTRINSICLIBS) -o $@
	cp CRDriver $(exec_dir)

clean:
	rm -rf CRDriver $(PA_OBJECTS) $(CR_OBJECTS) $(OTHER_OBJECTS) CRDriver.o RScaLAPACK.so $(exec_dir)/CRDriver
	(cd scalapack-2.0.2; $(MAKE) clean)
