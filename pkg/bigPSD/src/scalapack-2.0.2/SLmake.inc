CDEFS = -DAdd_
CC = llvm-gcc-4.2 -arch x86_64 -std=gnu99
FC = gfortran-4.2 -arch x86_64
NOOPT = -O0 
FCFLAGS = -O3 -I/opt/openmpi-1.6.4/include
CCFLAGS = -O3 -I/opt/openmpi-1.6.4/include
FCLOADER = gfortran-4.2 -arch x86_64
CCLOADER = llvm-gcc-4.2 -arch x86_64 -std=gnu99
FCLOADFLAGS = -O3 -I/opt/openmpi-1.6.4/include
CCLOADFLAGS = -O3 -I/opt/openmpi-1.6.4/include
ARCH = ar
ARCHFLAGS= cr
RANLIB = ranlib
SCALAPACKLIB = libscalapack.a
