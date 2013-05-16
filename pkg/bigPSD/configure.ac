dnl Process this file with autoconf to produce a configure script
AC_INIT(DESCRIPTION)

dnl Determine R_HOME.
: ${R_HOME=‘R RHOME‘}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

R_SCMD="${R_HOME}/bin/Rscript -e"

dnl Set EXT_LDFLAGS, otherwise default LDFLAGS, OBJS, and LIBS will be set.
AC_ARG_VAR(EXT_LDFLAGS, [ external linking flags of external libraries ])
AS_VAR_SET_IF(EXT_LDFLAGS, [], [ AS_VAR_SET(EXT_LDFLAGS, "") ])
if test "X$EXT_LDFLAGS" = "X" ; then
  EXT_LDFLAGS="\$(LAPACK_LIBS) \$(BLAS_LIBS)"
fi

dnl C and Fortran interface
dnl CDEFS should be one of "-DAdd_", "-DNoChange", "-Df77IsF2C", or "-DUpCase"
dnl which is determined by "../inst/test_intface/00_make.sh".
dnl For GNU make, gcc, and gfortran, the default is "-DAdd_".
AC_ARG_VAR(CDEFS, [ intface of C and Fortran ])
AS_VAR_SET_IF(CDEFS, [], [ AS_VAR_SET(CDEFS, "") ])
if test "X$CDEFS" = "X" ; then
  CDEFS=-D`${R_SCMD} "invisible(system('cd ./inst/test_intface/;./00_make.sh'))"`
  CDEFS_CHECK=TRUE
fi
if test "X$CDEFS" = "X" ; then
  CDEFS=-DAdd_
  CDEFS_CHECK=FALSE
fi

dnl Check MPI in path.
AC_CHECK_PROG(MPIRUN, mpirun, mpirun, "F")
AC_CHECK_PROG(MPIEXEC, mpiexec, mpiexec, "F")
AC_CHECK_PROG(ORTERUN, orterun, orterun, "F")
AC_CHECK_PROG(SED, sed, sed, "F")
AC_CHECK_PROG(MPICC, mpicc, mpicc, "F")
AC_CHECK_PROG(OMPI_INFO, ompi_info, ompi_info, "F")

MPITYPE="UNKNOWN"
dnl For OpenMPI
if test -z "${MPI_INCLUDE_PATH}" ; then
  if test "$MPITYPE" = "UNKNOWN" -o "$MPITYPE" = "OPENMPI" ; then
    if test "$SED" != "F" -a "$MPICC" != "F" -a "$OMPI_INFO" != "F" ; then
      echo "found sed, mpicc, and ompi_info ..."
      TMP_INC="F"
      TMP_LIB="F"

      TMP_INC_DIRS=`mpicc --showme:incdirs`
      echo ">> TMP_INC_DIRS = ${TMP_INC_DIRS}"

      TMP_IFS="${IFS}"; IFS=" "
      for dir in ${TMP_INC_DIRS}; do
        echo "checking ${dir} ..."
        if test -f "${dir}/mpi.h" ; then
          echo "found ${dir}/mpi.h ..."
          TMP_INC=${dir}
          break
        fi
      done
      IFS="${TMP_IFS}"

      TMP_LIB_DIRS=`mpicc --showme:libdirs`
      echo ">> TMP_LIB_DIRS = ${TMP_LIB_DIRS}"

      TMP_IFS="${IFS}"; IFS=" "
      for dir in ${TMP_LIB_DIRS}; do
        echo "checking ${dir} ..." 
        if test -f "${dir}/libmpi.so" ; then
          echo "found ${dir}/libmpi.so ..." 
          TMP_LIB=${dir}
          break
        fi
        dnl For Mac OS X
        if test -f "${dir}/libmpi.dylib" ; then
          echo "found ${dir}/libmpi.dylib ..."
          TMP_LIB=${dir}
          break
        fi
      done
      IFS="${TMP_IFS}"

      if test "${TMP_INC}" != "F" -a "${TMP_LIB}" != "F" ; then
        echo "found mpi.h and libmpi.so ..."
        MPITYPE="OPENMPI"
        MPI_INCLUDE_PATH="-I${TMP_INC}"
        MPI_LIB_PATH="-L${TMP_LIB}"
        AC_CHECK_LIB(util, openpty, [ MPI_LIBS="$MPI_LIBS -lutil" ])
	AC_CHECK_LIB(pthread, main, [ MPI_LIBS="$MPI_LIBS -lpthread" ])
        LIBS="${MPI_LIB_PATH} -lmpi ${MPI_LIBS}"
      fi
    fi
  fi
fi

dnl Report
echo " "
echo "****************** Results of bigPSD package configure *****************"
echo " "
echo ">> EXT_LDFLAGS = ${EXT_LDFLAGS}"
echo ">> CDEFS = ${CDEFS}"
echo ">> CDEFS_CHECK = ${CDEFS_CHECK}"
echo " "
echo "*************************************************************************"
echo " "

dnl Start to substitute templates
dnl AC_SUBST(MPI_INCLUDE_PATH)
AC_SUBST(LIBS)
AC_SUBST(MPI_INCLUDE_PATH)
AC_SUBST(CDEFS)
AC_OUTPUT(src/Makefile)

