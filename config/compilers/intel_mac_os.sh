omp=
incs=

mpi_libs=
misc_libs=

echo "Intel Fortran Compiler on Apple MacOS"

# MUST RUN Environment Setup prior to compilation	
# source /opt/intel/composerxe/bin/compilervars.csh intel64
# for Intel MKL library
# source /opt/intel/composerxe-2011.0.085/mkl/bin/mklvars.csh intel64	  

MODULE_CODE=0

# Add some additinal flags to the object directory
DPO=${DPO_BASE}/${DPO}_INTEL_MacOS/
if test ! -d ${DPO}; then mkdir ${DPO}; fi

# Set OpenMP flags.
if test ${USE_SMP} = 1; then omp="-openmp"; fi

# Set MPI flags.
if test ${USE_DMP} = 1; then
  FORTRAN_CMD=mpif90
  LINK_CMD=mpif90
else
  FORTRAN_CMD=ifort
  LINK_CMD=ifort
fi

# Set the MPI include path.  This must occur after
# the Fortran command is defined.
SET_MPI_INCLUDE
if test ${USE_DMP} = 1; then
  incs=${incs}" -I${MPI_INCLUDE_PATH}"
fi


# --> Verify compiler is in $PATH <-- #


# Set library paths:
ode="${DPO}odepack.a"
blas="${DPO}blas90.a"
dgtsv="${DPO}dgtsv90.a"

## Check if Intel MKL is available and accordingly set library options
if [ -z "$MKLROOT" ]
  then
    #echo "No argument supplied"
    mkl_libs="${blas} ${dgtsv}"
  else
    mkl_libs="-L ${MKLROOT}/lib -lmkl_intel -lmkl_sequential -lmkl_core"
fi

LIB_FLAGS="${ode} ${mkl_libs} ${misc_libs} ${mpi_libs} "

# Debug flags for Intel Fortran
dbg=
if test ${USE_DEBUG} = 1; then dbg="-g"; fi

# Base flags for Intel Fortran Linux compiler
#common="-c -I. -convert big_endian -assume byterecl -diag-disable remark"
#AIKE from make_mfix.macos
common="-c -I. -w -w95 -i_dynamic -convert big_endian -assume byterecl -diag-disable remark -cpp"
#common_d="-convert big_endian -assume byterecl -O0 -fpe0 -CB -traceback"

case $OPT in
  0)echo "Setting flags for debugging."
    dbg="-traceback -check all -fpe:0 -fp-model precise"
    FORT_FLAGS="${omp} ${incs} ${common} -FR ${dbg} -O0 -g "
    FORT_FLAGS3="${common} ${dbg} ${incs} -O0 -g "
    LINK_FLAGS="${omp} -g";;

  1)echo "Setting flags for low optimization."
    FORT_FLAGS="${omp} ${incs} ${common} -FR -O1 ${dbg} "
    FORT_FLAGS3="${common} ${incs} -O1 ${dbg} "
    LINK_FLAGS="${omp} ${dbg}";;

  2)echo "Setting flags for medium optimization."
    FORT_FLAGS="${omp} ${common} ${incs} -FR -O2 ${dbg}"
    FORT_FLAGS3="${common} -O1 ${dbg}"
    LINK_FLAGS="${omp} ${dbg}";;

  3)echo "Setting flags for high optimization."
    FORT_FLAGS="${omp} ${incs} ${common} -FR -O3"
    FORT_FLAGS3="${common} -O2"
    LINK_FLAGS="${omp}";;

  *)echo "Unsupported optimization level."
    exit;;
esac
