omp=
incs=

mpi_libs=
misc_libs=

echo "HOPPER@NERSC :: GNU gfortran Compiler on Cray XE6"
  
MODDIRPREFIX="-J"

# Add some additinal flags to the object directory
DPO=${DPO_BASE}/${DPO}_NERSC_HOPPER_GNU/
if test ! -d ${DPO}; then mkdir ${DPO}; fi

# Set OpenMP flags.
if test ${USE_SMP} = 1; then omp="-fopenmp"; fi

# Set MPI flags.
if test ${USE_DMP} = 1; then
  FORTRAN_CMD=ftn
  LINK_CMD=ftn
else
  FORTRAN_CMD=ftn
  LINK_CMD=ftn
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

# Debug flags for PGI Fortran
dbg=""
if test ${USE_DEBUG} = 1; then dbg="-g"; fi

# Base flags for Intel Fortran Linux compiler
#common="-c -I. -convert big_endian -assume byterecl -diag-disable remark"
common="-c -fconvert='big-endian' -cpp"

# Optimization flags for level 4
#optim="-V -fast -Mipa=fast,inline -Msmartalloc -Mfprelaxed -Mstack_arrays"
optim=""

case $OPT in
  0)echo "Setting flags for debugging."
    dbg="-g -fbounds-check"
    FORT_FLAGS="${omp} ${incs} ${common} ${dbg} -O0 -ffree-form -ffree-line-length-0"
    FORT_FLAGS3="${common} ${dbg} ${incs} -O0 "
    LINK_FLAGS="${omp} -g";;

  1)echo "Setting flags for low optimization."
    FORT_FLAGS="${omp} ${incs} ${common} -O1 ${dbg} -ffree-form -ffree-line-length-0"
    FORT_FLAGS3="${common} ${incs} -O1 ${dbg} "
    LINK_FLAGS="${omp} ${dbg}";;

  2)echo "Setting flags for medium optimization."
    FORT_FLAGS="${omp} ${common} ${incs} -O2 ${dbg} -ffree-form -ffree-line-length-0"
    FORT_FLAGS3="${common} -O1 ${dbg}"
    LINK_FLAGS="${omp} ${dbg}";;

  3)echo "Setting flags for high optimization."
    FORT_FLAGS="${omp} ${incs} ${common} -O3 -ffree-form -ffree-line-length-0"
    FORT_FLAGS3="${common} -O2"
    LINK_FLAGS="${omp}";;

  4)echo "Setting flags for Polyhedron based optimization flags."
    FORT_FLAGS="${omp} ${incs} ${common} ${optim} -ffree-form -ffree-line-length-0"
    FORT_FLAGS3="${common} -O2"
    LINK_FLAGS="${omp}";;
    
  *)echo "Unsupported optimization level."
    exit;;
esac
