omp=
incs=

mpi_libs=
misc_libs=

echo "SBEUC@NETL :: GCC Fortran Compiler"

MODDIRPREFIX="-J"

# Add some additinal flags to the object directory
DPO=${DPO_BASE}/${DPO}_GNU_SBEUC/
if test ! -d $DPO; then  mkdir $DPO; fi

# Add DPO to the includes list.
incs="${incs} -I${DPO}"

# Set OpenMP flags.
if test $USE_SMP = 1; then omp="-fopenmp"; fi

# Set MPI flags.
if test $USE_DMP = 1; then
  FORTRAN_CMD=mpif90
  LINK_CMD=mpif90
else
  FORTRAN_CMD=gfortran
  LINK_CMD=gfortran
fi

# Set the MPI include path.  This must occur after
# the Fortran command is defined.
SET_MPI_INCLUDE
if test ${USE_DMP} = 1; then
  incs="${incs} -I${MPI_INCLUDE_PATH}"
fi


# --> Verify compiler is in $PATH <-- #


# Set library paths:
ode="${DPO}odepack.a"
blas="${DPO}blas90.a"
dgtsv="${DPO}dgtsv90.a"

mkl_libs="${blas} ${dgtsv}"
LIB_FLAGS="${ode} ${mkl_libs} ${misc_libs} ${mpi_libs} "

# Debug flags for GNU Fortran
dbg=
if test ${USE_DEBUG} = 1; then dbg="-g"; fi

# Base flags for GNU Fortran:
common="-c -I. ${incs} -fconvert='big-endian' -cpp"
# GCC Fortran syntax flags.
sflag="-ffree-form -ffree-line-length-0"
# GCC hardware-specific flags.
hdwf="-mtune=corei7-avx -march=corei7-avx -masm=intel"

# Set debugging and optimization flags.
case $OPT in
  0)echo "Setting flags for debugging."
    dbg="-fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow"
    FORT_FLAGS="${omp} ${common} ${sflag} ${dbg} -g -O0"
    FORT_FLAGS3="${common} -g -O0"
    LINK_FLAGS="${omp} -g";;

  1)echo "Setting flags for low optimization."
    FORT_FLAGS="${omp} ${common} ${sflag} -O2 ${dbg}"
    FORT_FLAGS3="${common} -g -O1 ${dbg}"
    LINK_FLAGS="${omp} ${dbg}";;

  2)echo "Setting flags for medium optimization."
    FORT_FLAGS="${omp} ${common} ${sflag} -g -O2 ${dbg}"
    FORT_FLAGS3="${common} -O1 ${dbg}"
    LINK_FLAGS="${omp} ${dbg}";;

  3)echo "Setting flags for high optimization."
    FORT_FLAGS="${omp} ${common} ${sflag} ${hdwf} -Ofast ${dbg}"
    FORT_FLAGS3="${common} -O2 ${dbg}"
    LINK_FLAGS="${omp} ${dbg}";;

  *)echo "Unsupported optimization level."
    exit;;
esac
