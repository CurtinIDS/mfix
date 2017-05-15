# Initialize some local variables
omp=
incs=
mkl=

mkl_libs=
mpi_libs=
misc_libs=

echo "Portland Group Fortarn Compiler"

MODDIRPREFIX="-module "

# Add some additinal flags to the object directory
DPO=${DPO_BASE}/${DPO}_PGI/
if test ! -d ${DPO}; then mkdir ${DPO}; fi

# Add DPO to the includes list.
incs=${incs}" -I${DPO}"

# Set OpenMP flags.
if test ${USE_SMP} = 1; then omp="-mp"; fi


# Set MPI flags.
if test ${USE_DMP} = 1; then
  FORTRAN_CMD=mpif90
  LINK_CMD=mpif90
else
  FORTRAN_CMD=pgfortran
  LINK_CMD=pgfortran
fi


# --> Verify compiler is in $PATH <-- #


# Set the MPI include path.  This must occur after
# the Fortran command is defined.
SET_MPI_INCLUDE
if test ${USE_DMP} = 1; then
  incs=${incs}" -I${MPI_INCLUDE_PATH}"
fi


# Set library paths:
ode="${DPO}odepack.a"
blas="${DPO}blas90.a"
dgtsv="${DPO}dgtsv90.a"
mpi_libs=
misc_libs=

mkl_libs="${blas} ${dgtsv}"
LIB_FLAGS="${ode} ${mkl_libs} ${misc_libs} ${mpi_libs} "

# Base flags for GNU Fortran
common="-c -I. ${incs} -Mnosave -byteswapio -cpp"

case $OPT in
  0)echo "Setting flags for debugging."
    dbg="-Mbounds -Mchkptr -Mchkfpstk -Mchkstk -Ktrap=divz,fp,ovf,unf"
    FORT_FLAGS="${omp} ${common} -Mfreeform ${dbg} -O0 -g "
    FORT_FLAGS3="${common} ${dbg} ${incs} -O0 -g "
    LINK_FLAGS="${omp} -g";;

  1)echo "Setting flags for low optimization."
    FORT_FLAGS="${omp} ${common} -Mfreeform -O1 -g "
    FORT_FLAGS3="${common} ${incs} -O1 -g "
    LINK_FLAGS="${omp} -g";;

  2)echo "Setting flags for medium optimization."
    FORT_FLAGS="${omp} ${common} -Mfreeform -O2 "
    FORT_FLAGS3="${common} -O1 "
    LINK_FLAGS="${omp}";;

  3)echo "Setting flags for high optimization."
    FORT_FLAGS="${omp} ${common} -Mfreeform -O4"
    FORT_FLAGS3="${common} -O2"
    LINK_FLAGS="${omp}";;

  *)echo "Unsupported optimization level."
    exit;;
esac
