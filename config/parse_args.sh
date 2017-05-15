#!/bin/bsh -f

show_version() {
echo
echo "======================================================================"
echo ""
echo "            MFIX: Multiphase Flow with Interphase eXchanges"
echo ""
echo "                            mfix.netl.doe.gov"
echo ""
echo "                            Version ${MFIX_VERSION}"
echo ""
echo "======================================================================"
echo
echo "MFIX source: ${MFIX_SRC}"
echo "Operating System: ${opsys}"
echo "Processor type: ${proctyp}"
echo

}


show_usage() {
echo ""
echo "======================================================================"
echo ""
echo "Usage: mfix/model/make_mfix [OPTIONS]... [VAR=VALUE]..."
echo ""
echo "Invoke make_mfix without arguments for a guided build."
echo ""
echo "Configuration:"
echo ""
echo " -h, --help           Display this help and exit"
echo " -V, --version        Display version information and exit"
echo " -l, --long           Display all build options"
echo " -c, --clean          Remove previous build directory"
echo " -d, --default        Enable default build (serial, GCC compiler)"
echo ""
echo "Optional Features:"
echo ""
echo " --force-recompile    Forced recompile of modified files"
echo " --serial             Disable shared and distributed memory parallel"
echo " --smp                Enable shared memory parallel (OpenMP)"
echo " --dmp                Enable distributed memory paralle (MPI)"
echo " --debug              Compile with debugging information (e.g., -g flag)"
echo ""
echo " --exe=name - Specify name of the executable (default: mfix.exe)"
echo ""
echo " --opt=level - Specify compiler optimization level"
echo "     O0     -   No optimization with debugging flags"
echo "     O1     -   O1 optimization - Low"
echo "     O2     -   O2 optimization - Moderate"
echo "     O3     -   O3 optimization - Aggressive"
echo ""
echo " --compiler=option - Specify compiler flag file. There are three"
echo "                     generic options that cover most user needs:"
echo "    gcc            - GCC compiler (gfortran)"
echo "    intel          - Intel compiler (ifort)"
echo "    portland       - Portland Groupe (pgf90)"
echo ""
echo " Specific hardware configuration files: This script will first check for"
echo " the specified file in the run directory, then within mfix/config directory."
echo ""
echo "    gcc_default.sh         - same as gcc"
echo "    intel_default.sh       - same as intel"
echo "    portland_default.sh    - same as portland"
echo ""
echo "    Apple MacOS system files:"
echo "    gnu_mac_os.sh         - same as gcc"
echo "    intel_mac_os.sh       - same as intel"
echo "    portland_mac_os.sh    - same as portland"
echo ""
echo "    SBEUC system files:"
echo "      sbeuc_gcc-46.sh      - GCC 4.6 and Open MPI 1.5.5"
echo "      sbeuc_intel-131.sh   - Intel 13.1 and Intel MPI"
#echo ""
#echo "    OLCF system files:"
#echo "      olcf_xt4_pg.sh       - Cray XT4 - Portland Group (pgf90)"
#echo ""
#echo "    ALCF system files:"
#echo "      alcf_bgp_ibm.sh      - IBM BG/P (xlf90)"
#echo "      alcf_bgq_ibm.sh      - IBM BG/Q (xlf90)"
echo ""
echo "    Hopper@NERSC - Cray XE6:"
echo "      nersc_hopper_cray.sh - Cray compiler (ftn)"
echo "      nersc_hopper_intel.sh- Intel compiler (ftn)"
echo "      nersc_hopper_pgi.sh  - PGI compiler (ftn)"
echo "      nersc_hopper_gnu.sh  - GCC compiler (ftn)"
echo ""
echo ""
echo "Advanced Options:"
echo ""
echo " --codecov            Enable code coverage utility (Intel compiler)"
echo " --mkl                Enable Intel Math Kernel Library (requires building with Intel Fortran compiler)"
echo " --mic                Enable Intel MIC flags"
echo " --skip-rxn           Skip the reaction preprocessor"
echo ""
echo " -j                   Enable parallel build (GNU make utility only)"
echo ""
echo " --enable-tau         Enable TAU profiling support. Specify the location"
echo "                      of the TAU library by setting the TAUROOT"
echo "                      environment variable."
echo ""
echo " --enable-netcdf      Enable NetCDF support. Set NETCDF_INCLUDE to directory"
echo "                      of NetCDF include files and set NETCDF_LIB to directory of"
echo "                      NetCDF shared library"
echo ""
echo " --mpi=PATH           MPI instalation directory"
echo " --mpi_include=PATH   MPI include directory"
echo " --mpi_lib=PATH       MPI library directory"
echo ""
echo " Some influential environment variables:"
echo ""
echo " NETCDF_HOME          Specifies the location of the NetCDF library."
echo " TAUROOT              Specifies the location of the TAU library"
#echo "  FORTRAN_CMD"
#echo "  FORT_FLAGS"
echo ""
echo "======================================================================"
echo ""
echo ""
exit $1
}

# Display version.
#-------------------------------------------------------------------------->>
show_version

# If version information is explictly requested, exit.
for arg in $input; do
  case ${arg} in
    "--version"|"-V") exit 0;;
    *) echo "" > /dev/null 2>&1;
  esac
done



# Loop the the argument list.
for arg in $input; do

  case ${arg} in

# Display help.
#-------------------------------------------------------------------------->>
    "--help"|"-h")
      show_usage 0;;


# Clean out the last build.
#-------------------------------------------------------------------------->>
    "--clean" | "-c") CLEAN_OBJS=1;;

# Specify default compile options.
#-------------------------------------------------------------------------->>
    "--default"|"-d" )

      EXPERT=0
      FORCE_COMPILE=0
      COMP_FILE="${MFIX_CONFIG}/compilers/gcc_default.sh"
      REQ_COMP=0

      if test -z ${USE_SMP}; then USE_SMP=0; fi
      if test -z ${USE_DMP}; then USE_DMP=0; fi
      REQ_MODE=0

      OPT=3
      REQ_OPT=0;;


# Show all available options.
#-------------------------------------------------------------------------->>
    "--long" | "-l" )
      echo "All compilation options will be shown"
      EXPERT=1;;


# A recompile command is executed. Invoke the make file
# that was previously generated.
#-------------------------------------------------------------------------->>
    "--repeat" | "-r" )
      if test ! -f "${MFIX_SRC}/${MAKEFILE}"; then
        echo "Cannot locate Makefile in model directory."
        echo "Unable to repeat last compile."
        exit -1
      fi
      mfile=${MFIX_SRC}/${MAKEFILE}
      echo "Using last compile settings."
      DPO=$(grep "DPO=" ${mfile} | cut -d "=" -f2)
      AUTOCOMPILE=1;;

# Option for multiple make jobs
#-------------------------------------------------------------------------->>
    "-j" ) MAKE_ARGS="-j";;

# Enable specify debug flags.
#-------------------------------------------------------------------------->>
    "--debug" ) USE_DEBUG="1";;

# Code coverage flag
#-------------------------------------------------------------------------->>
    "--codecov" ) USE_CODECOV="1";;

# Specify optimization level.
#-------------------------------------------------------------------------->>
    "--opt="* )
      OPT=$(echo ${arg} | cut -d "=" -f2)
      case ${OPT} in
        O0) OPT=0; USE_DEBUG=1;;
        O1) OPT=1;;
        O2) OPT=2;;
        O3) OPT=3;;
        O4) OPT=4;;
        *) echo "Error unknown optimization level: ${arg}"
           echo "First character is upper case letter O, not zero"
           echo "Aborting."
           exit -1;;
      esac
      echo "Specified optimization level: ${OPT}"
      REQ_OPT=0;;


# Enable SMP model
#-------------------------------------------------------------------------->>
    "--smp" )
      USE_SMP=1
      echo "SMP build mode enabled."
      if test -z ${USE_DMP}; then USE_DMP=0; fi
      REQ_MODE=0;;


# Enable DMP model
#-------------------------------------------------------------------------->>
    "--dmp" )
      USE_DMP=1
      echo "DMP build mode enabled."
      if test -z ${USE_SMP}; then USE_SMP=0; fi
      REQ_MODE=0;;


# Enable Serial model
#-------------------------------------------------------------------------->>
    "--serial" )
      USE_SMP=0
      USE_DMP=0
      REQ_MODE=0;;

# Enable MIC model
#-------------------------------------------------------------------------->>
    "--mic" )
      USE_MIC=1
      REQ_MODE=0;;

# Bypass chemical reaction preprocessing.
#-------------------------------------------------------------------------->>
    "--skip-rxn" )
      if test ${REQ_RXNS} = 1; then
        echo "WARNING: Skipping chemical reaction preprocessing."
      fi
      REQ_RXNS=0;;


# Compiler selection
#-------------------------------------------------------------------------->>
    "--compiler="*)
      compiler=$(echo ${arg} | cut -d "=" -f2)
      case ${compiler} in
        gcc) COMP_FILE=${MFIX_CONFIG}/compilers/gcc_default.sh;;
        intel) COMP_FILE=${MFIX_CONFIG}/compilers/intel_default.sh;;
        portland) COMP_FILE=${MFIX_CONFIG}/compilers/portland_default.sh;;
        *) comp=$(echo ${arg} | cut -d "=" -f2)
          if test -f ${RUN_DIR}/${comp}; then
            COMP_FILE=${RUN_DIR}/${comp}
          elif test -f ${MFIX_CONFIG}/compilers/${comp}; then
            COMP_FILE=${MFIX_CONFIG}/compilers/${comp}
          elif test -f ${MFIX_CONFIG}/${comp}; then
            COMP_FILE=${MFIX_CONFIG}/${comp}
          elif test -f ${comp}; then
            COMP_FILE=${comp}
          else
            echo "  Error: Unable to locate compiler file!"
            echo "   >>> ${arg}"
            exit -1
          fi ;;
      esac
      REQ_COMP=0;;

# Enable MKL
#-------------------------------------------------------------------------->>
    "--mkl" )
      USE_MKL=1;;

# Enable force source file recompiling.
#-------------------------------------------------------------------------->>
    "--force-recompile"* )
      FORCE_COMPILE=1;;


# Use specified build directory.
#-------------------------------------------------------------------------->>
    "--build="*)
      dir=$(echo ${arg} | cut -d "=" -f2)
      if test ! -d ${dir}; then
        echo "  Specified build directory not found!"
        echo "   >>> ${dir}"
        exit -1
      fi
      cd ${dir}
      set `pwd` ; DPO_BASE=$1
      cd ${MFIX_SRC}
      echo "Build directory: ${DPO_BASE}";;


# Use specified executable name. (mfix.exe is the default)
#-------------------------------------------------------------------------->>
    "--exe="*)
      EXEC_FILE=$(echo ${arg} | cut -d "=" -f2)
      if test -z ${EXEC_FILE}; then
        echo "Specified executable name is empty!"
        echo "Aborting make_mfix."
        exit -1
      fi
      echo "User specified executable name: ${EXEC_FILE}";;


# MPI installation directory. This folder should contain the include
# and library directories.
#-------------------------------------------------------------------------->>
    "--mpi="*)
      dir=$(echo ${arg} | cut -d "=" -f2)
      if test ! -d ${dir}; then
        echo "   Specified MPI path not found!"
        echo "   >>> ${dir}"
        exit -1
      fi
      cd ${dir}
      set `pwd` ; MPI_PATH=$1
      cd ${MFIX_SRC}
      echo "MPI path: ${MPI_PATH}";;


# User specified path to MPI include directory.
#-------------------------------------------------------------------------->>
    "--mpi_include="*)
      dir=$(echo ${arg} | cut -d "=" -f2)
      if test ! -d ${dir}; then
        echo "   Specified MPI include path not found!"
        echo "   >>> ${dir}"
        exit -1
      fi
      cd ${dir}
      set `pwd` ; MPI_INCLUDE_PATH=$1
      if test ! -f "mpif.h"; then
        echo "  Specified mpi_include does not contain mpif.h"
        echo "   >>> $PWD"
        exit -1
      fi
      cd ${MFIX_SRC}
      echo "MPI include path: ${MPI_INCLUDE_PATH}";;


# MPI Library directory
#-------------------------------------------------------------------------->>
    "--mpi_lib="*)
      dir=$(echo ${arg} | cut -d "=" -f2)
      if test ! -d ${dir}; then
        echo "  Specified mpi library path not found!"
        echo "   >>> ${dir}"
        exit -1
      fi
      cd ${dir}
      set `pwd` ; MPI_LIB_PATH=$1
      cd ${MFIX_SRC}
      echo "MPI library path: ${MPI_LIB_PATH}";;


# Enable TAU profiling
#-------------------------------------------------------------------------->>
    "--enable-tau" )
      if test -z ${TAU_MAKEFILE}; then
        echo "  Fatal Error: TAU_MAKEFILE not set!"
        exit -1
      fi
      if test -z ${TAU_OPTIONS}; then
        echo "  Fatal Error: TAU_OPTIONS not set!"
        exit -1
      fi
      USE_TAU=1
      DPO=${DPO}_TAU;;


# Enable NetCDF output
#-------------------------------------------------------------------------->>
    "--enable-netcdf" )
      if test -z ${NETCDF_INCLUDE}; then
        echo "Error: Building with NetCDF support requires NETCDF_INCLUDE set to directory of NetCDF Fortran include files (netcdf.f90 or netcdf.mod)"
        exit -1
      fi
      if test -z ${NETCDF_LIB}; then
        echo "Error: Building with NetCDF support requires NETCDF_LIB set to directory of libnetcdff.so"
        exit -1
      fi
      USE_NETCDF=1
      DPO=${DPO}_NCDF;;


# An unknown build command was specified. Print the usage
# information and exit.
#-------------------------------------------------------------------------->>
    *)echo "Unknown flag: ${arg}"
      show_usage -1;;

  esac
done




#echo "Forcing a stop in parse_args.sh"
#exit
