if test ${REQ_COMP} = 1; then

  case "$opsys" in
    *Linux*)  mfix_os="LINUX";;
    *CYGWIN*) mfix_os="WINDOWS";;
    *MING*)   mfix_os="WINDOWS";;
    *Darwin*) mfix_os="MacOS";;    
    *) echo "Error: Unsupported OS. $opsys"
    exit;;
  esac

# Differentiate MFIX and POSTMFIX compilations
  case "$EXEC_FILE" in
    mfix.exe)  exec_label="MFIX";;
    post_mfix) exec_label="POSTMFIX";;
    *) echo "Error: Unknown executable file: $EXEC_FILE"
    exit;;
  esac

  case $mfix_os in
    LINUX)
      echo
      echo
      echo "=============================================================="
      echo " $exec_label Compilation directives for following compilers:"
      echo "=============================================================="
      echo "  [1] GCC (gfortran) version 4.4 and above "
      echo "  [2] Portland Group (pgf90) version 11.7 and above"
      echo "  [3] Intel (ifort) version 11.1 and above "

      if test $EXPERT = 1; then
        echo
        echo " <--- @NETL  -------------------------------------->"
        echo " [10] SBEUC :: Intel (ifort) version 14.0"
        echo " [11] SBEUC :: GCC (gfrtran) version 4.9"
        echo
        echo " <--- @OLCF  -------------------------------------->"
        echo " [20] Cray XT4 :: Portland Group (pgf90)"
        echo
        echo " <--- @ALCF  -------------------------------------->"
        echo " [30] Blue Gene/P :: IBM (xlf90)"
        echo " [31] Blue Gene/Q :: IBM (xlf90)"
        echo
        echo " <--- @NERSC  ------------------------------------->"
        echo " [40] Hopper :: CRAY (ftn)"
        echo " [41] Hopper :: PGI (ftn)"
        echo " [42] Hopper :: Intel (ftn)"
        echo " [43] Hopper :: GCC (ftn)"
        echo
      fi
      echo " "
      echo -n "Select the compiler to compile $exec_label? [1] "
      read compiler

      case $compiler in
# Generic Flags
        1)COMP_FILE="gcc_default.sh";;
        2)COMP_FILE="portland_default.sh";;
        3)COMP_FILE="intel_default.sh";;
# NETL Systems
       10)COMP_FILE="sbeuc_intel-131.sh";;
       11)COMP_FILE="sbeuc_gcc-46.sh";;
# OLCF Systems
       20)COMP_FILE="olcf_xt4_pg.sh";;
# ALCF Systems
       30)COMP_FILE="alcf_bgp_ibm.sh";;
       31)COMP_FILE="alcf_bgq_ibm.sh";;
# NERSC Systems
       40)COMP_FILE="nersc_hopper_cray.sh";;
       41)COMP_FILE="nersc_hopper_pgi.sh";;
       42)COMP_FILE="nersc_hopper_intel.sh";;
       43)COMP_FILE="nersc_hopper_gcc.sh";;

        *)COMP_FILE="gcc_default.sh";;
      esac

# Change from file name to absolute reference.
      COMP_FILE=${MFIX_CONFIG}/compilers/${COMP_FILE}
    ;;

#####################  WINDOWS   ##########################
    WINDOWS)
      echo
      echo
      echo " $opsys on Windows Detected"
      echo
      echo "=============================================================="
      echo "$exec_label Compilation directives available for following compilers:"
      echo "=============================================================="
      echo " [1] g95"
      echo " [2] gfortran"
      echo " "
      echo -n "Select the compiler to compile $exec_label? [2] "
      read compiler

      case $compiler in
        1)COMP_FILE=win_g95.sh;;
        2)COMP_FILE="win_gnu.sh";;
        *)COMP_FILE="win_gnu.sh";;
      esac
# Change from file name to absolute reference.
      COMP_FILE=${MFIX_CONFIG}/compilers/${COMP_FILE}
    ;;

#####################  Apple MacOS   ##########################
    MacOS)
      echo
      echo
      echo " $opsys on Apple MacOS Detected"
      echo
      echo "=============================================================="
      echo "MFIX Compilation directives available for following compilers:"
      echo "=============================================================="
      echo " [1] GNU Fortran (gfortran)"
      echo " [2] Portland Group (pgf90)"  
      echo " [3] Intel (ifort)"               
      echo " "
      echo -n "Select the compiler to compile MFIX? [1] "
      read compiler

      case $compiler in
        1)COMP_FILE="gnu_mac_os.sh";;
        2)COMP_FILE="portland_mac_os.sh";;
        3)COMP_FILE="intel_mac_os.sh";;      
        *)COMP_FILE="gnu_mac_os.sh";;
      esac

#ALPEMI MacOS system type written out so rxn_preproc.sh handles special sed exception  
      export mfix_os="MacOS"
#      echo "Machine 0" ${mfix_os}
      MACHINE_TYPE_SCRIPT=${MFIX_CONFIG}/compilers/Load_MACHINE_TYPE.sh
      echo "#!/bin/bash" > $MACHINE_TYPE_SCRIPT
      echo "export mfix_os=MacOS" >> $MACHINE_TYPE_SCRIPT
      chmod ug+rwx $MACHINE_TYPE_SCRIPT

# Change from file name to absolute reference.
      COMP_FILE=${MFIX_CONFIG}/compilers/${COMP_FILE}      
  esac
fi
