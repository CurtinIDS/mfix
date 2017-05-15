# Collect user build options:

if test $REQ_MODE = 1; then

  echo
  echo "=============================================================="
  echo " Execution modes:"
  echo "=============================================================="
  echo "[1] Serial"
  echo "[2] Parallel, Shared Memory (SMP)"
  echo "[3] Parallel, Distributed Memory (DMP)"
  echo "[4] Parallel, Hybrid (SMP+DMP)"
  echo " "
  echo -n "Select the mode of execution [1] : "
  read mode_of_execution

  case $mode_of_execution in
    1) USE_SMP=0; USE_DMP=0;;
    2) USE_SMP=1; USE_DMP=0;;
    3) USE_SMP=0; USE_DMP=1;;
    4) USE_SMP=1; USE_DMP=1;;
    *) USE_SMP=0; USE_DMP=0;;
  esac
fi




if test ${REQ_OPT} = 1; then
  echo ""
  echo "=============================================================="
  echo " Compiler optimization levels:"
  echo "=============================================================="
  echo "[0] Level 0: No optimization - Debug mode"
  echo "[1] Level 1: O1 optimization - Low"
  echo "[2] Level 2: O2 optimization - Moderate"
  echo "[3] Level 3: O3 optimization - Aggressive"
  echo " "
  echo -n "Select the level of optimization [3] : "
  read level_of_optimization

  case $level_of_optimization in
    0)OPT=0;;
    1)OPT=1;;
    2)OPT=2;;
    3)OPT=3;;
    *)OPT=3;;
  esac
fi

if test ${EXPERT} = 1; then
  echo ""
  echo "=============================================================="
  echo " Option to re-compile source files in run directory:"
  echo "=============================================================="
  echo "[1] Do not force re-compilation"
  echo "[2] Force re-compilation"
  echo " "
  echo -n "Select Option to re-compile source files in run directory [1] : "
  read re_compile_option

  case $re_compile_option in
    1) FORCE_COMPILE=0;;
    2) FORCE_COMPILE=1;;
    *) FORCE_COMPILE=0;;
  esac
fi


# Add mode flags and optimization to object directory name.
if test $USE_DMP = 1; then DPO=${DPO}_DMP; fi
if test $USE_SMP = 1; then DPO=${DPO}_SMP; fi
DPO=${DPO}_OPT${OPT}
