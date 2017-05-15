# Ensure that we are in the run directory.
cd "${RUN_DIR}"

# MFIX_SRC is relative to RUN_DIR.  Make it absolute.
cd "${MFIX_SRC}"
MFIX_SRC=$(pwd)

# Up to this point, MFIX_CONFIG is relative to MFIX_SRC
cd "${MFIX_SRC}/../config"
MFIX_CONFIG=$(pwd)

# Set the path to the tools directory
cd "${MFIX_SRC}/../tools"
MFIX_TOOLS=$(pwd)

# Set the path to the tools directory
cd "${MFIX_SRC}/../post_mfix"
MFIX_POST=$(pwd)

# Grab some system information.
opsys=$(uname -s)
proctyp=$(uname -p)
if test $proctyp = "unknown"; then
  proctyp=$(uname -m)
fi

cd "${RUN_DIR}"
# Only execute if the run directory contains usr_rates.f
if test -f "usr_rates.f" || test -f "des/usr_rates_des.f"; then
#  if test -f "usr_rates.f"; then echo -n "usr_rates.f"; fi
#  if test -f "des/usr_rates_des.f"; then echo -n "usr_rates_des.f"; fi
  REQ_RXNS=1
else
  REQ_RXNS=0
fi
# Get back into the model directory.
cd "${MFIX_SRC}"
