echo "Building Makefile."

# Change the name to the full path.
MAKEFILE=${MFIX_SRC}/${MAKEFILE}
tmpMFILE=${MFIX_SRC}/tmp.make

# Remove the previous tmp file if it exists.
if test -f ${tmpMFILE}; then rm ${tmpMFILE}; fi

# Include any TAU definitions.
if test ${USE_TAU} = 1; then
  echo "include ${TAU_MAKEFILE}" >> ${tmpMFILE}
  echo "TAU_OPTIONS=${TAU_OPTIONS}" >> ${tmpMFILE}
fi

FORT_FLAGS="${FORT_FLAGS} -I${RUN_DIR}"
FORT_FLAGS="${FORT_FLAGS} -I${MFIX_SRC}/include"

# Include any NetCDF definitions.
if test ${USE_NETCDF} = 1; then
  FORT_FLAGS="${FORT_FLAGS} -I${NETCDF_INCLUDE} -DNETCDF"
  LIB_FLAGS="${LIB_FLAGS} -L${NETCDF_LIB} -lnetcdff"
fi

# Include any Trilinos definitions.
#if test ${USE_TRILINOS} = 1; then
#  echo "using Trilinos ?"
#fi

# Include the base definitions:
echo "AR=${AR}" >> ${tmpMFILE}
echo "DPO=${DPO}" >> ${tmpMFILE}
echo "OBJ_EXT=${OBJ_EXT}" >> ${tmpMFILE}
echo "FORTRAN_EXT=${FORTRAN_EXT}" >> ${tmpMFILE}
echo "LINK_FLAGS=${LINK_FLAGS}" >> ${tmpMFILE}
echo "LIB_FLAGS=${LIB_FLAGS}" >> ${tmpMFILE}
echo "EXEC_FILE=${EXEC_FILE}" >> ${tmpMFILE}
echo "MODDIRPREFIX=${MODDIRPREFIX}" >> ${tmpMFILE}

# Include any TAU definitions.
if test ${USE_TAU} = 1; then
  echo "FORT_FLAGS=\$(TAU_OPTIONS) ${FORT_FLAGS}" >> ${tmpMFILE}
  echo "FORT_FLAGS3=\$(TAU_OPTIONS) ${FORT_FLAGS3}" >> ${tmpMFILE}
  echo "FORTRAN_CMD=\$(TAU_COMPILER) ${FORTRAN_CMD}" >> ${tmpMFILE}
  echo "LINK_CMD=\$(TAU_COMPILER) ${LINK_CMD}" >> ${tmpMFILE}
else
  echo "FORT_FLAGS=${FORT_FLAGS}" >> ${tmpMFILE}
  echo "FORT_FLAGS3=${FORT_FLAGS3}" >> ${tmpMFILE}
  echo "FORTRAN_CMD=${FORTRAN_CMD}" >> ${tmpMFILE}
  echo "LINK_CMD=${LINK_CMD}" >> ${tmpMFILE}
fi

echo "Updating file list."
# Build a list of all files under the model directory.
# All sub folders with "donothing" in the name are skipped"
mlist=${MFIX_SRC}/file.list
if test -f "${mlist}"; then rm ${mlist}; fi

for dir in $(find . -type d); do
  tdir=$(echo $dir | grep 'donothing')
  if test -z "${tdir}"; then
    cd "${dir}"
#    echo "Adding files from $dir"
    list=
    if (eval ls *.f) > /dev/null 2>&1; then :
      list=$(eval ls *.f)
    fi

    for file in $(echo ${list}); do
      echo "${dir}/${file}" >> ${mlist}
    done
    if test "${dir}" != "."; then
#      echo "leaving $(pwd)"
      cd ..
    fi
#  else
#    echo "skipping ${dir}"
  fi
done

echo "Building make utility."
# build the MFIX make executable.
$(${FORTRAN_CMD} -o ${MFIX_TOOLS}/mms-auto.exe ${MFIX_TOOLS}/mms.f90)
if test $? -ne 0; then
  echo "Error building make utility. Aborting."
  exit
elif test ! -e ${MFIX_TOOLS}/mms-auto.exe; then
  echo "Error locating make utility. Aborting."
  exit
fi

echo "Running make utility."
$(${MFIX_TOOLS}/mms-auto.exe)
if test $? -ne 0; then
  echo "Error reported in make utility. Aborting."
  exit
fi

# Remove the file list.
rm ${mlist}

# Only copy over the tmp make file if it differs from
# the existing make file.
if test -f ${MAKEFILE}; then
  cmp -s ${tmpMFILE} ${MAKEFILE}
  if test ! $? = 0; then
    /bin/cp -f ${tmpMFILE} ${MAKEFILE}
#    /bin/rm ${tmpMFILE}
    echo "Makefile was updated."
  else
    echo "Makefile is up to date."
  fi
else
  /bin/cp -f ${tmpMFILE} ${MAKEFILE}
#  /bin/rm ${tmpMFILE}
  echo "Makefile was created."
fi

# Ensure that the object directory has write permission.
if test ! -w ${DPO}; then
  echo "Write permission added to object directory"
  chmod -R u+w ${DPO}
else
  echo "Object directory has write permission"
fi
