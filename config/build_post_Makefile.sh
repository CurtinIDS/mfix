echo "Building Makefile."

cd ${MFIX_POST}

# Change the name to the full path.
MAKEFILE=${MFIX_POST}/${MAKEFILE}
tmpMFILE=${MFIX_POST}/tmp.make

# Remove the previous tmp file if it exists.
if test -f ${tmpMFILE}; then rm ${tmpMFILE}; fi

FORT_FLAGS="${FORT_FLAGS} -I./include"
FORT_FLAGS="${FORT_FLAGS} -I${MFIX_SRC}/include"

# Include any NetCDF definitions.
if test ${USE_NETCDF} = 1; then
  FORT_FLAGS="${FORT_FLAGS} -I${NETCDF_INCLUDE}"
  LIB_FLAGS="${LIB_FLAGS} -L${NETCDF_LIB} -lnetcdff"
fi

# Include the base definitions:
echo "DPO=${DPO}" >> ${tmpMFILE}
echo "OBJ_EXT=${OBJ_EXT}" >> ${tmpMFILE}
echo "FORTRAN_EXT=${FORTRAN_EXT}" >> ${tmpMFILE}
echo "FORT_FLAGS=${FORT_FLAGS}" >> ${tmpMFILE}
echo "FORT_FLAGS3=${FORT_FLAGS3}" >> ${tmpMFILE}
echo "FORTRAN_CMD=${FORTRAN_CMD}" >> ${tmpMFILE}
echo "LINK_FLAGS=${LINK_FLAGS}" >> ${tmpMFILE}
echo "LINK_CMD=${LINK_CMD}" >> ${tmpMFILE}
echo "LIB_FLAGS=" >> ${tmpMFILE}
echo "EXEC_FILE=${EXEC_FILE}" >> ${tmpMFILE}
echo "MODDIRPREFIX=${MODDIRPREFIX}" >> ${tmpMFILE}


echo "Building make utility."
# build the MFIX make executable.
$(${FORTRAN_CMD} -o ${MFIX_TOOLS}/mms_post-auto.exe ${MFIX_TOOLS}/mms_post.f90)
if test $? -ne 0; then
  echo "Error building make utility. Aborting."
  exit
elif test ! -e ${MFIX_TOOLS}/mms_post-auto.exe; then
  echo "Error locating make utility. Aborting."
  exit
fi

echo "Running make utility."
$(${MFIX_TOOLS}/mms_post-auto.exe)
if test $? -ne 0; then
  echo "Error reported in make utility. Aborting."
  exit
fi

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
