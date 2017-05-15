MFIX_NETCDF=${MFIX_SRC}/netcdf
netcdf_msg=

if test ${USE_NETCDF} = 1; then
  echo "NetCDF is enabled."
  netcdf_msg="Updated"
else
  MFIX_NETCDF=${MFIX_NETCDF}/donothing
  netcdf_msg="Reverted"
fi

list=
if (eval ls *.fi) > /dev/null 2>&1; then :
  list=$(eval ls *.fi)
fi

for file in ${list}; do
  cmp -s ${MFIX_SRC}/${file} ${MFIX_NETCDF}/${file}
  if test $? = 1; then
    echo "${netcdf_msg} MFIX NetCDF file: ${file}"
    cp ${MFIX_NETCDF}/${file} ${MFIX_SRC}/.
  fi
done
