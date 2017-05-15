#######################################################################
# Function: set_MPI_INCLUDE                                           #
#                                                                     #
# Try to determine the location of the mpif.h file.                   #
#                                                                     #
#######################################################################
SET_MPI_INCLUDE(){

# Serial runs uses the fake mpi include file.
  if test ${USE_DMP} = 0; then
    MPI_INCLUDE_PATH=
  else
    echo "Searching for mpif.h file"

# This may have been provided via input arguments. Thus, we have to
# check that MPI_INCLUDE_PATH isn't set for all checks.

    if test -z ${MPI_INCLUDE_PATH}; then
# User supplied mpi path but not the include path. Maybe the include
# path is a subdirectory of the mpi directory?
      if test ! -z ${MPI_PATH}; then
        inc=${MPI_PATH}/include
        if test -f ${inc}; then  MPI_INCLUDE_PATH=${inc}; fi
      fi
    fi

# We -might- be able to get it with showme: (Open MPI)
    if test -z ${MPI_INCLUDE_PATH}; then
      if (eval ${FORTRAN_CMD} -showme:incdirs) > /dev/null 2>&1; then
        MPI_INCLUDE_PATH=$(eval ${FORTRAN_CMD} -showme:incdirs)
      fi
    fi

# We -might- be able to get it with show (Intel MPI)
    if test -z ${MPI_INCLUDE_PATH}; then
      if (eval ${FORTRAN_CMD} -show) > /dev/null 2>&1; then
        out=$(eval ${FORTRAN_CMD} -show)
        for arg in ${out}; do
          dir=$(echo ${arg} | grep include | sed -e 's/-I\(.*\)\/include/\1/')
          if test ! -z ${dir}; then
            inc=${dir}/include
            if test -f ${inc}/mpif.h; then MPI_INCLUDE_PATH=${inc}; fi
          fi
        done
      fi
    fi

# We -might- be able to get it with compile-info (MPICH)
    if test -z ${MPI_INCLUDE_PATH}; then
      if (eval ${FORTRAN_CMD} -compile-info) > /dev/null 2>&1; then
        out=$(eval ${FORTRAN_CMD} -compile-info)
        for arg in ${out}; do
          dir=$(echo ${arg} | grep include | sed -e 's/-I\(.*\)\/include/\1/')
          if test ! -z ${dir}; then
            inc=${dir}/include
            if test -f ${inc}/mpif.h; then MPI_INCLUDE_PATH=${inc}; fi
          fi
        done
      fi
    fi


# I've checked all the placed I know where to look. Ask the user.
    if test -z ${MPI_INCLUDE_PATH}; then
      echo "Unable to locate mpif.h"
      echo -n "Please provide the location of the mpif.h file: "
      read MPI_INCLUDE_PATH
    else
      echo "Using mpif.h: ${MPI_INCLUDE_PATH}"
      export MPI_INCLUDE_PATH=${MPI_INCLUDE_PATH}
    fi

    if test -f ${MPI_INCLUDE_PATH}/mpif.h; then
      if test -f ${MFIX_SRC}/mpif.h; then
        /bin/rm ${MFIX_SRC}/mpif.h
        if test ! $? = 0; then
          echo "Unable to remove existing mpif.h file in source directory!"
          echo "Aborting."
          exit
        fi
      fi
      ln -sf ${MPI_INCLUDE_PATH}/mpif.h ${MFIX_SRC}/dmp_modules/mpif.h
    else
      echo "Fatal Error: Unable to locate the mpif.h header file!"
      echo "Aborting."
      exit
    fi
  fi
} # END set_MPI_INCLUDE
