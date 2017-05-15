update_src()
{

#---------------------------- Step 1 ---------------------------------#
# Create default bckup copy *.f and *.inc files in model/subdir       #
# directory, if they do not exist.                                    #
#---------------------------------------------------------------------#
  subdir=$1

  run=${RUN_DIR}/$subdir
  src=${MFIX_SRC}/$subdir

# Loop over file extenions (.f, .F, .inc)
  for ext in $(echo ${EXT_LIST}); do

# Move into model sub directories located in the run directory.
# Note: The run directory (./) is considered the base model directory.
    if [ -d "${run}" ]; then
# Move into the run (sub)directory.
      cd "${run}"
# Make a list of all the files with the (ext) file extension. This list
# is empty if there are no files.
      list=
      if (eval ls *.${ext}) > /dev/null 2>&1; then :
        list=$(eval ls *.${ext})
      fi
# Loop over files in the run (sub)directory.
      for file in $(echo ${list}); do
# Check that the file exists in source directory.
        if test -r ${src}/${file}; then
# Construct a filename with a zero before the extension.
# (e.g., file.inc --> file.0inc)
          backup=`echo ${file} | sed 's/\./\.0/'`
# Copy the original file over to the backup extension.
          if test ! -r ${src}/$backup; then
            echo "  > Creating a backup: ${subdir}/${backup}"
            /bin/cp -f ${src}/$file ${src}/$backup
          fi
        fi
      done
    fi
  done



#---------------------------- Step 2 ----------------------------------#
# After generated the backup of the original file, copy the new user   #
# file from the run directory to the model directory.                  #
#----------------------------------------------------------------------#

# Move into the source (sub)directory.
  cd ${src}

# Make a list of existing backup files.
  list=
  if (eval ls *.0*) > /dev/null 2>&1; then
    list=$(eval ls *.0*)
  fi

  for backup in ${list};  do

# Verify that the user has a corresponding file.
    file=`echo ${backup} | sed 's/\.0/\./'`
    if test -r "${run}/${file}"; then

# Check if the file in the run directory differs from the file in the
# source (sub)direcotry.

# Check if the files have changed.
      cmp -s "${run}/${file}" "${file}"

# The file in the run (sub)directory is different than the one in the
# source directory. Copy the modified file into the source tree.
      if test $? = 1; then
        echo "  > Updating: ${subdir}/${file}"
        chmod u+w "${file}"
        /bin/cp -f "${run}/${file}" "${file}"
# The files are the same so only copy in the file if FORCED_COMPILE.
      elif test ${FORCE_COMPILE} = 1; then
        echo "  > Forced updated: ${subdir}/${file}"
        chmod u+w "${file}"
        /bin/cp -f "${run}/${file}" "${file}"
      fi

# User does not have a corresponding file. This assumption is that the
# backup file should be restored.
    else

# Compare the backup file and current source file.
      cmp -s "${backup}" "${file}"
# If the files are different, restore the source file with the backup.
      if test $? = 1; then
        /bin/mv -f "${backup}" "${file}"
        touch "${file}"
      else
# Remove the backup copy.
        /bin/rm -f "${backup}"
      fi
      echo "  > Restored: ${subdir}/${file}"
    fi
  done
}    # eof function update_dir


cd "${RUN_DIR}"
echo "Checking for user modified files."


for folder in $(echo ${FOLDER_LIST}); do
  update_src ${folder}
done

cd "${MFIX_SRC}"
