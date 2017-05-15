# Ensure that we are in the post directory.
echo "Copying model directory incs"

/bin/cp -f ${MFIX_SRC}/namelist.inc ${MFIX_POST}/.

if test ! -d ${MFIX_POST}/des; then
  echo "Copying des subfolder incs"
  mkdir ${MFIX_POST}/des
fi
/bin/cp -f ${MFIX_SRC}/des/desnamelist.inc ${MFIX_POST}/des/desnamelist.inc


if test ! -d ${MFIX_POST}/qmomk; then
  echo "Copying qmomk subfolder incs"
  mkdir ${MFIX_POST}/qmomk
fi
/bin/cp -f ${MFIX_SRC}/qmomk/qmomknamelist.inc ${MFIX_POST}/qmomk/qmomknamelist.inc

