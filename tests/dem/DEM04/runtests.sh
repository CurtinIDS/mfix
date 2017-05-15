#!/bin/bash -lex

RUN_NAME="DEM04"

DES_IM=ADAMS_BASHFORTH
for DES_MEW in 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0; do
  rm -f ${RUN_NAME}* &> /dev/null
  time -p ./mfix DES_INTG_METHOD=\"${DES_IM}\" \
    MEW=${DES_MEW} MEW_W=${DES_MEW}
done

#diff -q POST_posvel.dat AUTOTEST/POST_posvel.dat

numdiff \
    -a 0.000001 -r 0.05 \
    --exclude=1:4 --exclude=2:4 \
    AUTOTEST/POST_AVEL.dat POST_AVEL.dat

numdiff \
    -a 0.000001 -r 0.05 \
    --exclude=1:4 --exclude=2:4 \
    AUTOTEST/POST_TVEL.dat POST_TVEL.dat

numdiff \
    -a 0.000001 -r 0.05 \
    --exclude=1:4 --exclude=2:4 \
    AUTOTEST/POST_TIME.dat POST_TIME.dat
