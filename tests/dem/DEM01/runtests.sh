#!/bin/bash -lex

RUN_NAME="DEM01"

DES_IM=EULER
for DES_IM in EULER ADAMS_BASHFORTH; do
  DES_KN=10000
  for DES_ETA in 0.9 0.8 0.7 0.6; do
    rm -f ${RUN_NAME}* &> /dev/null
    time -p ./mfix DES_INTG_METHOD=\"${DES_IM}\" \
      DES_EN_INPUT=${DES_ETA} DES_EN_WALL_INPUT=${DES_ETA} \
      KN=${DES_KN} KN_W=${DES_KN}
  done

  for DES_KN in 25000 50000 100000; do
    for DES_ETA in 1.0 0.9 0.8 0.7 0.6; do
      rm -f ${RUN_NAME}* &> /dev/null
      time -p ./mfix DES_INTG_METHOD=\"${DES_IM}\" \
        DES_EN_INPUT=${DES_ETA} DES_EN_WALL_INPUT=${DES_ETA} \
        KN=${DES_KN} KN_W=${DES_KN}
    done
  done
done
#diff -q POST_posvel.dat AUTOTEST/POST_posvel.dat

numdiff \
    -a 0.000001 -r 0.05 \
    --exclude=1:5 --exclude=2:5 \
    AUTOTEST/POST_POS.dat POST_POS.dat

numdiff \
    -a 0.000001 -r 0.05 \
    --exclude=1:5-6 --exclude=2:5-6 \
    AUTOTEST/POST_VEL.dat POST_VEL.dat
