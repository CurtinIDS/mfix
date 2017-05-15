#!/bin/bash -lex

MFIX=${MFIX_HOME-"../../../"}

rm -f POST_*.dat &> /dev/null

RUN_NAME="DEM05"

rm -f ${RUN_NAME}* &> /dev/null
time -p ./mfix DES_COLL_MODEL=\"LSD\" \
  KN=1.72d7 KT_FAC="@(1.48/1.72)" KN_W=1.72d7 KT_W_FAC="@(1.48/1.72)" \
  "DES_EN_INPUT(1:3)=3*1.0 DES_EN_WALL_INPUT(1:2)=2*1.0"

rm -f ${RUN_NAME}* &> /dev/null
time -p ./mfix DES_COLL_MODEL=\"HERTZIAN\" \
  "E_YOUNG(1)=380.0d9 E_YOUNG(2)=70.0d9 Ew_YOUNG=70.0d9" \
  "V_POISSON(1)=0.23  V_POISSON(2)=0.25 Vw_POISSON=0.25" \
  "DES_EN_INPUT(1:3)=3*1.0 DES_EN_WALL_INPUT(1:2)=2*1.0" \
  "DES_ET_INPUT(1:3)=3*1.0 DES_ET_WALL_INPUT(1:2)=2*1.0"

##diff -q POST_posvel.dat AUTOTEST/POST_posvel.dat

numdiff \
    -a 0.000001 -r 0.05 \
    AUTOTEST/POST_ALPHA.dat POST_ALPHA.dat || echo "file POST_ALPHA.dat does not match"

numdiff \
    -a 0.000001 -r 0.05 \
    AUTOTEST/POST_COEFF.dat POST_COEFF.dat || echo "file POST_COEFF.dat does not match"

numdiff \
    -a 0.000001 -r 0.05 \
    AUTOTEST/POST_OMEGA.dat POST_OMEGA.dat || echo "file POST_OMEGA.dat does not match"
