#!/bin/bash -lx

MFIX=${MFIX_HOME-"../../../"}

rm POST_*.dat &> /dev/null

RUN_NAME="DEM06"

rm -f ${RUN_NAME}* &> /dev/null
time -p ./mfix DES_ONEWAY_COUPLED=.T. \
    DES_INTERP_ON=.F. DES_INTERP_MEAN_FIELDS=.F.

rm -f ${RUN_NAME}* &> /dev/null
time -p ./mfix DES_ONEWAY_COUPLED=.T. \
    DES_INTERP_ON=.T. DES_INTERP_MEAN_FIELDS=.T. \
    DES_INTERP_SCHEME=\'GARG_2012\'

rm -f ${RUN_NAME}* &> /dev/null
time -p ./mfix DES_ONEWAY_COUPLED=.T. \
    DES_INTERP_ON=.T. DES_INTERP_MEAN_FIELDS=.T. \
    DES_INTERP_SCHEME=\'SQUARE_DPVM\' DES_INTERP_WIDTH=2.0d-3

rm -f ${RUN_NAME}* &> /dev/null
time -p ./mfix DES_ONEWAY_COUPLED=.F. \
    DES_INTERP_ON=.F. DES_INTERP_MEAN_FIELDS=.F.

rm -f ${RUN_NAME}* &> /dev/null
time -p ./mfix DES_ONEWAY_COUPLED=.F. \
    DES_INTERP_ON=.T. DES_INTERP_MEAN_FIELDS=.T. \
    DES_INTERP_SCHEME=\'GARG_2012\'

rm -f ${RUN_NAME}* &> /dev/null
time -p ./mfix DES_ONEWAY_COUPLED=.F. \
    DES_INTERP_ON=.T. DES_INTERP_MEAN_FIELDS=.T. \
    DES_INTERP_SCHEME=\'SQUARE_DPVM\' DES_INTERP_WIDTH=3.0d-3

rm -f ${RUN_NAME}* &> /dev/null
time -p ./mfix DES_ONEWAY_COUPLED=.F. \
    DES_INTERP_ON=.T. DES_INTERP_MEAN_FIELDS=.T. \
    DES_INTERP_SCHEME=\'SQUARE_DPVM\' DES_INTERP_WIDTH=4.0d-3

numdiff \
    -a 0.000001 -r 0.05 \
    --exclude=1:4 --exclude=2:4 \
    AUTOTEST/POST_VEL.dat POST_VEL.dat

numdiff \
    -a 0.000001 -r 0.05 \
    --exclude=1:4 --exclude=2:4 \
    AUTOTEST/POST_POS.dat POST_POS.dat
