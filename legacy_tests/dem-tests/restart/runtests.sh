#!/bin/bash -lex

rm -f DEM_RESTART* POST_posvel.dat
time -p ./mfix TSTOP=1.0 RUN_TYPE=\'NEW\'
diff -q POST_posvel.dat AUTOTEST/POST_posvel.dat

rm -f DEM_RESTART* POST_posvel.dat
time -p ./mfix TSTOP=0.33549 RUN_TYPE=\'NEW\'
time -p ./mfix TSTOP=1.0 RUN_TYPE=\'RESTART_1\'
diff -q POST_posvel.dat AUTOTEST/POST_posvel.dat
