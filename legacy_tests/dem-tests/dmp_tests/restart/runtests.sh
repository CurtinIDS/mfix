#!/bin/bash -lex

rm -f DEM_RESTART* POST_posvel.dat
time -p mpirun -np 4 ./mfix TSTOP=1.0 NODESI=2 NODESJ=2 RUN_TYPE=\'NEW\'
diff -q POST_posvel.dat AUTOTEST/POST_posvel.dat

rm -f DEM_RESTART* POST_posvel.dat
time -p mpirun -np 4 ./mfix TSTOP=0.33549 NODESI=2 NODESJ=2 RUN_TYPE=\'NEW\'
time -p mpirun -np 4 ./mfix TSTOP=1.0 NODESI=2 NODESJ=2 RUN_TYPE=\'RESTART_1\'
diff -q POST_posvel.dat AUTOTEST/POST_posvel.dat
