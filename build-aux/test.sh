#!/bin/bash -e

cases=$( find . -name mfix.dat | grep -v variable_ep_g_2d |grep -v "tests/dem"|grep -v "tests/mms"|grep -v "tests/tfm")

for case in ${cases}; do
        rundir=$(dirname $case)
        chmod -R u+w ${rundir}
        cp build-aux/Makefile.usr ${rundir}/Makefile
        make -C ${rundir} run_t0
done
