#!/bin/bash -lex

if [ -f runtests.sh ]; then
    exec ./runtests.sh
fi

post_script=AUTOTEST/post.script.NEW

if [ -n "${MPIRANKS}" ]; then
    mpirun -np ${MPIRANKS} ./mfix${EXEEXT}
else
    ./mfix${EXEEXT}
fi
if [ -e ${post_script} ]; then
    OMP_NUM_THREADS=1 ./postmfix${EXEEXT} < ${post_script}
fi

post_dats=AUTOTEST/POST*.dat

for test_post_file in ${post_dats}; do
	  numdiff -a 0.000001 -r 0.05 ${test_post_file} $(basename ${test_post_file}) || echo "Post results differ"
done
