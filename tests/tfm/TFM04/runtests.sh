#!/bin/bash -exl

# remove old result files
rm -f error_summary.dat
rm -f solution_tec_block.dat

# Run case
echo "******** Running simulation..."
#mpirun -np 1 ./mfix.exe nodesi=1 nodesj=1 nodesk=1
./mfix Discretize=9*0
rm -f TFM04.*

./mfix Discretize=9*3
rm -f TFM04.*

./mfix Discretize=9*2
rm -f TFM04.*

./mfix Discretize=9*5
rm -f TFM04.*

#./mfix Discretize=9*4
#rm -f TFM04.*

./mfix Discretize=9*7
rm -f TFM04.*

./mfix Discretize=9*6
rm -f TFM04.*

./mfix Discretize=9*8
rm -f TFM04.*

./mfix Discretize=9*9
rm -f TFM04.*

rm ./mfix

echo "******** Done."

numdiff \
    -a 0.000001 -r 0.05 \
    AUTOTEST/error_summary.dat error_summary.dat

# uncomment the following to generate plots:
#echo "******** Generating plots..."
#python plot_results.py &
