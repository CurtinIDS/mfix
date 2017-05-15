#!/bin/bash -exl

# Run case
echo "******** Running simulation..."
./mfix > out.log
rm -f $CASE_DIR/{TFM02.*,out.log}
#rm -f $CASE_DIR/de_norms.dat
rm -f $CASE_DIR/mfix

echo "******** Done."

numdiff \
    -a 0.000001 -r 0.05 \
    AUTOTEST/de_norms.dat de_norms.dat

# uncomment the following to generate plots:
#echo "******** Generating plots..."
#python plot_results.py &
