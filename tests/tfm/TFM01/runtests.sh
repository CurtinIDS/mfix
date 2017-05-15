#!/bin/bash -exl

# Run mesh_32 (i.e., 32x32 for 2D, 32x32x32 for 3D)
echo "******** Running mesh_32..."
./mfix imax=32 jmax=32 > out.log
rm -f TFM01.* out.log
#rm -f $CASE_DIR/de_norms.dat

echo "******** Done."

numdiff \
    -a 0.000001 -r 0.05 \
    AUTOTEST/de_norms.dat de_norms.dat

# uncomment the following to generate plots:
#echo "******** Generating plots..."
#python plot_results.py &
