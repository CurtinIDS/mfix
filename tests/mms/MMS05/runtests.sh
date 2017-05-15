#!/bin/bash -elx

##### Change into the current working directory
####$ -cwd
#####
##### The name for the job. It will be displayed this way on qstat
####$ -N MMS05
#####
##### Number of cores to request
####$ -pe dev 64
#####
####$ -r n
#####
##### Queue Name
####$ -q dev

# set case directory
#setenv CASE_DIR `pwd`
export CASE_DIR=`pwd`

# copy common files
cp ../usr_common/usr_mod.f ./usr_mod.f
cp ../usr_common/usr3.f ./usr3.f

cp ../nonuniform_grids_3d/mesh_*.dat .

# compile MFIX
echo "******** Compiling MFIX..."
cd $CASE_DIR
../../../configure_mfix --enable-dmp FC=mpif90 FCFLAGS="-O0 -g -fcheck=all -DFLAG_MMS"
make -j

# remove these files if exist:
echo "******** Removing old files..."
# if [ -e [de_norms_collected.dat] ]; then rm de_norms_collected.dat; fi
rm -f de_norms_collected.dat

# create backup before adding user-defined grid spacing to input file
echo "******** Creating backup for mfix.dat..."
cp $CASE_DIR/mfix.dat $CASE_DIR/mfix_backup.dat

# Run mesh_8 (i.e., 8x8 for 2D, 8x8x8 for 3D)
echo "******** Running mesh_8..."
cat $CASE_DIR/mfix_backup.dat mesh_8.dat > mfix.dat
mpirun -np 1 ./mfix imax=8 jmax=8 kmax=8 DISCRETIZE=9*2 UR_FAC=9*0.05 > out.log
cat $CASE_DIR/de_norms.dat >> $CASE_DIR/de_norms_collected.dat
rm -f $CASE_DIR/{MMS05.*,de_norms.dat,out.log}
rm -f $CASE_DIR/solution_*.dat
#mkdir $CASE_DIR/mesh_8
#mv $CASE_DIR/solution_* $CASE_DIR/mesh_8/

# Run mesh_16 (i.e., 16x16 for 2D, 16x16x16 for 3D)
echo "******** Running mesh_16..."
cat $CASE_DIR/mfix_backup.dat mesh_16.dat > mfix.dat
mpirun -np 1 ./mfix imax=16 jmax=16 kmax=16 > out.log
cat $CASE_DIR/de_norms.dat >> $CASE_DIR/de_norms_collected.dat
rm -f $CASE_DIR/{MMS05.*,de_norms.dat,out.log}
rm -f $CASE_DIR/solution_*.dat
#mkdir $CASE_DIR/mesh_16
#mv $CASE_DIR/solution_* $CASE_DIR/mesh_16/

## Run mesh_32 (i.e., 32x32 for 2D, 32x32x32 for 3D)
#echo "******** Running mesh_32..."
#cat $CASE_DIR/mfix_backup.dat mesh_32.dat > mfix.dat
#mpirun -np 4 mfix imax=32 jmax=32 kmax=32 nodesi=2 nodesj=2 nodesk=1 > out.log
#cat $CASE_DIR/de_norms.dat >> $CASE_DIR/de_norms_collected.dat
#rm -f $CASE_DIR/{MMS05.*,de_norms.dat,out.log}
#rm -f $CASE_DIR/solution_*.dat
##mkdir $CASE_DIR/mesh_32
##mv $CASE_DIR/solution_* $CASE_DIR/mesh_32/

## Run mesh_64 (i.e., 64x64 for 2D, 64x64x64 for 3D)
#echo "******** Running mesh_64..."
#cat $CASE_DIR/mfix_backup.dat mesh_64.dat > mfix.dat
#mpirun -np 64 mfix imax=64 jmax=64 kmax=64 nodesi=4 nodesj=4 nodesk=4 > out.log
#cat $CASE_DIR/de_norms.dat >> $CASE_DIR/de_norms_collected.dat
#rm -f $CASE_DIR/{MMS05.*,de_norms.dat,out.log}
#rm -f $CASE_DIR/solution_*.dat
##mkdir $CASE_DIR/mesh_64
##mv $CASE_DIR/solution_* $CASE_DIR/mesh_64/
#
## Run mesh_128 (i.e., 128x128 for 2D, 128x128x128 for 3D)
#echo "******** Running mesh_128..."
#cat $CASE_DIR/mfix_backup.dat mesh_128.dat > mfix.dat
#mpirun -np 64 mfix imax=128 jmax=128 kmax=128 nodesi=4 nodesj=4 nodesk=4 > out.log
#cat $CASE_DIR/de_norms.dat >> $CASE_DIR/de_norms_collected.dat
#rm -f $CASE_DIR/{MMS05.*,de_norms.dat,out.log}
#rm -f $CASE_DIR/solution_*.dat
##mkdir $CASE_DIR/mesh_128
##mv $CASE_DIR/solution_* $CASE_DIR/mesh_128/

# Evaluate observed orders
cp ../usr_common/ooa_test.f95 $CASE_DIR
echo "******** Calculating observed orders..."
gfortran -o ooa_test ooa_test.f95
./ooa_test
rm -f $CASE_DIR/{ooa_test,ooa_test.f95}
rm -f $CASE_DIR/de_norms_collected.dat
#mv $CASE_DIR/de_norms_collected.dat $CASE_DIR/AUTOTEST/de_norms_collected.dat
#mv $CASE_DIR/de_l2.dat $CASE_DIR/AUTOTEST/de_l2.dat
#mv $CASE_DIR/de_linf.dat $CASE_DIR/AUTOTEST/de_linf.dat
#mv $CASE_DIR/ooa_l2.dat $CASE_DIR/AUTOTEST/ooa_l2.dat
#mv $CASE_DIR/ooa_linf.dat $CASE_DIR/AUTOTEST/ooa_linf.dat

rm -f $CASE_DIR/{usr_mod.f,usr3.f,mfix}
rm -f $CASE_DIR/mesh_*.dat

mv mfix_backup.dat mfix.dat

echo "******** Done."

ndselect \
  -b 1 -e 3 -s 1 AUTOTEST/de_l2.dat > tmp.dat
numdiff \
  -a 0.000001 -r 0.05 \
  de_l2.dat tmp.dat

ndselect \
  -b 1 -e 3 -s 1 AUTOTEST/de_linf.dat > tmp.dat
numdiff \
  -a 0.000001 -r 0.05 \
  de_linf.dat tmp.dat

ndselect \
  -b 1 -e 2 -s 1 AUTOTEST/ooa_l2.dat > tmp.dat
numdiff \
  -a 0.001 -r 0.05 \
  ooa_l2.dat tmp.dat

ndselect \
  -b 1 -e 2 -s 1 AUTOTEST/ooa_linf.dat > tmp.dat
numdiff \
  -a 0.001 -r 0.05 \
  ooa_linf.dat tmp.dat
