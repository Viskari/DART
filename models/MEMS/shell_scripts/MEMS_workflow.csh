#!/bin/tcsh

set model_dir = 'advance_model'

# Create a clean temporary directory and go there
rm -rf  $model_dir  || exit 1
mkdir -p $model_dir  || exit 1
cd       $model_dir  || exit 1

# Copy/link the associated programs

cp ../../ModelData/LUCAS.csv .
cp ../../ModelData/defaults_and_priors_MEMS.csv .

ln -s ../../ModelData/Initial_State.R .
ln -s ../../ModelData/MEMS_R_annual.R .
ln -s ../../ModelData/MEMS_steady.R .
ln -s ../../ModelData/NCinput.R .

Rscript Initial_State.R
ncgen filter_input.cdl -o filter_input.nc

cd ..
cp $model_dir/filter_input.nc .
cp ../Obs/obs_seq.out .
