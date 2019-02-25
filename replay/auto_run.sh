#!/bin/bash

run_num=$1
ln -sf db_sbs.gems.dat_oneMPD_noped db_sbs.gems.dat
analyzer -b -q 'replay.C('$run_num')' ;
../beam -t ped  -r $run_num -c ../beam.conf
ln -sf ../DBfiles/db_sbs.gems.dat_run_$run_num db_sbs.gems.dat
analyzer -b -q 'replay.C('$run_num')' ;
../beam -t rms -r $run_num -c ../beam.conf


