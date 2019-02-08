#!/bin/bash

run_num=$1
ln -sf DBfiles/db_sbs.gems.dat_oneMPD_noped db_sbs.gems.dat
analyzer -b -q 'replay.C('$run_num')' ;
./beam -t ped  -r $run_num;
ln -sf DBfiles/db_sbs.gems.dat_run_$run_num db_sbs.gems.dat
analyzer -b -q 'replay.C('$run_num')' ;
./beam -t rms -r $run_num;
cp DBfiles/analysis_rms.h_run_$run_num table_rms.h
./beam -r $run_num

