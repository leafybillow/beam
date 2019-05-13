#!/bin/bash

run_num=$1
ln -sf db_sbs.gems.dat_oneMPD_noped db_sbs.gems.dat
# ln -sf db_sbs.gems.dat_twoMPD_noped db_sbs.gems.dat
analyzer -b -q 'replay.C('$run_num')' ;

../beam -t ped  -r $run_num -c oneMPD.conf
#../beam -t ped  -r $run_num -c twoMPD.conf

ln -sf $BEAM_DB_PATH\db_sbs.gems.dat_run_$run_num db_sbs.gems.dat
analyzer -b -q 'replay.C('$run_num')' ;


