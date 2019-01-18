#!/bin/bash

run_num=$1
file=rootfiles/run_$run_num.root

ln -sf DBfiles/db_sbs.gems.dat_noped db_sbs.gems.dat ;
analyzer -b -q 'replay.C('$run_num')' ;
analyzer -b -q 'analysis_ped.C( "'$file'")';
analyzer -b -q 'replay.C('$run_num')' ;
analyzer -b -q 'analysis_rms.C("'$file'")' ;
cp DBfiles/analysis_rms.h_run_$run_num analysis_rms.h
analyzer -b -q 'analysis.C("'$file'")';