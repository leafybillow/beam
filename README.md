
# GEM Tracking Analysis Software for Beam Test 
Author : Tao Ye
       <tao.ye@stonybrook.edu> 

## Usage:
 beam [OPTIONS] -r <run_number>
 
 beam [OPTIONS] -f <input_file>

### Options:

* -h : Print help info 
* -t <analysis_type>:  Optional, available options: ped, rms, ana. 
	* ana: Default mode, reconstruction and  tracking analysis
	* ped: Generate a pedestal DB file for specific run
* -c <config_file>: optional, use beam.config by default 
* -r <run_num>: mandatory, if input file is not specfied 
* -f <input_filename>: mandatory, if <run_num> is not specfied 
* -o <output_filename>: optional, specify output rootfile/plot name
* -e <number_of_events>: optional, analyze the first <number of events> only
* -S <number_of_event_shift>: optional, number event shift applied to QDC.
* -P:  output plots ONLY, no rootfile will be created. 

### Examples
```
	beam -r 411
	beam -r 411 -P
	beam -f rootfiles/run_411.root
	beam -r 411 -c user.conf	
	beam -t ped -r 411
	beam -r 411 -e 200
	beam -r 411 -S 32 
```
