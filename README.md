
# Beam Test Analysis Software 
author : Tao Ye <tao.ye@stonybrook.edu> 

## Usage:
### beam [-t] [-r] [-f] [-c] [-o] [-h] 
### Options:
-	-h : Print help info 
-	-t <analysis_type>:  Optional, available options: ped, rms, ana, plot. 
	-	 -- *ana*: Default mode, reconstruction and  tracking analysis
	-	 -- ped: Generate a pedestal DB file for specific run
	-	 -- rms: Generate a pedestal RMS table file
	-	 -- plot: output plots only, useful fo debug 
-	-c <config_file>: optional, use beam.config by default 
-	-r <run_num>: mandatory, if input file is not specfied 
-	-f <input_filename>: mandatory, if <run_num> is not specfied 
-	-o <output_filename>: optional, specify output rootfile/plot name 
### Examples
```
	beam -r 411
	beam -r 411 -c user.conf	
	beam -t ped -r 411
```
