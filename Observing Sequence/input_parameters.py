#file to input parameters for modelvlbi


#specify the address of vex file
vex_file_address = 'input_files/n21q1.vex'

#specify the address of sum file
sum_file_address = 'input_files/n21q1.sum'

#specify the address of directory containing stations antabfs files
antabfs_files_address = 'input_files/'

#specify the address of SCHED source catalog file
sources_gsfc2016_files_address = 'input_files/sources.gsfc2016.txt'

#specify the address of the file containing source fluxes
source_fluxes_address = 'source_fluxes.txt'

#specify the name of the array (e.g. 'GMVA', 'VLA', 'EVN')
array = 'EVN'

#specify the first and last scan to be simulated. To simulate all the scans, set
# start_scan = 1 and end_scan = 'until_last_scan'
start_scan = 2
end_scan = 3

