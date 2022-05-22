# This is the file which should be executed within CASA. It reads companion
# files included as part of model-vlbi.


# Importing libraries
import os
import sys
import glob
import pprint
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, timedelta
from simutil import simutil
u = simutil()



# Specifying the working directory
current_dir = os.getcwd()
print('Current working directory: ' + current_dir + '\n')

# Reading the address of input files
execfile('input_parameters.py')

# Reading station diameters and SEFDs
execfile('station_diameters_and_SEFDs.py')

# Reading observation sequence
execfile('observing_sequence_reader.py')

# Loading functions of model_vlbi
execfile('model_vlbi_functions.py')




image_dim = 1024

# Giving python indices for scans list
start_scan = start_scan - 1
# if type(end_scan) == int:
#     end_scan = end_scan - 1

# Setting index if user wants to simulate until last scan
if end_scan == 'until_last_scan':
    end_scan = len(scans) - 1



# vlbi_plotelevation(scans, start_scan, end_scan, sources, usehourangle=True,
#                    elev_limit=0)

# vlbi_vp(scans, modes, station_abbrev_dic, dia_dic, start_scan, end_scan, blockage_diam='0.75m', max_rad='1.000deg')
# #
#
#
# # %%
perform_observation(scans, station_names, modes, sources, integration_time,
                    start_scan, end_scan)
# # # # #
# # # # #
# # # # # # %%
mylengths, beam_max, lambd, pix_res = compute_beam_max_and_pix_res(x_adj_dic,
                                                                   y_adj_dic,
                                                                   z_adj_dic,
                                                                   modes)
# # #
# # # # # # %%
# # source_fluxes_address = '../../source_fluxes.txt'
# flux_file_existance_check(sources, start_date, source_fluxes_address)
# #
# #
# # # # %%
flux_sources = flux_file_completion_check(sources, source_fluxes_address)
# # # # #
# # # # #
# #
# # # fits_header_check(flux_sources)
# # user_model_images(flux_sources, pix_res)
#
# create_input_models(sources, pix_res, modes, flux_sources)
# #


# # #Create noisy ms.
# create_noisy_ms(scans, freq_setup, integration_time, npol, flux_sources,
#                 start_scan, end_scan)
# #
# #
# #
# # #%%
corrupt_model_data(scans, pix_res, start_scan, end_scan, modes)
# #
# # #%%
# combine_measurement_sets(scans, start_scan, end_scan)


print('\a')

