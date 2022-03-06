# This is a wrapper file for model-vlbi.

#Importing libraries
import os, sys
import glob
import matplotlib.pyplot as plt

#%%

#Specifying the working directory
current_dir = os.getcwd()
print('Current working directory: ' + current_dir + '\n')

if 'small_ms' or 'per_scan_ms' in current_dir:
    file_address = '../../input_files/c171a'

else:
    file_address = 'input_files/c171a'

print('Directory for GMVA session files: ' + file_address + '\n')

#%%
#Executing the python file to extract observing sequence
execfile('../../observing_sequence_reader.py')

execfile('../../model_vlbi_functions.py')


# Specifying the name of the array (e.g. 'GMVA')
array = 'GMVA'
# array = 'VLA'
image_dim = 1024

start_scan = 86
end_scan = 90

# modelvlbi_elevation(scans, x_adj_dic, y_adj_dic, z_adj_dic, sources,
#                         start_scan, end_scan)

# vlbi_plotelevation(scans, start_scan, end_scan, sources, usehourangle=True,
#                        elev_limit=0)




# vlbi_vp(scans, modes, station_abbrev_dic, dia_dic, start_scan, end_scan, blockage_diam='0.75m', max_rad='1.000deg')
# #
#
#
# # %%
# perform_observation(scans, station_names, modes, sources, integration_time,
#                     start_scan, end_scan)
# # # #
# # # #
# # # # # %%
mylengths, beam_max, lambd, pix_res = compute_beam_max_and_pix_res(x_adj_dic,
                                                                   y_adj_dic,
                                                                   z_adj_dic,
                                                                   modes)
# #
# #
# #
# # # # # %%
source_fluxes_address = '../../source_fluxes.txt'
flux_file_existance_check(sources, start_date, source_fluxes_address)


# # %%
flux_sources = flux_file_completion_check(sources, source_fluxes_address)
# # #
# # #

# fits_header_check(flux_sources)
user_model_images(flux_sources)

# create_input_models(sources, pix_res, modes, flux_sources)
# #
# # #Create noisy ms.
# # create_noisy_ms(scans, freq_setup, integration_time, npol, flux_sources,
# #                 start_scan, end_scan)
# #
# #
# #
# # #%%
# corrupt_model_data(scans, pix_res, start_scan, end_scan, modes)
# #
# # #%%
# combine_measurement_sets(scans, start_scan, end_scan)


    

print('\a')





          



