# This is a wrapper file for model-vlbi.

#Importing libraries
import os, sys
import glob

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
# array = 'GMVA'
array = 'VLA'


vlbi_vp(scans, modes, station_abbrev_dic, dia_dic, blockage_diam='0.75m', max_rad='1.000deg')



# %%
perform_observation(scans, dia, station_names, pos_obs, modes, sources, integration_time)


# %%
mylengths, beam_max, lambd, pix_res = compute_beam_max_and_pix_res(x_adj_dic, y_adj_dic, z_adj_dic, modes)


# %%
source_fluxes_address = '../../source_fluxes.txt'
flux_file_existance_check(sources, start_date, source_fluxes_address)


# %%
flux_sources = flux_file_completion_check(sources, source_fluxes_address)


create_input_models(sources, pix_res, modes, flux_sources)

#Create noisy ms.
create_noisy_ms(scans, freq_setup, integration_time, npol, flux_sources)



#%%
corrupt_model_data(scans, pix_res)

#%%
combine_measurement_sets(scans)


    






          



