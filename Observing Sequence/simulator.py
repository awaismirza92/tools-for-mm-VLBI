# This is a wrapper file for model-vlbi.

#Importing libraries
import os, sys
import glob

#%%

#Specifying the working directory
current_dir = os.getcwd()
print(current_dir)

if 'small_ms' or 'per_scan_ms' in current_dir:
    file_address = '../../input_files/c171a'

else:
    file_address = 'input_files/c171a'


print(file_address)
#%%

#Executing the python file to extract observing sequence
execfile('../../observing_sequence_reader.py')

execfile('../../model_vlbi_functions.py')

# import model_vlbi_functions

# %%



   
# perform_observation(scans, dia, station_names, pos_obs, modes, integration_time)


# %%

# mylengths, beam_max, lambd, pix_res = compute_beam_max_and_pix_res(x_adj_dic, y_adj_dic, z_adj_dic, modes)



# %%
source_fluxes_address = '../../source_fluxes.txt'
flux_file_existance_check(sources, start_date, source_fluxes_address)


# %%


        
# flux_sources = flux_file_completion_check(sources, start_date, source_fluxes_address)
         
# create_input_models(sources, pix_res, modes, flux_sources)            
        
# #Create noisy ms.
create_noisy_ms(scans, freq_setup, integration_time, npol)
    

    
#%%

# corrupt_model_data(scans)     
    


#%%

# combine_measurment_sets(scans)


    






          



