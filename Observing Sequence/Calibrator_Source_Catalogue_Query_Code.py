from datetime import datetime
import io


from astropy.io.votable import parse
from urllib.request import urlopen


file_address = 'input_files/c171a'

sources = {}
start_date = 'undefined'

with open(file_address + '.vex.difx', 'r') as vex_file:

    for line in vex_file.readlines():
 
        if 'source_name' in line:
            source_name = line.split('= ')[1][:-2]
            direction = [source_name]

        if '*' not in line and 'ra =' in line:
            coord_parts = line.split('= ')
            
            epoch = coord_parts[-1][:-2]
            direction.append(epoch)
            
            RA = coord_parts[1][:-6]
            direction.append(RA)
            
            Dec = coord_parts[2][:-18].replace("'", 'm').replace('"', 's')
            direction.append(Dec)
            
            print(direction)
            print('\n')
            
            sources[source_name] = [epoch, RA, Dec]
            
        if 'start=' in line and 'exper_nominal' not in line and start_date == 'undefined':

            line_parts = line.split('=')

            start_date = line_parts[1][:-7]
            start_date = datetime.strptime(start_date, '%Yy%jd%Hh%Mm%S')
            start_date = start_date.strftime('%d-%b-%Y')
            
            print(start_date)
            print('\n')

#%%
file_object = open("source_fluxes.txt", "w")

file_object.write('Source'.ljust(14) +  
                  'Date'.ljust(14) + 
                  'Freq(Hz)'.ljust(14) +
                  'ModelShape'.ljust(14) +
                  'Flux(Jy)'.ljust(14) +
                  'FluxError(Jy)'.ljust(14) + 
                  'SpectralIndex'.ljust(18) +  
                  'SpectralIndexError'.ljust(22) + 
                  'Warning'.ljust(10) +  
                  'NearestMeasurementDate'.ljust(25) +  
                  'FluxEstimatorVersion'.ljust(10) +  
                  '\n')

freq = '86E+09'

for source_name in sources.keys():
    
    query_url = f'https://almascience.eso.org/sc/flux?DATE={start_date}&FREQUENCY={freq}&NAME={source_name}'
    url_response = urlopen(query_url)
    io_bytes_object = io.BytesIO(url_response.read())
    vo_table = parse(io_bytes_object)
    first_table = vo_table.get_first_table()
    
    flux_density = first_table.array['FluxDensity'][0].round(5)
    flux_density_error = first_table.array['FluxDensityError'][0].round(5)
    
    spectral_index = first_table.array['SpectralIndex'][0].round(5)
    spectral_index_error = first_table.array['SpectralIndexError'][0].round(5)
    
    data_conditions = first_table.array['DataConditions'][0]
    near_measure_date = first_table.array['Nearest Measurement Date'][0]
    
    version = first_table.array['Version'][0]
    
    print(source_name.ljust(14), str(flux_density).ljust(14))
    
    if flux_density != '--':
        file_object.write(source_name.ljust(14) +  
                          start_date.ljust(14) + 
                          freq.ljust(14) +
                          'point'.ljust(14) +
                          str(flux_density).ljust(14) +
                          str(flux_density_error).ljust(14) +
                          str(spectral_index).ljust(18) + 
                          str(spectral_index_error).ljust(22) +
                          str(data_conditions).ljust(10) + 
                          str(near_measure_date).ljust(25) + 
                          str(version).ljust(10) + 
                          '\n')
    else:
        file_object.write(source_name.ljust(14) +  
                          start_date.ljust(14) + 
                          freq.ljust(14) +
                          'point'.ljust(14) +
                          str(flux_density).ljust(14) +
                          '\n')

    
file_object.close()


print("\nThe retrieved data as been saved in the file 'source_fluxes.txt'. " +
      "Please note that the symbol '--' means that the flux of " +
      "corresponding source has not been found at Calibrator Source Catalogue.")
