#Code to query ALMA Source Catalog to obtain the flux of sources

from datetime import datetime
import io
from astropy.io.votable import parse
from urllib.request import urlopen

#reading vex file address
from input_parameters import *

print("The address of vex file taken from 'input_parameters.py' is: " +
      vex_file_address)
print('\n')

sources = {}
start_date = 'undefined'
mode_name = 'undefined'
mode_freq = 'undefined'

print('The source names and coordinates read from the vex file are as follows:')

#reading vex file
with open(vex_file_address, 'r') as vex_file:
    # vex_file_lines = vex_file.readlines()
    for i, line in enumerate(vex_file.readlines()):

        #reading source names and coordinates
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
            # print('\n')
            
            sources[source_name] = [epoch, RA, Dec]


        #reading start date of observations
        if 'start=' in line and 'exper_nominal' not in line and start_date == 'undefined':

            line_parts = line.split('=')

            start_date = line_parts[1][:-7]
            start_date = datetime.strptime(start_date, '%Yy%jd%Hh%Mm%S')
            start_date = start_date.strftime('%d-%b-%Y')

            print('\n')
            print('The start date of observations is: ' + start_date + \
                  '. This will be used to obtain flux values.')
            print('\n')


#reading mode and frequency information
with open(vex_file_address, 'r') as vex_file:
    vex_file_lines = vex_file.readlines()
    for i, line in enumerate(vex_file_lines):

        if 'ref $PROCEDURES' in line and mode_name == 'undefined':

            if 'def' in vex_file_lines[i - 1]:
                mode_name = vex_file_lines[i - 1][4:-2]

            if 'def' in vex_file_lines[i - 2]:
                mode_name = vex_file_lines[i - 2][4:-2]

        if 'ref $FREQ' in line and mode_freq == 'undefined':

            line_parts_equal = line.split('=')
            mode_freq = line_parts_equal[1][1:12]
            if 'MHz' in mode_freq:
                mode_freq = mode_freq.split('.')[0]
                mode_freq = mode_freq + 'E+06'

            print('\nFor mode ' + mode_name +
                  ', the following frequency has been found: ' + mode_freq +
                  ' Hz')
            print('\n')


file_object = open("source_fluxes.txt", "w")

#specifying the header row in source_fluxes.txt
file_object.write('Source'.ljust(14) +  
                  'Date'.ljust(14) + 
                  'Freq(Hz)'.ljust(14) +
                  'ModelShape'.ljust(14) +
                  'Flux(Jy)'.ljust(14) +
                  'FluxError(Jy)'.ljust(14) +
                  'MajorAxis(arcsec)'.ljust(18) +
                  'MinorAxis(arcsec)'.ljust(18) +
                  'SpectralIndex'.ljust(18) +  
                  'SpectralIndexError'.ljust(22) + 
                  'Warning'.ljust(10) +  
                  'NearestMeasurementDate'.ljust(25) +  
                  'FluxEstimatorVersion'.ljust(10) +  
                  '\n')


# freq = '86E+09'


for source_name in sources.keys():

    #querying the ALMA source catalog
    query_url = ('https://almascience.eso.org/sc/flux?DATE={}&FREQUENCY={}&NAME={}'.
                 format(start_date, mode_freq, source_name))
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

    #writing source information if it is available
    if flux_density != '--':
        file_object.write(source_name.ljust(14) +  
                          start_date.ljust(14) + 
                          mode_freq.ljust(14) +
                          'point'.ljust(14) +
                          str(flux_density).ljust(14) +
                          str(flux_density_error).ljust(14) +
                          '0.001'.ljust(18) +
                          '0.001'.ljust(18) +
                          str(spectral_index).ljust(18) + 
                          str(spectral_index_error).ljust(22) +
                          str(data_conditions).ljust(10) + 
                          str(near_measure_date).ljust(25) + 
                          str(version).ljust(10) + 
                          '\n')

    #writing null source information if it is not available
    else:
        file_object.write(source_name.ljust(14) +  
                          start_date.ljust(14) + 
                          mode_freq.ljust(14) +
                          'point'.ljust(14) +
                          str(flux_density).ljust(14) +
                          str(flux_density_error).ljust(14) +
                          '0.001'.ljust(18) +
                          '0.001'.ljust(18) +
                          '\n')

file_object.close()


print("\nThe retrieved data as been saved in the file 'source_fluxes.txt'. " +
      "Please note that the symbol '--' means that the flux of " +
      "corresponding source has not been found at Calibrator Source Catalogue.")
