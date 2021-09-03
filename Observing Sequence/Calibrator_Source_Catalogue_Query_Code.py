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

file_object = open("calibrator_fluxes.txt", "w")

file_object.write('Source'.ljust(10) +  'Flux'.ljust(10) 
                  + 'Date'.ljust(14) + 'Freq'.ljust(10) + '\n')

freq = '86E+09'
for source_name in sources.keys():
    
    query_url = f'https://almascience.eso.org/sc/flux?DATE={start_date}&FREQUENCY={freq}&NAME={source_name}'
    url_response = urlopen(query_url)
    io_bytes_object = io.BytesIO(url_response.read())
    vo_table = parse(io_bytes_object)
    first_table = vo_table.get_first_table()
    flux_density = first_table.array['FluxDensity'][0]
    print(source_name, flux_density)
    
    if flux_density != '--':
        flux_density = round(flux_density, 2)
        file_object.write(source_name.ljust(10) +  str(flux_density).ljust(10) 
                          + start_date.ljust(14) + freq.ljust(10) + '\n')

    
file_object.close()