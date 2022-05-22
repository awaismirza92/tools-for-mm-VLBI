# Code to obtain the flux of sources for low-frequency observations. Measurements
# of fluxes at 2.20 GHz & 8.40 GHz in SCHED source catalog are interpolated to get
# flux at desired frequency

from datetime import datetime
import numpy as np


# reading vex file address
from input_parameters import *
print("The address of vex file taken from 'input_parameters.py' is: " +
      vex_file_address + '\n')


# reading vex file
sources = {}
start_date = 'undefined'

print('The source names and coordinates read from the vex file are as follows:')

with open(vex_file_address, 'r') as vex_file:
    for i, line in enumerate(vex_file.readlines()):

        # reading source names and coordinates
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

            sources[source_name] = [epoch, RA, Dec]


# reading mode and frequency information
mode_name = 'undefined'
mode_freq = 'undefined'

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
                mode_freq = float(mode_freq)/1e3

            print("\nFor mode '" + mode_name +
                  "', the following frequency has been found: " + str(mode_freq) +
                  ' GHz')


# reading the date of generation of SCHED catalog
with open(sources_gsfc2016_files_address, 'r') as sched_catalog:
    sched_catalog_lines = sched_catalog.readlines()
    for i, line in enumerate(sched_catalog_lines):

        if 'Generated' in line:
            gen_time = line.split()[2]
            gen_time = datetime.strptime(gen_time, '%Y-%b-%d.').date()
            gen_time = gen_time.strftime("%d-%b-%Y")
            print(f"\nThe date of the generation of the SCHED catalog is: {gen_time}.")


file_object = open("source_fluxes.txt", "w")

# specifying the header row in source_fluxes.txt
file_object.write('Source'.ljust(14) +
                  'Date'.ljust(14) +
                  'Freq(Hz)'.ljust(14) +
                  'ModelShape'.ljust(14) +
                  'Flux(Jy)'.ljust(14) +
                  'FluxError(Jy)'.ljust(14) +
                  'MajorAxis(arcsec)'.ljust(18) +
                  'MinorAxis(arcsec)'.ljust(18) +
                  '\n')


# extracting fluxes from SCHEd catalog and using it for extrapolation
freqs = np.zeros(2)
fluxes = np.zeros(2)

with open(sources_gsfc2016_files_address, 'r') as sched_catalog:
    lines = sched_catalog.readlines()
    for source_name in sources.keys():
        found = False
        for i, line in enumerate(lines):
            if source_name in line:
                flux_line = lines[i+2]
                line_parts = flux_line.split(',')

                freqs[0] = line_parts[0].split('=')[1]
                fluxes[0] = line_parts[1]
                freqs[1] = line_parts[3]
                fluxes[1] = line_parts[4]

                fit_param = np.polyfit(freqs, fluxes, 1)
                fit_func = np.poly1d(fit_param)

                interpol_flux = fit_func(mode_freq)

                file_object.write(source_name.ljust(14) +
                                  gen_time.ljust(14) +
                                  str(mode_freq).ljust(14) +
                                  'point'.ljust(14) +
                                  str(round(interpol_flux, 3)).ljust(14) +
                                  'unknown'.ljust(14) +
                                  '0.001'.ljust(18) +
                                  '0.001'.ljust(18) +
                                  '\n')

                found = True

            if i == (len(lines) - 1) and found == False:
                file_object.write(source_name.ljust(14) +
                                  gen_time.ljust(14) +
                                  'point'.ljust(14) +
                                  '--'.ljust(14) +
                                  'unknown'.ljust(14) +
                                  '0.001'.ljust(18) +
                                  '0.001'.ljust(18) +
                                  '\n')

file_object.close()

print("\nThe retrieved data as been saved in the file 'source_fluxes.txt'. " +
      "Please note that the symbol '--' (if present) means that the flux of " +
      "corresponding source has not been found in SCHED Source Catalog.")
