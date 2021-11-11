#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import numpy as np



def vlbi_vp(scans, modes, station_abbrev_dic, dia_dic,
            start_scan, end_scan,
               blockage_diam='0.75m',
               max_rad='1.000deg'):
    """ Function to compute and save the voltage pattern of a heterogeneous
    array of telescopes. Currently, only an airy disk voltage pattern can be
    computed with this function. simvlbi implements the table computed with this
    function to corrupt the visibilities. This function uses the
    vpmanager.setpbairy() tool in CASA.

    Parameters
    -----------
    configfile: str
            Configuration file containing coordinates (in Earth-centred
            left-handed coordinates) and diameter and SEFD of the telescopes in
            the array configuration. Eg: '<path-to-file>/<filename>.confg'.
    array: {'EHT', 'GMVA', other}
            Array of telescopes which is used for simulating the data. If using
            the expert mode with a different configuration of telescopes than
            'EHT' or 'GMVA', please specify your own array name and not 'EHT'
            or 'GMVA'.
    ref_freq: str
            Central frequency in GHz of spectral window of observation.
            Default: '230.609 GHz' for EHT and '86.0 GHz' for GMVA.
    projectname: str
            Name of the project with .ms file extension. The function computes
            voltage patterns for each telescope in the configuration file and
            saves them into a CASA table with the name '<projectname>_vp.tab'.
            Eg: 'simulation_<arrayname>.ms'
    blockage_diam: str, optional
            The effective diameter (in meters) of subreflector blockage that is
            used to compute the airy disk voltage pattern. '0.75 m' is the
            default value considered.
    max_rad: str, optional
            Maximum radial extent of the voltage pattern/primary beam (scales
            as 1/frequency) in degrees. This parameter is specified in the format
            '1.0 deg' for a 1 degree maximum radial extent. Default: '1.0 deg'.
    """

    # Read in the configuration file to obtain diameters of the antennae

    print('Creating voltage patterns for each scan configuration.')

    os.system('mkdir -p vp_tables')

    for scan in scans[start_scan:end_scan]:

        print(scan['scan_no'])

        os.system('rm -rf ' + 'vp_tables/scan_' + scan['scan_no'] + '_vp.tab')

        # Clear any voltage patterns previously determined in terminal.
        vp.reset()

        array = 'GMVA'

        # Create Voltage patterns for each telescope.
        # print("Computing voltage patterns..")
        for station_abbrev in scan['participating_stations']:
            station_name = station_abbrev_dic[station_abbrev]
            station_dia = dia_dic[station_name]
            vp.setpbairy(telescope=array, dishdiam=station_dia,
                         blockagediam=blockage_diam, maxrad=max_rad,
                         reffreq=modes[0]['freq'][0], dopb=True)

        # Save the estimated voltage patterns in a CASA table which will be accessed
        # to corrupt simulated data with primary beams
        # print("Voltage pattern is saved in a CASA table under the name "
        #       + projectname + "_vp.tab and is used to corrupt the visibilities.")
        vp.saveastable('vp_tables/scan_' + scan['scan_no'] + '_vp.tab')

    print('Done. \n')





#%%
def calculate_ref_position(x_part_stations, y_part_stations, z_part_stations):

    # print('Calculating reference position of the antenna configuration.')
    #Calculating reference position of the antenna configuration

    cofa_x = pl.average(x_part_stations)
    cofa_y = pl.average(y_part_stations)
    cofa_z = pl.average(z_part_stations)
    cofa_lat,cofa_lon,cofa_alt = u.xyz2long(cofa_x, cofa_y, cofa_z, 'WGS84')
    pos_obs = me.position("WGS84",qa.quantity(cofa_lon, "rad"), qa.quantity(cofa_lat, "rad"), qa.quantity(cofa_alt, "m"))

    # print('Done. \n')

    return pos_obs



def perform_observation(scans, dia, station_names, modes, sources,
                        integration_time, start_scan, end_scan):

    print('Perfoming observation for each scan with sm.observe. Scan number mentioned below.')

    for scan in scans[start_scan:end_scan]:
        # Creating new MS.
        ms_name = 'scan_' + scan['scan_no'] + '.ms'
        os.system('rm -rf ' + ms_name)
        sm.open(ms_name)

        vp.setpbantresptable(telescope=array, dopb=True,
                             antresppath='vp_tables/scan_' + scan['scan_no'] + '_vp.tab')


        x_part_stations = []
        y_part_stations = []
        z_part_stations = []
        dia_part_stations = []
        names_part_stations = []
        for station_abbrev in scan['participating_stations']:
            station_name = station_abbrev_dic[station_abbrev]
            x_part_stations.append(x_adj_dic[station_abbrev])
            y_part_stations.append(y_adj_dic[station_abbrev])
            z_part_stations.append(z_adj_dic[station_abbrev])
            dia_part_stations.append(dia_dic[station_name])
            names_part_stations.append(station_abbrev_dic[station_abbrev])

        scan_pos_obs = calculate_ref_position(x_part_stations, y_part_stations,
                                              z_part_stations)

        # Specifying the array configuration
        sm.setconfig(
            telescopename=array,
            x=x_part_stations,
            y=y_part_stations,
            z=z_part_stations,
            dishdiameter=dia_part_stations,
            mount='alt-az',
            antname=names_part_stations,
            coordsystem='global',
            referencelocation=scan_pos_obs)

        # Specifying spectral windows
        sm.setspwindow(
            spwname=scans[0]['mode'],
            freq=modes[0]['freq'][0],
            deltafreq=modes[0]['freq_resolution'][0],
            freqresolution=modes[0]['freq_resolution'][0],
            nchannels=modes[0]['n_channels'][0],
            stokes='RR LL')

        # Setting limits to flag data such as shadowing or when source is below
        # certain elevation
        elev_limit = 10.0
        sm.setlimits(shadowlimit=0.001, elevationlimit=str(elev_limit) + 'deg')

        # Specing feed and autocorrelation parameters
        sm.setfeed('perfect R L')
        sm.setauto(autocorrwt=0.0)

        print(scan['scan_no'])

        # Specifing the field to be observed
        sm.setfield(sourcename=scan['source_name'],
                    sourcedirection=sources[scan['source_name']])

        # Setting the integration time and reference time for observation
        sm.settimes(integrationtime=integration_time,
                    usehourangle=False,
                    referencetime=me.epoch('TAI', scan['start_time']))

        # Performing observations
        sm.observe(scan['source_name'], scans[0]['mode'],
                   starttime='0 s',
                   stoptime=str(scan['observation_time'][0]) + ' s')

    sm.close()
    print('Done.\n')


def compute_beam_max_and_pix_res(x_adj_dic, y_adj_dic, z_adj_dic, modes):
    # Computing the baseline lengths between the different antennae.
    # [lat,lon,el] (x, y, z) [m] are relative to the center
    # of the array, oriented with x and y tangent to the closest point at the
    # COFA (latitude, longitude) on the WGS84 reference ellipsoid, with z
    # normal to the ellipsoid and y pointing north.

    print('Compting maximum baseline and pixel resolution.')

    cx = np.mean(x_adj_dic.values())
    cy = np.mean(y_adj_dic.values())
    cz = np.mean(z_adj_dic.values())
    n_antennas = len(x_adj_dic)
    lat, lon, el = u.itrf2loc(x_adj_dic.values(), y_adj_dic.values(),
                              z_adj_dic.values(), cx, cy, cz)
    mylengths = np.zeros(n_antennas * (n_antennas - 1) / 2)
    k = 0
    for i in range(n_antennas):
        for j in range(i + 1, n_antennas):
            x = lat[i] - lat[j]
            y = lon[i] - lon[j]
            z = el[i] - el[j]
            mylengths[k] = np.sqrt((x ** 2 + y ** 2 + z ** 2))
            k = k + 1

    # Calculating the largest baseline length in m
    beam_max = np.round(np.amax(mylengths), 4)

    # Computing the wavelenght of observations
    lambd = np.round(300 / qa.quantity(modes[0]['freq'][0], 'MHz')['value'], 4)  # in m

    # Calculating pixel resolution
    pix_res = lambd / beam_max * 180 / np.pi * 3600
    pix_res = '{}'.format(pix_res) + 'arcsec'

    print('Done. \n')

    return mylengths, beam_max, lambd, pix_res


def flux_file_existance_check(sources, start_date, source_fluxes_address,
                              freq='86E+09'):

    print("Finding the file 'source_fluxes.txt'.")

    try:
        with open(source_fluxes_address, 'r') as flux_file:
            print("Found. \n")

    except:
        print("Error: The file 'source_fluxes.txt' was not found in the current " +
              "working directory: %s. \n" % os.getcwd())

        template_file = open("source_fluxes_template.txt", "w")
        template_file.write('Source'.ljust(14) +
                            'Date'.ljust(14) +
                            'Freq(Hz)'.ljust(14) +
                            'ModelShape'.ljust(14) +
                            'Flux(Jy)'.ljust(14) +
                            '\n')
        for source_name in sources.keys():
            template_file.write(source_name.ljust(14) +
                                start_date.ljust(14) +
                                freq.ljust(14) +
                                'point'.ljust(14) +
                                '99'.ljust(14) +
                                '\n')
        template_file.close()

        print("Please use the code 'Calibrator_Source_Catalogue_Query_Code.py' to" +
              " obtain 'source_fluxes.txt' file for your sources. Otherwise, please " +
              "make the file manually yourself. A template file " +
              "'source_fluxes_template.txt' has been generated for your reference. " +
              "Please note that the template has incorrect flux values (99 Jy) in it.")
        print('\nThe code is terminating now. \n')
        sys.exit()


def flux_file_completion_check(sources, source_fluxes_address):

    print("Checking the completion of the file 'source_fluxes.txt'.")

    initial_print = True

    flux_sources = []
    flux_sources_names = []
    with open(source_fluxes_address, 'r') as flux_file:
        lines_list = flux_file.readlines()

        for line in lines_list[1:]:

            line_parts = line.split()

            flux_sources_names.append(line_parts[0])

            flux_source = {'source_name': line_parts[0]}
            flux_source['date'] = line_parts[1]
            flux_source['freq'] = line_parts[2]
            flux_source['ModelShape'] = line_parts[3]
            flux_source['flux'] = float(line_parts[4])

            flux_sources.append(flux_source)

            if initial_print:
                print("The following sources have been found in " +
                      "'source_fluxes.txt': ")
                initial_print = False

            print(flux_source)

            if flux_source['flux'] == '--':
                print("The flux for % is -- which" % flux_source['source_name'] +
                      " cannot be correct. Please correct it." +
                      "\n\nThe code is terminating now\n")
                sys.exit()

    count = 0
    for observe_source in sources.keys():
        if observe_source not in flux_sources_names:
            if count == 0:
                print("\nThe following sources are missing in 'source_fluxes.txt':")
                count += 1
            print(observe_source)

    if count == 1:
        print("Please mention the flux for these sources in " +
              "'source_fluxes.txt'.")
        print('\nThe code is terminating now.\n')
        sys.exit()

    else:
        print("All the observation sources are present in 'source_fluxes.txt'")
        print('Done. \n')
        return flux_sources




def create_input_models(sources, pix_res, modes, flux_sources):

    print('Creating models of sources.')

    os.system('mkdir -p model_images')

    # Create model image for the sources
    for source_name, source_direction in sources.items():

        print("Creating the model image with the filename: {} ."
              .format(source_name + '.modelimage.im'))

        cl.done()

        vp.setpbantresptable(telescope=array, dopb=True,
                         antresppath='vp_tables/scan_' + scan['scan_no'] + '_vp.tab')

        input_model = 'point'  # A string that could be either point,
        # Gaussian, disk, or limbdarkeneddisk

        for flux_source in flux_sources:
            if flux_source['source_name'] == source_name:
                flux = flux_source['flux']
                input_model = flux_source['ModelShape']

        cl.addcomponent(dir=source_direction, flux=flux, fluxunit='Jy',
                        freq=modes[0]['freq'][0], shape=input_model)

        ia.fromshape(outfile='model_images/' + source_name + '.modelimage.im',
                     shape=[512, 512, 1, 1],
                     overwrite=True)
        cs = ia.coordsys()
        cs.setunits(['rad', 'rad', '', 'Hz'])

        cell_rad = qa.convert(pix_res, "rad")['value']
        cs.setincrement([-cell_rad, cell_rad], 'direction')
        cs.setreferencevalue([qa.convert(source_direction[1], 'rad')['value'],
                              qa.convert(source_direction[2], 'rad')['value']],
                             type="direction")

        cs.setreferencevalue(modes[0]['freq'][0], 'spectral')
        cs.setincrement(qa.convert(modes[0]['freq_resolution'][0], 'GHz'),
                        'spectral')
        ia.setcoordsys(cs.torecord())
        ia.setbrightnessunit("Jy/pixel")
        ia.modify(cl.torecord(), subtract=False)
        exportfits(imagename='model_images/' + source_name + '.modelimage.im',
                   fitsimage='model_images/' + source_name + '.modelimage.fits',
                   overwrite=True)

        ia.close()
        cl.done()

    print('Done. \n')




def create_noisy_ms(scans, freq_setup, integration_time, npol, flux_sources,
                    start_scan, end_scan):
    # print("Visibilities added to Data column of MS. Now corrupting visibilities "
    #       + "with previously computed rms noise. New noisy MS will be created.")

    print('Creating noisy measurement set for each scan.')

    nbits = 2.0  # It's standard to have two bits
    eta_c = (0.88 if nbits == 2 else 0.64)  # Correlation efficient remains 0.88 if nbits = 2
    # nbase = len(mylengths)                  # Number of baselines

    sigma_array = []
    for scan in scans[start_scan:end_scan]:

        print('\n')
        print(scan['scan_no'])
        print(scan['participating_stations'])

        vp.setpbantresptable(telescope=array, dopb=True,
                             antresppath='vp_tables/scan_' + scan['scan_no'] + '_vp.tab')

        for flux_source in flux_sources:
            if scan['source_name'] in flux_source.values():
                print('Source flux: ' + str(flux_source['flux']))


        participating_stations_sefds = []
        for station in scan['participating_stations']:
            participating_stations_sefds.append(SEFD_dic[station])

        # Calculating overall SEFD
        participating_stations_sefds = np.array(participating_stations_sefds)
        sefd_reciprocal = 1. / np.outer(participating_stations_sefds,
                                        participating_stations_sefds)
        sefd_reciprocal = np.tril(sefd_reciprocal)
        np.fill_diagonal(sefd_reciprocal, 0)
        sefd_star = 1 / np.sqrt(np.sum(sefd_reciprocal))

        # Calculating baselines in the scan
        number_participating_stations = len(scan['participating_stations'])
        nbase = (number_participating_stations *
                 (number_participating_stations - 1) / 2)

        # Making Casa quantities of bandwidth per channel and integration time
        bandwidth_per_channel = qa.quantity(freq_setup[-2], 'MHz')
        integration_time_qa = qa.quantity(integration_time, 's')

        bandwidth_per_channel = qa.convert(bandwidth_per_channel, 'Hz')

        # Calculating simple noise for the scan
        sigma_simple = ((1 / eta_c) * sefd_star * np.sqrt(npol * nbase) /
                        np.sqrt(integration_time_qa['value'] *
                                bandwidth_per_channel['value'] * nbits))
        print('Sigma per scan: ', sigma_simple)
        sigma_array.append(sigma_simple)
        sigma_simple = qa.quantity(sigma_simple, 'Jy')

        print(scan['source_name'], scan['start_time'], scan['observation_time'])

        # Specifying name of noisy measurement set
        ms_name = 'scan_' + scan['scan_no'] + '.ms'
        noisy_ms_name = ms_name[:-3] + '.noisy.ms'

        # Reading files present in the current working directory
        current_files = glob.glob('*')

        # Deleting existing (if any) noisy measurement set
        if noisy_ms_name in current_files:
            os.system('rm -fr ' + noisy_ms_name)
            print('Removed previously present: ' + noisy_ms_name)

        # Creating a noisy measurement set
        os.system('cp -r ' + ms_name + ' ' + noisy_ms_name)
        sm.openfromms(noisy_ms_name)
        sm.setnoise(mode='simplenoise', simplenoise=sigma_simple)
        scan_vptable = 'vp_tables/scan_' + scan['scan_no'] + '_vp.tab'
        sm.setvp(dovp=True, usedefaultvp=False, vptable=scan_vptable)
        # To invoke uv-domain primary beam convolution for heterogeneous arrays we
        # set the fourier transform machine to mosaic.

        sm.setoptions(ftmachine="mosaic")
        sm.corrupt()

        sm.done()
        sm.close()

    print('Done. \n')

    # print(ms_name
    #       + ".noisy.ms has been created which contains the visibilities of "
    #       + ms_name + ".ms corrupted with thermal noise.")





def corrupt_model_data(scans, pix_res, start_scan, end_scan):
    # Corrupting model image with voltage pattern, copy it to data column of
    # Measurement Set.

    print('Adding model data into noisy measurment sets.')

    for scan in scans[start_scan:end_scan]:
        # print("Corrupting model data with voltage pattern which is estimated "
        #       + "separately for each antenna in the array.")

        print(scan['scan_no'])

        vp.setpbantresptable(telescope=array, dopb=True,
                             antresppath='vp_tables/scan_' + scan['scan_no'] + '_vp.tab')

        ms_name = 'scan_' + scan['scan_no'] + '.ms'
        noisy_ms_name = ms_name[:-3] + '.noisy.ms'

        im.open(noisy_ms_name)
        im.selectvis()
        im.defineimage(nx=512, ny=512, cellx=pix_res, celly=pix_res,
                       facets=1)
        im.setoptions(ftmachine='mosaic')
        im.ft(model='model_images/' + scan['source_name'] + '.modelimage.im', incremental=False)
        im.close()
        im.open(thems=noisy_ms_name)
        uvsub(vis=noisy_ms_name, reverse=True)
        # uvsub(vis=noisy_ms_name)
        im.close()

    print('Done. \n')

def combine_measurement_sets(scans, start_scan, end_scan):

    print("Combining individual scan noisy measurment sets into a single " +
          " measurment set named 'all_scans_combined.ms'.")

    noisy_ms_name_list = []
    for scan in scans[start_scan:end_scan]:
        ms_name = 'scan_' + scan['scan_no'] + '.ms'
        noisy_ms_name = ms_name[:-3] + '.noisy.ms'

        noisy_ms_name_list.append(noisy_ms_name)

    combined_ms_name = 'all_scans_combined.ms'

    # Reading files present in the current working directory
    current_files = glob.glob('*')

    # Deleting existing (if any) noisy measurement set
    if combined_ms_name in current_files:
        os.system('rm -r ' + combined_ms_name)
        print('Removed previously present: ' + combined_ms_name)

    concat(vis=noisy_ms_name_list, concatvis=combined_ms_name)

    print('Done. \n')
