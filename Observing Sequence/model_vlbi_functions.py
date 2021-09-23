#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys


def perform_observation(scans, dia, station_names, pos_obs, modes, integration_time):
    
    # Specifying the name of the array (e.g. 'GMVA') for sm.setconfig later on
    # array = 'GMVA'
    array = 'VLA'

    for scan in scans[80:100]:

        # Creating new MS.        
        ms_name = 'scan_' + scan['scan_no'] + '.ms'        
        os.system('rm -rf ' + ms_name)  
        sm.open(ms_name)
        
        #Specifying the array configuration
        sm.setconfig(
        telescopename = array,
        x = x_adj_dic.values(),
        y = y_adj_dic.values(),
        z = z_adj_dic.values(),
        dishdiameter = dia,
        mount = 'alt-az',
        antname = station_names,
        coordsystem = 'global',
        referencelocation = pos_obs)
    
        #Specifying spectral windows
        sm.setspwindow(
        spwname = scans[0]['mode'], 
        freq = modes[0]['freq'][0],                      
        deltafreq = modes[0]['freq_resolution'][0],      
        freqresolution = modes[0]['freq_resolution'][0], 
        nchannels = modes[0]['n_channels'][0],           
        stokes = 'RR LL')             
    
 
        #Setting limits to flag data such as shadowing or when source is below 
        #certain elevation
        elev_limit = 10.0      
        sm.setlimits(shadowlimit = 0.001, elevationlimit = str(elev_limit)+'deg')   
        
        
        #Specing feed and autocorrelation parameters
        sm.setfeed('perfect R L')         
        sm.setauto(autocorrwt = 0.0)  
        

        
        print('\n')
        print(scan['scan_no'])
        
        #Specifing the field to be observed
        sm.setfield(sourcename = scan['source_name'], 
                    sourcedirection = sources[scan['source_name']])
          
    
        #Setting the integration time and reference time for observation
        sm.settimes(integrationtime = integration_time, 
                    usehourangle = False, 
                    referencetime = me.epoch('TAI', scan['start_time']))   
    
    	#Performing observations
        sm.observe(scan['source_name'], scans[0]['mode'], 
                    starttime = '0 s', 
                    stoptime = str(scan['observation_time'][0]) + ' s')
    
        
    sm.close()
    
    
    
    
def compute_beam_max_and_pix_res(x_adj_dic, y_adj_dic, z_adj_dic, modes):

    # Computing the baseline lengths between the different antennae.  
    # [lat,lon,el] (x, y, z) [m] are relative to the center
    # of the array, oriented with x and y tangent to the closest point at the
    # COFA (latitude, longitude) on the WGS84 reference ellipsoid, with z
    # normal to the ellipsoid and y pointing north.
    
    cx = np.mean(x_adj_dic.values())
    cy = np.mean(y_adj_dic.values())
    cz = np.mean(z_adj_dic.values())
    nAntennas = len(x_adj_dic)
    lat, lon, el = u.itrf2loc(x_adj_dic.values(), y_adj_dic.values(), 
                              z_adj_dic.values(), cx, cy, cz)
    mylengths = np.zeros(nAntennas * (nAntennas - 1) / 2)
    k = 0
    for i in range(nAntennas):
        for j in range(i + 1, nAntennas):
            x = lat[i] - lat[j]
            y = lon[i] - lon[j]
            z = el[i] - el[j]
            mylengths[k] = np.sqrt((x**2 + y**2 + z**2))
            k = k + 1
            
            
    # Calculating the largest baseline length in m
    beam_max = np.round(np.amax(mylengths),  4)  
    
    # Computing the wavelenght of observations
    lambd = np.round(300 / qa.quantity(modes[0]['freq'][0], 'MHz')['value'], 4)  # in m
    
    #Calculating pixel resolution
    pix_res = lambd / beam_max * 180 / np.pi * 3600 
    pix_res = '{}'.format(pix_res) + 'arcsec'
    
    return mylengths, beam_max, lambd, pix_res






def flux_file_existance_check(sources, start_date, source_fluxes_address, 
                              freq = '86E+09'):
    try:
        with open(source_fluxes_address, 'r') as flux_file:
            print("\nThe file 'source_fluxes.txt' has been found")
            
    except:
        print("Error: The file 'source_fluxes.txt' was not found in the current " +
             "working directory: %s. \n" % os.getcwd())
        
        template_file = open("source_fluxes_template.txt", "w")
        template_file.write('Source'.ljust(14) +  
                      'Date'.ljust(14) + 
                      'Freq(Hz)'.ljust(14) + 
                      'Flux(Jy)'.ljust(14) + 
                      'ModelShape'.ljust(14) + 
                      '\n')
        for source_name in sources.keys():
            template_file.write(source_name.ljust(14) +  
                              start_date.ljust(14) + 
                              freq.ljust(14) + 
                              '99'.ljust(14) + 
                              'point'.ljust(14) + 
                              '\n')
        template_file.close()
            
        print("Please use the code 'Calibrator_Source_Catalogue_Query_Code.py' to" +
              " obtain 'source_fluxes.txt' file for your sources. Otherwise, please " +
              "make the file manually yourself. A template file " + 
              "'source_fluxes_template.txt' has been generated for your reference. " +
              "Please note that the template has incorrect flux values (99 Jy) in it.")
        print('\nThe code is terminating now.')
        sys.exit()
        
   
        
   
    
def flux_file_completion_check(sources, start_date, source_fluxes_address):
        
    initial_print = True
    
    flux_sources = []
    flux_sources_names = []
    with open(source_fluxes_address, 'r') as flux_file:
        lines_list = flux_file.readlines()
        
        for line in lines_list[1:]:
            
            line_parts = line.split()
            
            flux_sources_names.append(line_parts[0])
            
            flux_source = {'source_name' : line_parts[0]}
            flux_source['date'] = line_parts[1]
            flux_source['freq'] = line_parts[2]
            flux_source['flux'] = float(line_parts[3])
            flux_source['ModelShape'] = line_parts[4]
            
            flux_sources.append(flux_source)
            
            if initial_print == True:
                print("The following sources have been found in " + 
                      "'source_fluxes.txt': ")
                initial_print = False
                
            print(flux_source)
            
            if flux_source['flux'] == '--':
                print("The flux for % is -- which" % flux_source['source_name'] +
                      " cannot be correct. Please correct it."+
                      "\n\nThe code is terminating now")
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
        print('\nThe code is terminating now.')
        sys.exit()
        
    else:
         print("All the observation sources are present in 'source_fluxes.txt'\n")
         return flux_sources
    
         
         
         
         
def create_input_models(sources, pix_res, modes, flux_sources):  



    #Create model image for the sources
    for source_name, source_direction in sources.items():
        
        print("Creating the model image with the filename: {} ."
              .format(source_name + '.modelimage.im'))
        cl.done()
        
        input_model = 'point' #A string that could be either point, 
                              #Gaussian, disk, or limbdarkeneddisk
                              
        for flux_source in flux_sources:
            if flux_source['source_name'] == source_name:
                flux = flux_source['flux']
                input_model = flux_source['ModelShape']
                              
        cl.addcomponent(dir = source_direction, flux = flux, fluxunit = 'Jy', 
                        freq = modes[0]['freq'][0], shape = input_model)
        
        ia.fromshape(outfile = 'model_images/' + source_name+'.modelimage.im', 
                     shape = [512,512,1,1], 
                     overwrite = True)
        cs = ia.coordsys()
        cs.setunits(['rad','rad','','Hz'])
        
        
        cell_rad = qa.convert(pix_res,"rad")['value']
        cs.setincrement([-cell_rad, cell_rad], 'direction')
        cs.setreferencevalue([qa.convert(source_direction[1],'rad')['value'],
                              qa.convert(source_direction[2],'rad')['value']], 
                             type = "direction")
        
        cs.setreferencevalue(modes[0]['freq'][0], 'spectral')
        cs.setincrement(qa.convert(modes[0]['freq_resolution'][0], 'GHz'), 
                        'spectral')
        ia.setcoordsys(cs.torecord())
        ia.setbrightnessunit("Jy/pixel")
        ia.modify(cl.torecord(), subtract=False)
        exportfits(imagename = 'model_images/' + source_name+'.modelimage.im', 
                   fitsimage = 'model_images/' + source_name+'.modelimage.fits',
                   overwrite=True)
        
        ia.close()
        cl.done()         
        
        
        
        
def create_noisy_ms(scans, freq_setup, integration_time, npol):
    
    print("Visibilities added to Data column of MS. Now corrupting visibilities "
          + "with previously computed rms noise. New noisy MS will be created.")
    
    
    nbits = 2.0                             #It's standard to have two bits
    eta_c = (0.88 if nbits == 2 else 0.64)  #Correlation efficient remains 0.88 if nbits = 2
    # nbase = len(mylengths)                  # Number of baselines
       
    
    sigma_array = []
    for scan in scans[80:100]:
        
        print('\n')
        print(scan['scan_no'])        
        print(scan['participating_stations'])
        
        participating_stations_SEFDs = []
        for station in scan['participating_stations']:
            participating_stations_SEFDs.append(SEFD_dic[station])            
        
        #Calculating overall SEFD
        participating_stations_SEFDs = np.array(participating_stations_SEFDs)       
        SEFD_reciprocal = 1. / np.outer(participating_stations_SEFDs, 
                                        participating_stations_SEFDs)
        SEFD_reciprocal = np.tril(SEFD_reciprocal)
        np.fill_diagonal(SEFD_reciprocal, 0)
        SEFD_star = 1 / np.sqrt(np.sum(SEFD_reciprocal))
        print('SEFD_star:', SEFD_star)
            
        #Calculating baselines in the scan
        number_participating_stations = len(scan['participating_stations'])
        nbase = (number_participating_stations * 
                 (number_participating_stations - 1) / 2)
        print('Baselines: ', nbase)
        
        #Making Casa quantities of bandwidth per channel and integration time
        bandwidth_per_channel = qa.quantity(freq_setup[-2], 'MHz')
        integration_time_qa = qa.quantity(integration_time, 's')
        
        bandwidth_per_channel = qa.convert(bandwidth_per_channel, 'Hz')
        print(bandwidth_per_channel)
        
        #Calculating simple noise for the scan
        sigma_simple = ((1/eta_c) * SEFD_star * np.sqrt(npol * nbase) / 
                        np.sqrt(integration_time_qa['value'] * 
                                bandwidth_per_channel['value'] * nbits))
        print('Sigma per scan: ', sigma_simple)        
        sigma_array.append(sigma_simple)
        sigma_simple = qa.quantity(sigma_simple, 'Jy')
        
        print(scan['source_name'], scan['start_time'], scan['observation_time'])
        
        #Specifying name of noisy measurement set
        ms_name = 'scan_' + scan['scan_no'] + '.ms'
        noisy_ms_name = ms_name[:-3] + '.noisy.ms'   
        
        #Reading files present in the current working directory
        current_files = glob.glob('*')  
        
        #Deleting existing (if any) noisy measurement set
        if noisy_ms_name in current_files:
            os.system('rm -r ' + noisy_ms_name)
            print('Removed previously present: ' + noisy_ms_name)
        
        
        #Creating a noisy measurement set
        os.system('cp -r ' + ms_name + ' ' + noisy_ms_name)
        sm.openfromms(noisy_ms_name)
        sm.setnoise(mode = 'simplenoise', simplenoise = sigma_simple)
        sm.corrupt()
        
        
        # To invoke uv-domain primary beam convolution for heterogeneous arrays we
        # set the fourier transform machine to mosaic.
        
        sm.setoptions(ftmachine="mosaic")
        sm.done()
        sm.close()
        

    print(sigma_array)
        
    
    print(ms_name
    + ".noisy.ms has been created which contains the visibilities of "
    + ms_name + ".ms corrupted with thermal noise.")
    
    
    
    
    
    
def corrupt_model_data(scans):

    # Corrupting model image with voltage pattern, copy it to data column of
    # Measurement Set.
        
    for scan in scans[80:100]:
    
        print("Corrupting model data with voltage pattern which is estimated "
              + "separately for each antenna in the array.")
        
        ms_name = 'scan_' + scan['scan_no'] + '.ms'
        noisy_ms_name = ms_name[:-3] + '.noisy.ms'
        
        im.open(noisy_ms_name)
        im.selectvis()
        im.defineimage(nx = 512, ny = 512, cellx = pix_res, celly = pix_res, 
                       facets = 2)
        im.setvp(dovp = True, usedefaultvp = True)
        im.setoptions(ftmachine = 'mosaic')
        # im.ft(model = 'model_images/' + 'Einstein.im')
        im.ft(model = 'model_images/' + scan['source_name'] + '.modelimage.im')
        im.close()
        im.open(thems = noisy_ms_name)
        uvsub(vis = noisy_ms_name, reverse = True)
        im.close()
        
        
def combine_measurment_sets(scans):
    
    noisy_ms_name_list = []
    for scan in scans[80:100]:
        ms_name = 'scan_' + scan['scan_no'] + '.ms'
        noisy_ms_name = ms_name[:-3] + '.noisy.ms'
        
        noisy_ms_name_list.append(noisy_ms_name)
        
    print(noisy_ms_name_list)
    
    combined_ms_name = 'all_scans_combined.ms'
    
    #Reading files present in the current working directory
    current_files = glob.glob('*')  
    
    #Deleting existing (if any) noisy measurement set
    if noisy_ms_name in current_files:
        os.system('rm -r ' + combined_ms_name)
        print('Removed previously present: ' + combined_ms_name)
    
    concat(vis = noisy_ms_name_list, concatvis = combined_ms_name)
    