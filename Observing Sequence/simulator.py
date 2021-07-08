import os

#%%

current_dir = os.getcwd()
print(current_dir)

# os.chdir('/mnt/e_drive/MEGA/MSc_Astrophysics/5th_Semester/Job/Codes/reading_files/')

current_dir = os.getcwd()
print(current_dir)

if 'small_ms' or 'per_scan_ms' in current_dir:
    file_address = '../../input_files/c171a'

else:
    file_address = 'input_files/c171a'


# file_address = 'input_files/mb007'

print(file_address)
#%%

# import observing_sequence_reader

execfile('../../observing_sequence_reader.py')
execfile('../../simvlbiv2.3_borrows.py')


# %%
# wrapper

# name of the array ('GMVA', 'EHT' or a user-defined name)
# array = 'GMVA'
array = 'VLA'

# project_name = project_name + '_10_scans'

# ms_name = project_name + '.ms'

# %%
# simulation code

perform_observation = 'Yes'
create_input_models = 'Yes'
corrupt_model_data = 'Yes'
create_noisy_ms = 'Yes'
combine_measurment_sets = 'Yes'


perform_observation = 'No'
create_input_models = 'No'
corrupt_model_data = 'No'
create_noisy_ms = 'No'
# combine_measurment_sets = 'No'
   
if perform_observation == 'Yes':

    # os.system('rm -rf ' + ms_name)
    
    # # Create new MS.
    # sm.open(project_name+ '.ms')
    # sm.setconfig(
    #     telescopename = array,
    #     x = x_adj_dic.values(),
    #     y = y_adj_dic.values(),
    #     z = z_adj_dic.values(),
    #     dishdiameter = dia,
    #     mount = 'alt-az',
    #     antname = station_names,
    #     coordsystem = 'global',
    #     referencelocation = pos_obs)
    
    # #Initialize spectral windows. Here we get to set a name for each band.
    # sm.setspwindow(
    #     spwname = scans[0]['mode'], 
    #     freq = modes[0]['freq'][0],                      #"{}GHz".format(float(scan['freq'][0][:-3])/1e3),    #Converting from MHz to GHz
    #     deltafreq = modes[0]['freq_resolution'][0],      #scan['delta_freq'], 
    #     freqresolution = modes[0]['freq_resolution'][0], #scan['freq_resolution'][0], 
    #     nchannels = modes[0]['n_channels'][0],           #scan['n_channels'][0], 
    #     stokes = 'RR LL')             
    
    
    # elev_limit = 10.0 
    
    # #Set Limits to flag data such as when source is below certain elevation.
    # sm.setlimits(shadowlimit = 0.001, elevationlimit = str(elev_limit)+'deg') 
    
    # sm.setfeed('perfect R L')         
    # sm.setauto(autocorrwt = 0.0)  



    for scan in scans[80:100]:
        
        ms_name = 'scan_' + scan['scan_no'] + '.ms'
        
        os.system('rm -rf ' + ms_name)
    
        # Create new MS.
        sm.open(ms_name)
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
    
        #Initialize spectral windows. Here we get to set a name for each band.
        sm.setspwindow(
        spwname = scans[0]['mode'], 
        freq = modes[0]['freq'][0],                      #"{}GHz".format(float(scan['freq'][0][:-3])/1e3),    #Converting from MHz to GHz
        deltafreq = modes[0]['freq_resolution'][0],      #scan['delta_freq'], 
        freqresolution = modes[0]['freq_resolution'][0], #scan['freq_resolution'][0], 
        nchannels = modes[0]['n_channels'][0],           #scan['n_channels'][0], 
        stokes = 'RR LL')             
    
    
        elev_limit = 10.0 
        
        #Set Limits to flag data such as when source is below certain elevation.
        sm.setlimits(shadowlimit = 0.001, elevationlimit = str(elev_limit)+'deg') 
        
        sm.setfeed('perfect R L')         
        sm.setauto(autocorrwt = 0.0)  
        
        
        
        
        print('\n')
        print(scan['scan_no'])
        
        #Initialize a source
        sm.setfield(sourcename = scan['source_name'], 
                    sourcedirection = sources[scan['source_name']])
          
    
        
        #setLimits, setfeed, setauto can be out of loop?
        
        #Set integration time and reference time for observation
        sm.settimes(integrationtime = integration_time, 
                    usehourangle = False, 
                    referencetime = me.epoch('TAI', scan['start_time']))   #start times in UTC?
    
    	#Initialise observation
        sm.observe(scan['source_name'], scans[0]['mode'], 
                    starttime = '0 s', 
                    stoptime = str(scan['observation_time'][0]) + ' s')
    
    
    
    
    sm.close()
    
    # print((time.time() - t_start)/60)


# %%

# Now we compute the baseline lengths between the different antennae and
    # finally compute the largest baseline length to be used to compute
    # synthesised beam.  [lat,lon,el] (x, y, z) [m] are relative to the center
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
lambd = np.round(300 / qa.quantity(modes[0]['freq'][0], 'MHz')['value'], 4)  # in m
# print(lambd)
# lambd = np.round(3e8 / (float(modes[0]['freq'][0][:-3]) * 1e6), 4)  # in m
beam_max = np.round(np.amax(mylengths),  4)  # Longest baseline in m.

pix_res = lambd / beam_max * 180 / np.pi * 3600
pix_res = '{}'.format(pix_res) + 'arcsec'

print(qa.quantity(modes[0]['freq'][0], 'GHz')['value'], modes[0]['freq'][0], 
      modes[0]['freq'][0][:-3], lambd, beam_max, pix_res)

# %%



   
if create_input_models == 'Yes':


    #Create a model image with point source
    for source_name, source_direction in sources.items():
        
        print("Creating the model image with the filename: {} ."
              .format(source_name + '.modelimage.im'))
        cl.done()
        
        input_model = 'point' #A string that could be either point, 
                              #Gaussian, disk, or limbdarkeneddisk
                              
        cl.addcomponent(dir = source_direction, flux = 8.0, fluxunit = 'Jy', 
                        freq = modes[0]['freq'][0], shape = input_model)
        
        ia.fromshape(outfile = 'model_images/' + source_name+'.modelimage.im', 
                     shape = [512,512,1,1], 
                     overwrite = True)
        cs = ia.coordsys()
        cs.setunits(['rad','rad','','Hz'])
        
        # pix_res = '10 arcsec'   
        cell_rad = qa.convert(pix_res,"rad")['value']
        print(cell_rad)
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



   
if corrupt_model_data == 'Yes':

    # Corrupt model image with voltage pattern, copy it to data column of
    # Measurement Set.
    
    t_start = time.time()
    
    print("Corrupting model data with voltage pattern which is estimated "
          + "separately for each antenna in the array.")
    im.open(ms_name)
    im.selectvis()
    im.defineimage(nx = 512, ny = 512, cellx = pix_res, celly = pix_res, 
                   facets = 2)
    im.setvp(dovp = True, usedefaultvp = True)
    im.setoptions(ftmachine = 'mosaic')
    im.ft(model = 'model_images/' + 'Einstein.im')
    # im.ft(model = 'model_images/' + '3C84' + '.modelimage.im')
    im.close()
    im.open(thems = ms_name, usescratch = True)
    uvsub(vis = ms_name, reverse = True)
    im.close()
    print((time.time() - t_start)/60)
    
    
    
# Create noisy ms.
if create_noisy_ms == 'Yes':
    
    print("Visibilities added to Data column of MS. Now corrupting visibilities "
          + "with previously computed rms noise. New noisy MS will be created.")
    
    nbits = 2.0                             #It's standard to have two bits
    eta_c = (0.88 if nbits == 2 else 0.64)  #Correlation efficient remains 0.88 if nbits = 2
    nbase = len(mylengths)  # Number of baselines
    
    # data_rate = 2048                        
    data_rate = (2.0 * qa.quantity(freq_setup[-2], 'MHz')['value'] 
                    * npol * nbits * 2.0)
        
    
    sigma_array = []
    for scan in scans[80:100]:
        
        print('\n')
        print(scan['scan_no'])        
        # print(scan['participating_stations'])
        
        participating_stations_SEFDs = []
        for station in scan['participating_stations']:
            participating_stations_SEFDs.append(SEFD_dic[station])            
        # print(participating_stations_SEFDs)
        
        
        participating_stations_SEFDs = np.array(participating_stations_SEFDs)       
        SEFD_reciprocal = 1. / np.outer(participating_stations_SEFDs, participating_stations_SEFDs)
        print(SEFD_reciprocal)
        
        
        SEFD_star = 1 / np.sqrt(np.sum(SEFD_reciprocal))
        print('SEFD_star:', SEFD_star)
            
        
        print('Observation time: ', scan['observation_time'][0])
        
        
        noise = ( ((1/eta_c) * SEFD_star) / 
                np.sqrt(data_rate * scan['observation_time'][0] / 2))
        print('Noise per scan:', noise)
        
        
        integration_time_qa = qa.quantity(integration_time, 's')
        # print(integration_time_qa['value'])
        
        
        sigma = noise * np.sqrt(nchannels * npol * nbase * 
                                scan['observation_time'][0] / 
                                integration_time_qa['value'])
        print('Sigma per scan: ', sigma)
        
        sigma_array.append(sigma)
        
        
        
        
        ms_name = 'scan_' + scan['scan_no'] + '.ms'
        
        
        
        sigma = qa.quantity(sigma, 'Jy')
        
        noisy_ms_name = ms_name[:-3] + '.noisy.ms'
        
        os.system('cp -r ' + ms_name + ' ' + noisy_ms_name)
        sm.openfromms(noisy_ms_name)
        # sm.setfield()
        sm.setnoise(mode = 'simplenoise', simplenoise = sigma)
        sm.corrupt()
        
        
        # To invoke uv-domain primary beam convolution for heterogeneous arrays we
        # set the fourier transform machine to mosaic.
        
        sm.setoptions(ftmachine="mosaic")
        sm.done()
        sm.close()
        

    print(sigma_array)
        
    print(freq_setup[-2][:-3])
    print(data_rate)
    
    
        
    # sigma = qa.quantity(sigma, 'Jy')
    
    
    # os.system('cp -r ' + ms_name + ' ' + ms_name + '.noisy.ms')
    # sm.openfromms(ms_name + '.noisy.ms')
    # sm.setfield()
    # sm.setnoise(mode = 'simplenoise', simplenoise = sigma)
    # sm.corrupt()


    # To invoke uv-domain primary beam convolution for heterogeneous arrays we
    # set the fourier transform machine to mosaic.
    
    # sm.setoptions(ftmachine="mosaic")
    # sm.done()
    # sm.close()
    
    
    print(ms_name
    + ".noisy.ms has been created which contains the visibilities of "
    + ms_name + ".ms corrupted with thermal noise.")


print('\n')
    
if combine_measurment_sets == 'Yes':
    
    noisy_ms_name_list = []
    for scan in scans[80:100]:
        ms_name = 'scan_' + scan['scan_no'] + '.ms'
        noisy_ms_name = ms_name[:-3] + '.noisy.ms'
        
        noisy_ms_name_list.append(noisy_ms_name)
        
    print(noisy_ms_name_list)
    
    
    concat(vis = noisy_ms_name_list, concatvis = 'all_scans_combined.ms')
          
    
print('\n')
print('perform_observation: ', perform_observation)
print('create_input_models: ', create_input_models)
print('corrupt_model_data: ', corrupt_model_data)
print('create_noisy_ms: ', create_noisy_ms)


# %%


# import winsound

# frequency = 2500  # Set Frequency To 2500 Hertz
# duration = 1000  # Set Duration To 1000 ms == 1 second
# winsound.Beep(frequency, duration)