from datetime import datetime
import os
import pylab as pl
import time

t_start = time.time()

#%%
from simutil import simutil
u = simutil()

#%%

current_dir = os.getcwd()
print(current_dir)

if 'unsync_data' in current_dir:
    file_address = '../input_files/c171a'

else:
    file_address = 'input_files/c171a'


# file_address = 'input_files/mb007'

print(file_address)

#%%
#Specifying values of SEFD & diameter

SEFD_dic = {'EFLSBERG' : 1000,
            'EB_RDBE' : 1000,   
            'ONSALA60' : 5102,
            'YEBES40M' : 1667,
            'METSAHOV' : 17647,
            'KVNYS' : 3226,
            'KVNUS' : 3226,
            'KVNTN' : 3226,     
            'GBT_VLBA' : 137,
            'VLBA_NL' : 2500,           
            'VLBA_BR' : 2500,
            'VLBA_LA' : 2500,
            'VLBA_FD' : 2500,
            'VLBA_PT' : 2500,
            'VLBA_OV' : 2500,
            'VLBA_KP' : 2500,
            'PICOVEL' : 654,
            'VLBA_MK' : 2500,
            'ALMA' : 68,
            'LMT' : 1714,
            'PdBure' : 818}


# print(type(SEFD_dic.values()))

dia_dic = {'EFLSBERG' : 80,
            'EB_RDBE' : 80,       
            'ONSALA60' : 20,
            'YEBES40M' : 40,
            'METSAHOV' : 14,
            'KVNYS' : 21,
            'KVNUS' : 21,
            'KVNTN' : 21,            
            'GBT_VLBA' : 100,
            'VLBA_NL' : 25,           
            'VLBA_BR' : 25,
            'VLBA_LA' : 25,
            'VLBA_FD' : 25,
            'VLBA_PT' : 25,
            'VLBA_OV' : 25,
            'VLBA_KP' : 25,
            'PICOVEL' : 30,
            'VLBA_MK' : 25,
            'ALMA' : 12,
            'LMT' : 50,
            'PdBure' : 33.2}

station_full_name_dic = {'EFLSBERG' : 'Ef',
            'EB_RDBE' : 'Eb',   
            'ONSALA60' : 'On',
            'YEBES40M' : 'Ys',
            'METSAHOV' : 'Mh',
            'KVNYS' : 'Ky',
            'KVNUS' : 'Ku',
            'KVNTN' : 'Kt',     
            'GBT_VLBA' : 'Gb',
            'VLBA_NL' : 'Nl',           
            'VLBA_BR' : 'Br',
            'VLBA_LA' : 'La',
            'VLBA_FD' : 'Fd',
            'VLBA_PT' : 'Pt',
            'VLBA_OV' : 'Ov',
            'VLBA_KP' : 'Kp',
            'PICOVEL' : 'PV',
            'VLBA_MK' : 'Mk',
            'ALMA' : 'Aa',
            'LMT' : 'Lm',             #
            'PdBure' : 'PB'}          #

station_abbrev_dic = {'Ef' : 'EFLSBERG', 
            'Eb' : 'EB_RDBE' ,   
            'On' : 'ONSALA60',
            'Ys' : 'YEBES40M',
            'Mh' : 'METSAHOV',
            'Ky' : 'KVNYS',
            'Ku' : 'KVNUS',
            'Kt' : 'KVNTN',     
            'Gb' : 'GBT_VLBA',
            'Nl' : 'VLBA_NL',           
            'Br' : 'VLBA_BR',
            'La' : 'VLBA_LA',
            'Fd' : 'VLBA_FD',
            'Pt' : 'VLBA_PT',
            'Ov' : 'VLBA_OV',
            'Kp' : 'VLBA_KP',
            'PV' : 'PICOVEL',
            'Mk' : 'VLBA_MK',
            'Aa' : 'ALMA',
            'Lm' : 'LMT' ,             #
            'PB': 'PdBure'}          #

# %%
# reading stations

station = ['Station', 'Code', 'MJD0', 'Adjusted positions [X]', 
           'Adjusted positions [Y]', 'Adjusted positions [Z]']
print(station)

x_adj_dic = {}
y_adj_dic = {}
z_adj_dic = {}
station_names = []

line_counter = 0
with open(file_address + '.sum', 'r') as sum_file:
    lines_list = sum_file.readlines()
    
    for line in lines_list:
        line_counter += 1
        if line[0:24] == '   Plate tectonic motion':
            # print(line)
            table_start_line_number = line_counter + 2

        if line[0:16] == 'BASELINE LENGTHS':
            # print(line)
            table_end_line_number = line_counter - 4 

    for line in lines_list[table_start_line_number:table_end_line_number]:
        station_name = line[3:12].rstrip()
        station = [station_name]
        station_names.append(station_name)
        
        station_code = line[13:15]
        station.append(station_code)
        
        MJDO = line[42:47]
        station.append(MJDO)
        
        x_adj_dic[station_code] = float(line[49:61].lstrip())
        # x_adj.append(float(line[49:61].lstrip()))
        station.append(line[49:61].lstrip())
        
        y_adj_dic[station_code] = float(line[62:74].lstrip())
        # y_adj.append(float(line[62:74].lstrip()))
        station.append(line[62:74].lstrip())
        
        z_adj_dic[station_code] = float(line[76:-1])
        # z_adj.append(float(line[76:-1]))
        station.append(line[76:-1])
        
        print(station)
        
print(x_adj_dic)
print(type(x_adj_dic.values()))
print(x_adj_dic.values())        

# %%        

SEFD = []
dia = []

for name in station_names:
    SEFD.append(SEFD_dic[name])
    dia.append(dia_dic[name])

print(station_names)
print(dia) 
print(SEFD)     

#%%
#Calculating reference position of the antenna configuration

cofa_x = pl.average(x_adj_dic.values())
cofa_y = pl.average(x_adj_dic.values())
cofa_z = pl.average(x_adj_dic.values())
cofa_lat,cofa_lon,cofa_alt = u.xyz2long(cofa_x, cofa_y, cofa_z, 'WGS84')
pos_obs = me.position("WGS84",qa.quantity(cofa_lon, "rad"), qa.quantity(cofa_lat, "rad"), qa.quantity(cofa_alt, "m"))




# %%
# reading SEFD

# SEFD = []
# station_name = []
# dia = []
# with open(file_address[:-5] + 'gmva.cfg', 'r') as gmva_cfg_file:
#     for line in gmva_cfg_file.readlines():
#         line_parts=line.split()
        
#         if (line.startswith('#') == False):
#             station_name.append(line_parts[4])
#             SEFD.append(float(line_parts[5]))
#             dia.append(float(line_parts[3]))

# print(station_name)
# print(SEFD)
# print(dia)



#%%

# for i, station_i in enumerate(stations_list):
#     for j, station_j in enumerate(station_name):
#         if (station_i[0] == station_j.upper()):
#             stations_list[i].insert(1, SEFD[j])

# for station in stations_list:
#     print(station)


#%%
#reading project name and integration time

with open(file_address + '.vex.obs', 'r') as vex_file:

    for line in vex_file.readlines():

        if 'exper_name' in line:
            project_name = line.split('= ')[1][:-2]
            print(project_name)
 
        if 'integr_time' in line:
            # print(line)
            integration_time = line.split(':')[1].lstrip().rstrip()
            print(integration_time)
        
# %%
# reading targets

sources = {}

with open(file_address + '.vex.obs', 'r') as vex_file:

    for line in vex_file.readlines():
 
        if 'source_name' in line:
            # print(line)
            source_name = line.split('= ')[1][:-2]
            direction = [source_name]
            # source = {'name': line.split('= ')[1][:-2]}
            # print(direction)

        if '*' not in line and 'ra =' in line:
            # print(line)
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
            
            # sources.append(source)


# %%
# reading scans

scan_stations_abbrev = []
scan_stations_int_time = []
scans = []

with open(file_address + '.vex.obs', 'r') as vex_file:

    for line in vex_file.readlines():

        if 'scan No' in line:
            scan_no = line[7:11]
            scan = {'scan_no' : scan_no}
            # scan = [scan_no]
                        
        if 'REFERENCE_POINTING_DETERMINE' in line:
            intent = 'REFERENCE_POINTING_DETERMINE'
            scan['intent'] = intent
            # scan.append(intent)
            
        if 'start=' in line and 'exper_nominal' not in line:

            line_parts = line.split('=')

            start_time = line_parts[1][:-7]
            start_time = datetime.strptime(start_time, '%Yy%jd%Hh%Mm%S')
            start_time = start_time.strftime('%Y/%m/%d/%H:%M:%S')
            
            scan['start_time'] = start_time
            # scan.append(start_time)
            

            mode = line_parts[2][:8]
            if ';' in mode:
                mode = mode[:-1]
            scan['mode'] = mode
            # scan.append(mode)

            source_name = line_parts[3][:-2]
            scan['source_name'] = source_name
            # scan.append(source_name)

        if 'station=' in line and 'sec' in line:

            line_parts = line.split('=')

            station_abbrev = line_parts[1][:2]
            scan_stations_abbrev.append(station_abbrev)

            line_parts = line.split(':')

            # print(line_parts)

            int_time = line_parts[2][2:-4]
            scan_stations_int_time.append(int_time)

        if 'endscan' in line:
            
            scan['participating_stations'] = scan_stations_abbrev
            scan['integration_time'] = scan_stations_int_time
            
            if 'REFERENCE_POINTING_DETERMINE' not in scan.values():
                print('\n')
                print(scan)
                
                scans.append(scan)
                
            scan_stations_abbrev = []
            scan_stations_int_time = []

            # break

# %%

print(len(scans))

count = 0
for scan in scans:
    if scan['mode'] == '3mm_RDBE':
        count += 1

print(count)
# %%
# reading modes

mode_freq = []
freq_n_channels = []
freq_bandwidth = []
freq_stations = []

modes = []

with open(file_address + '.vex.obs', 'r') as vex_file:
    vex_file_lines = vex_file.readlines()
    for i, line in enumerate(vex_file_lines):

        if 'ref $PROCEDURES' in line:

            mode_name = vex_file_lines[i - 1][4:-2]
            mode = {'mode_name': mode_name}
            mode['stations'] = []

            line_parts = line.split('=')
            procedures = line_parts[1][1:-2]
            mode['procedures'] = procedures

        if 'ref $FREQ' in line:

            line_parts_equal = line.split('=')
            mode_freq.append(line_parts_equal[1][1:12])
            
            n_channels = line_parts_equal[1][12:14]
            if 'x' in n_channels:
                n_channels = int(n_channels[:-1])
            freq_n_channels.append(int(n_channels))
            
            
            line_parts_colon = line.split(':')
            
            band_width = line_parts_colon[0][-5:]
            if '#' in band_width:
                line_parts_hash = line.split('#')
                band_width = line_parts_hash[0][-5:]
                
            freq_bandwidth.append(band_width)
            
            for part in line_parts_colon[1:]:
                if '\n' in part:
                    part = part[:-2]
                freq_stations.append(part)
                
            mode['stations'].append(freq_stations)
            freq_stations = []

        if 'ref $PASS_ORDER' in line:
            mode['freq'] = mode_freq
            mode['n_channels'] = freq_n_channels
            mode['freq_resolution'] = freq_bandwidth
            # mode['stations'] = freq_stations
            
            modes.append(mode)
            print('\n')
            print(mode)
            
            mode_freq = []
            freq_n_channels = []
            freq_bandwidth = []
            freq_stations = []

# %%
# reading frequency setup

with open(file_address + '.vex.obs', 'r') as vex_file:
    vex_file_lines = vex_file.readlines()
    for i, line in enumerate(vex_file_lines):

        if ('* mode =' in line and 'stations =' in line and 
                'sample_rate' in vex_file_lines[i + 1]):

            freq_name = vex_file_lines[i - 1][4:-2]
            freq_setup = [freq_name]

            line_parts = line.split('=')

            mode_no = line_parts[1][1:3].lstrip()
            freq_setup.append(mode_no)

            stations = line_parts[2][:-1]
            freq_setup.append(stations)

        if 'chan_def' in line:

            line_parts = line.split(':')

            ch_no = line_parts[4][2:-1]
            freq_setup.append(ch_no)

            cent_freq = line_parts[1][1:-1]
            freq_setup.append(cent_freq)
            
            side_band = line_parts[2][1:-1].lstrip()
            freq_setup.append(side_band)

            bandwidth = line_parts[3][2:-1].lstrip()
            freq_setup.append(bandwidth)

            pol = line_parts[-1][-4:-1]
            freq_setup.append(pol)

        if 'enddef' in line and 'chan_def' in vex_file_lines[i - 1]:

            print(freq_setup)
            print('\n')
            
            
  # %% 
# participating_stations_names = []

# x_adj_scan = []
# y_adj_scan = []
# z_adj_scan = []
# SEFD_scan = []
# dia_scan = []
# for station_abbrev in scan['participating_stations']:
#     if station_abbrev in station_abbrev_dic.keys():
#         participating_stations_names.append(station_abbrev_dic[station_abbrev])
    
#     x_adj_scan.append(x_adj_dic[station_abbrev])
#     y_adj_scan.append(y_adj_dic[station_abbrev])
#     z_adj_scan.append(z_adj_dic[station_abbrev])
    
    
    


# for name in participating_stations_names:
#     SEFD_scan.append(SEFD_dic[name])
#     dia_scan.append(dia_dic[name])

# print(dia_scan) 
# print(SEFD_scan)    
    
# print(participating_stations_names)
# # scan['n_channels']
# print(scan['mode'])
# # print(mode['mode_name' = scan['mode']])

# scan['freq'] = []
# scan['freq_resolution'] = []
# scan['n_channels'] = []

# for mode in modes:
#     if mode['mode_name'] == scan['mode']:
#         for ind, station_set in enumerate(mode['stations']):
#             for station_i in station_set:
#                 for station_j in scan['participating_stations']:
#                     if station_i == station_j:
#                         if mode['freq'][ind] not in scan['freq']:
#                             scan['freq'].append(mode['freq'][ind])
                            
#                         if mode['freq_resolution'][ind] not in scan['freq_resolution']:
#                             scan['freq_resolution'].append(mode['freq_resolution'][ind])
                            
#                         if mode['n_channels'][ind] not in scan['n_channels']:
#                             scan['n_channels'].append(mode['n_channels'][ind])
            
            
#             # scan_freq = mode['freq']
# print(scan['freq'])

# # scan['freq'] = f"{float(scan['freq'][0][:-3])/1e3}GHz"

# print("{}GHz".format(float(scan['freq'][0][:-3])/1e3))
# print(scan['freq_resolution'])
# print(scan['n_channels'])

# scan['delta_freq'] = str(float(scan['freq_resolution'][0][:-3]) * scan['n_channels'][0]) + ' MHz'
# print(scan['delta_freq'])

  # %% 


# print(scans[0]['start_time'])
# print(type(scans[0]['start_time']))



# my_date = datetime.strptime(scans[0]['start_time'], '%Yy%jd%Hh%Mm%S')

# my_date = my_date.strftime('%Y/%m/%d/%H:%M:%S')
# print(my_date)
  # %% 

#Calculating reference position of the antenna configuration

# cofa_x = pl.average(x_adj_scan)
# cofa_y = pl.average(y_adj_scan)
# cofa_z = pl.average(z_adj_scan)
# cofa_lat,cofa_lon,cofa_alt = u.xyz2long(cofa_x, cofa_y, cofa_z, 'WGS84')
# pos_obs = me.position("WGS84",qa.quantity(cofa_lon, "rad"), qa.quantity(cofa_lat, "rad"), qa.quantity(cofa_alt, "m"))




# print([x_adj[station_abbrev] for station_abbrev in scan['participating_stations']])

# %%
# wrapper

# name of the array ('GMVA', 'EHT' or a user-defined name)
# array = 'GMVA'
array = 'VLA'


# %%
# simulation code

# os.system('cd unsync_data')


os.system('rm -rf ' + project_name + '.ms')

# Create new MS.
sm.open(project_name+ '-model.ms')
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
   

for scan in scans[:10]:
    
    print('\n')
    print(scan['scan_no'])
    
    #Initialize a source
    sm.setfield(sourcename = scan['source_name'], 
                sourcedirection = sources[scan['source_name']])
      

    
    #setLimits, setfeed, setauto can be out of loop?
    
    #Set integration time and reference time for observation
    sm.settimes(integrationtime = integration_time, 
                usehourangle = False, 
                referencetime = me.epoch('UTC', scan['start_time']))   #start times in UTC?

	#Initialise observation
    sm.observe(scan['source_name'], scans[0]['mode'], 
                starttime = '0 s', 
                stoptime = scan['integration_time'][0] + ' s')
    
    
    
    
sm.close()

print((time.time() - t_start)/60)


# %%

print("Creating the model image with the filename: {} ."
      .format('.modelimage.im'))

# %%



#Create a model image with point source
for source_name, source_direction in sources.items():
    
    print("Creating the model image with the filename: {} ."
          .format(source_name + '.modelimage.im'))
    cl.done()
    
    input_model = 'point' #A string that could be either point, 
                          #Gaussian, disk, or limbdarkeneddisk
                          
    cl.addcomponent(dir = source_direction, flux = 1.0, fluxunit = 'Jy', 
                    freq = modes[0]['freq'][0], shape = input_model)
    
    ia.fromshape(outfile = source_name+'.modelimage.im', shape = [512,512,1,1], 
                 overwrite = True)
    cs = ia.coordsys()
    cs.setunits(['rad','rad','','Hz'])
    
    pix_res = '10 arcsec'   
    cell_rad = qa.convert(pix_res,"rad")['value']
    print(cell_rad)
    cs.setincrement([-cell_rad, cell_rad], 'direction')
    cs.setreferencevalue([qa.convert(source_direction[1],'rad')['value'],
                          qa.convert(source_direction[2],'rad')['value']], 
                         type = "direction")
    
    cs.setreferencevalue(modes[0]['freq'][0], 'spectral')
    cs.setincrement(qa.convert(modes[0]['freq_resolution'][0], 'GHz'), 'spectral')
    ia.setcoordsys(cs.torecord())
    ia.setbrightnessunit("Jy/pixel")
    ia.modify(cl.torecord(), subtract=False)
    exportfits(imagename=project_name+'.modelimage.im', 
               fitsimage=project_name+'.modelimage.fits', overwrite=True)
    
    ia.close()
    cl.done()


# %%


# import winsound

# frequency = 2500  # Set Frequency To 2500 Hertz
# duration = 1000  # Set Duration To 1000 ms == 1 second
# winsound.Beep(frequency, duration)