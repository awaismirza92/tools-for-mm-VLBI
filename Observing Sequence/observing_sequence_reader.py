from datetime import datetime
import os
import pylab as pl
import time
import numpy as np

t_start = time.time()

#%%
from simutil import simutil
u = simutil()



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
            
            
 

