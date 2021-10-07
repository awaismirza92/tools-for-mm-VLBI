from datetime import datetime
import pylab as pl
import time
import numpy as np

t_start = time.time()

#%%
from simutil import simutil
u = simutil()


#%%
#Specifying values of SEFD & diameter
print('Specifying SEFDs, diameters & abreviations for stations. ')

SEFD_dic_name = {'EFLSBERG' : 1000,
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

SEFD_dic = {'Eb' : 1000, 
            'On' : 5102,
            'Ys' : 1667,
            'Mh' : 17647,
            'Ky' : 3226,
            'Ku' : 3226,
            'Kt' : 3226,     
            'Gb' : 137,
            'Nl' : 2500,           
            'Br' : 2500,
            'La' : 2500,
            'Fd' : 2500,
            'Pt' : 2500,
            'Ov' : 2500,
            'Kp' : 2500,
            'PV' : 654,
            'Mk' : 2500,
            'Aa' : 68,
            'Lm' : 1714,
            'PB' : 818}


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

print('Done. \n')

# %%
# reading stations
print('Reading names, codes, adjusted coordinates of stations.')

station = ['Station', 'Code', 'MJD0', 'Adjusted positions [X]', 
           'Adjusted positions [Y]', 'Adjusted positions [Z]']


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
            table_start_line_number = line_counter + 2

        if line[0:16] == 'BASELINE LENGTHS':
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
        station.append(line[49:61].lstrip())
        
        y_adj_dic[station_code] = float(line[62:74].lstrip())
        station.append(line[62:74].lstrip())
        
        z_adj_dic[station_code] = float(line[76:-1])
        station.append(line[76:-1])

print('Done. \n')


# %%        
print('Extracting SEFDs & diameters of participating stations.')

SEFD = []
dia = []

for name in station_names:
    SEFD.append(SEFD_dic_name[name])
    dia.append(dia_dic[name])

print('Done. \n')

#%%
print('Calculating reference position of the antenna configuration.')
#Calculating reference position of the antenna configuration

cofa_x = pl.average(x_adj_dic.values())
cofa_y = pl.average(x_adj_dic.values())
cofa_z = pl.average(x_adj_dic.values())
cofa_lat,cofa_lon,cofa_alt = u.xyz2long(cofa_x, cofa_y, cofa_z, 'WGS84')
pos_obs = me.position("WGS84",qa.quantity(cofa_lon, "rad"), qa.quantity(cofa_lat, "rad"), qa.quantity(cofa_alt, "m"))

print('Done. \n')



#%%
#reading project name and integration time

with open(file_address + '.vex.obs', 'r') as vex_file:

    for line in vex_file.readlines():

        if 'exper_name' in line:
            project_name = line.split('= ')[1][:-2]

        if 'integr_time' in line:
            integration_time = line.split(':')[1].lstrip().rstrip()

# %%
# reading targets
print('Reading science targets. ')
sources = {}

with open(file_address + '.vex.obs', 'r') as vex_file:

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
            
            sources[source_name] = [epoch, RA, Dec]
            
print('Done. \n')


# %%
# reading scans
print('Reading all scans of the observing session. ')

scan_stations_abbrev = []
scan_stations_obs_time = []
scans = []

with open(file_address + '.vex.obs', 'r') as vex_file:

    for line in vex_file.readlines():

        if 'scan No' in line:
            scan_no = line[7:11]
            scan = {'scan_no' : scan_no}

        if 'REFERENCE_POINTING_DETERMINE' in line:
            intent = 'REFERENCE_POINTING_DETERMINE'
            scan['intent'] = intent

        if 'start=' in line and 'exper_nominal' not in line:
            line_parts = line.split('=')

            start_time = line_parts[1][:-7]
            start_time = datetime.strptime(start_time, '%Yy%jd%Hh%Mm%S')
            
            start_date = start_time.strftime('%d-%b-%Y')
            
            start_time = start_time.strftime('%Y/%m/%d/%H:%M:%S')            
            scan['start_time'] = start_time

            mode = line_parts[2][:8]
            if ';' in mode:
                mode = mode[:-1]
            scan['mode'] = mode

            source_name = line_parts[3][:-2]
            scan['source_name'] = source_name

        if 'station=' in line and 'sec' in line:
            line_parts = line.split('=')

            station_abbrev = line_parts[1][:2]
            scan_stations_abbrev.append(station_abbrev)

            line_parts = line.split(':')

            obs_time = int(line_parts[2][2:-4])
            scan_stations_obs_time.append(obs_time)

        if 'endscan' in line:
            
            scan['participating_stations'] = scan_stations_abbrev
            scan['observation_time'] = scan_stations_obs_time
            
            if 'REFERENCE_POINTING_DETERMINE' not in scan.values():
                scans.append(scan)
                
            scan_stations_abbrev = []
            scan_stations_obs_time = []

print('Done. \n')

# %%
# reading modes
print('Reading scanning modes.')

mode_freq = []
freq_n_channels = []
freq_bandwidth = []
freq_stations = []

modes = []

with open(file_address + '.vex.difx', 'r') as vex_file:
    vex_file_lines = vex_file.readlines()
    for i, line in enumerate(vex_file_lines):

        if 'ref $PROCEDURES' in line:
            
            if 'def' in vex_file_lines[i - 1]:
                mode_name = vex_file_lines[i - 1][4:-2]
                
            if 'def' in vex_file_lines[i - 2]:
                mode_name = vex_file_lines[i - 2][4:-2]
            
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

            modes.append(mode)

            mode_freq = []
            freq_n_channels = []
            freq_bandwidth = []
            freq_stations = []

print('Done. \n')

# %%
# reading frequency setup
print('Reading frequency setups.')

freq_defs = []

with open(file_address + '.vex.difx', 'r') as vex_file:
    vex_file_lines = vex_file.readlines()
    for i, line in enumerate(vex_file_lines):

        if ('* mode =' in line and 'stations =' in line and 
                'sample_rate' in vex_file_lines[i + 1]):

            freq_name = vex_file_lines[i - 1][4:-2]
            freq_setup = [freq_name]
            freq_def = {'name': freq_name}

            line_parts = line.split('=')

            mode_no = line_parts[1][1:3].lstrip()
            freq_setup.append(mode_no)
            freq_def['mode_no'] = mode_no

            stations = line_parts[2][:-1]
            stations = stations.split(':')
            freq_setup.append(stations)
            freq_def['stations'] = stations
            
            ch_nos = []
            cent_freqs = []
            side_bands = []
            bandwidths = []
            pols = []

        if 'chan_def' in line:

            line_parts = line.split(':')

            ch_no = line_parts[4][2:-1]
            freq_setup.append(ch_no)
            ch_nos.append(ch_no)

            cent_freq = line_parts[1][1:-1]
            freq_setup.append(cent_freq)
            cent_freqs.append(cent_freq)
            
            side_band = line_parts[2][1:-1].lstrip()
            freq_setup.append(side_band)
            side_bands.append(side_band)

            bandwidth = line_parts[3][2:-1].lstrip()
            freq_setup.append(bandwidth)
            bandwidths.append(bandwidth)

            pol = line_parts[-1][-4:-1]
            freq_setup.append(pol)
            pols.append(pol)

        if 'enddef' in line and 'chan_def' in vex_file_lines[i - 1]:

            freq_def['ch_nos'] = ch_nos
            freq_def['cent_freqs'] = cent_freqs
            freq_def['side_bands'] = side_bands
            freq_def['bandwidths'] = bandwidths
            freq_def['pols'] = pols
            
            
            freq_def['no_lcp_ch'] = freq_def['pols'].count('Lcp') 
            freq_def['no_rcp_ch'] = freq_def['pols'].count('Rcp')
            
            # for item in freq_def.items():
            #     print(item)
            
            freq_defs.append(freq_def)
            
print('Done. \n')

#%%
#reading number of channels
print('Reading number of channels.')

with open(file_address + '.vex.difx', 'r') as vex_file:

    for line in vex_file.readlines():

        if 'number_channels' in line:
            nchannels = int(line.split(':')[1])
            # print('nchannels: ', nchannels)


#Number of polarizations has been hard coded to be equal to 2. If this parameter
#needs to be made flexible later on, changes may be done in the following lines.

        # if 'dual polarisation' in line:
npol = 2.0 
            # print('npol', npol)
 
print('Done. \n')
