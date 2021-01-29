

file_address = 'input_files/c171a'
file_address = 'input_files/mb007'


# %%
# reading stations

station = ['Station', 'Code', 'MJD0', 'Adjusted positions [X]', 
           'Adjusted positions [Y]', 'Adjusted positions [Z]']
print(station)

stations_list = []
x_adj = []
y_adj = []
z_adj = []

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
        
        station_code = line[13:15]
        station.append(station_code)
        
        MJDO = line[42:47]
        station.append(MJDO)
        
        x_adj = line[49:61].lstrip()
        station.append(x_adj)
        
        y_adj = line[62:74].lstrip()
        station.append(y_adj)
        
        z_adj = line[76:-1]
        station.append(z_adj)
        
        print(station)
        
        stations_list.append(station)

# %%        

print(stations_list[:][0])        
# %%
# reading SEFD

SEFD = []
station_name = []
with open(file_address[:-5] + 'gmva.cfg', 'r') as gmva_cfg_file:
    for line in gmva_cfg_file.readlines():
        line_parts=line.split()
        
        if (line.startswith('#') == False):
            station_name.append(line_parts[4])
            SEFD.append(float(line_parts[5]))

print(station_name)
print(SEFD)

SEFD_dic = {'EFLSBERG' : 1000,
            'EB_RDBE' : 1000,       #?
            'ONSALA60' : 5102,
            'YEBES40M' : 1667,
            'METSAHOV' : 17647,
            'KVNYS' : 3226,
            'KVNUS' : 3226,
            'KVNTN' : 3226,            #?
            'GBT_VLBA' : 137,
            'VLBA_NL' : 2500,           #?
            'VLBA_BR' : 2500,
            'VLBA_LA' : 2500,
            'VLBA_FD' : 2500,
            'VLBA_PT' : 2500,
            'VLBA_OV' : 2500,
            'VLBA_KP' : 2500,
            'PICOVEL' : 654,
            'VLBA_MK' : 2500,
            'ALMA' : 68}


print(type(SEFD_dic.values()))
#%%

for i, station_i in enumerate(stations_list):
    for j, station_j in enumerate(station_name):
        if (station_i[0] == station_j.upper()):
            stations_list[i].insert(1, SEFD[j])

for station in stations_list:
    print(station)
        
# %%
# reading targets

with open(file_address + '.vex.obs', 'r') as vex_file:

    for line in vex_file.readlines():

        if 'exper_name' in line:
            projectname = line.split('= ')[1][:-2]
            # print(projectname)

        if 'source_name' in line:
            # print(line)
            direction = [line.split('= ')[1][:-2]]
            # print(direction)

        if '*' not in line and 'ra =' in line:
            # print(line)
            coord_parts = line.split('= ')
            direction.append(coord_parts[-1][:-2])
            direction.append(coord_parts[1][:-6])
            direction.append(
                coord_parts[2][:-18].replace("'", 'm').replace('"', 's'))
            print(direction)

            print('\n')


# %%
# reading scans

with open(file_address + '.vex.obs', 'r') as vex_file:

    for line in vex_file.readlines():

        if 'scan No' in line:
            scan_no = line[7:11]
            scan = [scan_no]
                        
        if 'REFERENCE_POINTING_DETERMINE' in line:
            intent = 'REFERENCE_POINTING_DETERMINE'
            scan = [intent]
            
        if 'start=' in line and 'exper_nominal' not in line:

            line_parts = line.split('=')

            start_time = line_parts[1][:-7]
            scan.append(start_time)

            mode = line_parts[2][:8]
            scan.append(mode)

            source_name = line_parts[3][:-2]
            scan.append(source_name)

        if 'station=' in line and 'sec' in line:

            line_parts = line.split('=')

            station_abbrev = line_parts[1][:2]
            scan.append(station_abbrev)

            line_parts = line.split(':')

            # print(line_parts)

            int_time = line_parts[2][2:-4]
            scan.append(int_time)

        if 'endscan' in line:
            if 'REFERENCE_POINTING_DETERMINE' not in scan:
                print('\n')
                print(scan)

            # break


# %%
# reading modes

with open(file_address + '.vex.obs', 'r') as vex_file:
    vex_file_lines = vex_file.readlines()
    for i, line in enumerate(vex_file_lines):

        if 'ref $PROCEDURES' in line:

            mode_name = vex_file_lines[i - 1][4:-2]
            mode = [mode_name]

            line_parts = line.split('=')
            procedures = line_parts[1][1:-2]
            mode.append(procedures)

        if 'ref $FREQ' in line:

            line_parts = line.split('=')
            freq = line_parts[1][1:-2]
            mode.append(freq)

        if 'ref $PASS_ORDER' in line:

            print(mode)


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
# wrapper

# name of the array ('GMVA', 'EHT' or a user-defined name)
array = 'GMVA'


# %%
# simulation code

# Create new MS.
sm.open(projectname)
sm.setconfig(
    telescopename=array,
    x=xx,
    y=yy,
    z=zz,
    dishdiameter=diam,
    mount='alt-az',
    antname=pads,
    coordsystem='global',
    referencelocation=pos)
