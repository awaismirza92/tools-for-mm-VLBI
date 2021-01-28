import pandas as pd

file_address = 'input_files/c171a'
# file_address = 'input_files/mb007'

#%%
#reading stations


line_counter = 0
with open(file_address + '.sum', 'r') as sum_file:
    for line in sum_file.readlines():
        line_counter += 1
        if line[0:24]== '   Plate tectonic motion':
            # print(line)
            skip_lines_start = line_counter + 2
            
        if line[0:16]== 'BASELINE LENGTHS':
            print(line)
            table_end_line_number = line_counter - 4

skip_lines_end = line_counter - table_end_line_number
            
# print(skip_lines_start, table_end_line_number, skip_lines_end)



            
col_names = ['Station', 'Code', 'Station motions [X]', 'Station motions [Y]', 
              'Station motions [Z]', 'MJD0', 'Adjusted positions [X]', 'Adjusted positions [Y]', 
              'Adjusted positions [Z]']            

sum_file = pd.read_csv('input_files/c171a.sum', skiprows = skip_lines_start, 
                        skipfooter = skip_lines_end, 
                        engine = 'python', names = col_names, sep='\s+', index_col=0)

# print(sum_file.shape)
print(sum_file)



sum_adjusted_pos = sum_file[['Adjusted positions [X]', 'Adjusted positions [Y]', 
              'Adjusted positions [Z]']]

print(sum_adjusted_pos)





#%%
#reading targets

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
            direction.append(coord_parts[2][:-18].replace("'", 'm').replace('"', 's'))
            print(direction)
            
            print('\n')

            

#%%
#reading scans

with open(file_address + '.vex.obs', 'r') as vex_file:
    
    for line in vex_file.readlines():
        
        if 'scan No' in line:
            scan_no = line[7:11]
            scan = [scan_no]
            print('\n')    
            
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
            
            print(scan) 
            
            # break
    
    
    
#%%
#reading modes    
    
with open(file_address + '.vex.obs', 'r') as vex_file:
    vex_file_lines = vex_file.readlines()
    for i, line in enumerate(vex_file_lines):
        
        if 'ref $PROCEDURES' in line:
        
            mode_name = vex_file_lines[i-1][4:-2]
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
   
            
#%%           
#reading frequency setup

with open(file_address + '.vex.obs', 'r') as vex_file:
    vex_file_lines = vex_file.readlines()
    for i, line in enumerate(vex_file_lines):
        
        if '* mode =' in line and 'stations =' in line and 'sample_rate' in vex_file_lines[i+1] :
            
            freq_name = vex_file_lines[i-1][4:-2]
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
            
            bandwidth = line_parts[3][2:-1].lstrip()
            freq_setup.append(bandwidth)
            
            pol = line_parts[-1][-4:-1]
            freq_setup.append(pol)
            
                        
        if 'enddef' in line and 'chan_def' in vex_file_lines[i-1]:
            
            print(freq_setup)
            print('\n')
            
#%%        
#wrapper

# name of the array ('GMVA', 'EHT' or a user-defined name) 
array='GMVA'


#%%
#simulation code

#Create new MS.
sm.open(projectname)
sm.setconfig(telescopename=array, x=xx, y=yy, z=zz, dishdiameter=diam, mount='alt-az', antname=pads, coordsystem='global', referencelocation=pos)
