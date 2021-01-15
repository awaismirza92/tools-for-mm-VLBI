import pandas as pd

line_counter = 0
with open('input_files/c171a.sum', 'r') as sum_file:
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
# print(sum_file)

# print(sum_file['Station'])
# print(sum_file.loc[1:5, 'Station':'MJD0'])
# print(sum_file.loc[:, ['Station', 'Adjusted positions [X]', 'Adjusted positions [Y]', 
#              'Adjusted positions [Z]']])

# print(sum_file[['Station', 'Adjusted positions [X]', 'Adjusted positions [Y]', 
#              'Adjusted positions [Z]']])

sum_adjusted_pos = sum_file[['Adjusted positions [X]', 'Adjusted positions [Y]', 
              'Adjusted positions [Z]']]

print(sum_adjusted_pos)


# print(sum_file.loc['PICOVEL'])
# print(sum_file['Code'])


#%%

with open('input_files/c171a.vex.obs', 'r') as vex_file:
    
    for line in vex_file.readlines():
        
        if 'exper_name' in line:
            project_name = line.split('= ')[1][:-2]
            print(project_name)
            
        if 'source_name' in line:
            print(line)
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

with open('input_files/c171a.key', 'r') as key_file:
    for line in key_file.readlines():
        if '!' not in line and 'SOURCE' in line and 'EQUINOX' in line:
            line_parts = line.split('=')
            direction = [line_parts[-1][1:-4], line_parts[2][:-4].rstrip(), line_parts[3][1:-8].rstrip() ]
            print(direction)
            
            print('\n')
            

#%%

with open('input_files/c171a.vex.obs', 'r') as vex_file:
    
    for line in vex_file.readlines():
        
        if 'scan No' in line:
            scan_no = line[7:11]
            scan = [scan_no]
            print('\n')    
            
        if 'start=' in line and 'exper_nominal' not in line:
            
            line_parts = line.split('=')
            
            start_time = line_parts[1][:-7]
            scan.append(start_time)
            
            source_name = line_parts[3][:-2]
            scan.append(source_name)
                        
        if 'station=' in line and 'sec' in line:
            
            line_parts = line.split('=')
            
            station_abbrev = line_parts[1][:2]
            scan.append(station_abbrev)
            
            line_parts = line.split(':')
            
            print(line_parts)
            
            int_time = line_parts[2][2:-4]
            scan.append(int_time)
                        
        if 'endscan' in line:
            
            print(scan) 
            
            break
    
    
