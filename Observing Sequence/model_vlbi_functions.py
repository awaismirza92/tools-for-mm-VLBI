#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys


def flux_file_existance_check(sources, start_date, freq = '86E+09'):
    try:
        with open('source_fluxes.txt', 'r') as flux_file:
            print("The file 'source_fluxes.txt' has been found")
            
    except FileNotFoundError:
        print("Error: The file 'source_fluxes.txt' was not found in the current " +
              f"working directory ({os.getcwd()}). \n")
        
        template_file = open("source_fluxes_template.txt", "w")
        template_file.write('Source'.ljust(14) +  
                      'Date'.ljust(14) + 
                      'Freq(Hz)'.ljust(14) + 
                      'Flux(Jy)'.ljust(14) + 
                      '\n')
        for source_name in sources.keys():
            template_file.write(source_name.ljust(14) +  
                              start_date.ljust(14) + 
                              freq.ljust(14) + 
                              '99'.ljust(14) + 
                              '\n')
        template_file.close()
            
        print("Please use the code 'Calibrator_Source_Catalogue_Query_Code.py' to" +
              " obtain 'source_fluxes.txt' file for your sources. Otherwise, please " +
              "make the file manually yourself. A template file " + 
              "'source_fluxes_template.txt' has been generated for your reference. " +
              "Please note that the template has incorrect flux values (99 Jy) in it.")
        print('\nThe code is terminating now.')
        sys.exit()
        
   
        
   
    
def flux_file_completion_check(sources, start_date):
        
    initial_print = True
    
    flux_sources = []
    flux_sources_names = []
    with open('source_fluxes.txt', 'r') as flux_file:
        lines_list = flux_file.readlines()
        
        for line in lines_list[1:]:
            
            line_parts = line.split()
            
            flux_sources_names.append(line_parts[0])
            
            flux_source = {'source_name' : line_parts[0]}
            flux_source['date'] = line_parts[1]
            flux_source['freq'] = line_parts[2]
            flux_source['flux'] = line_parts[3]
            
            flux_sources.append(flux_source)
            
            if initial_print == True:
                print("The following sources have been found in " + 
                      "'source_fluxes.txt': ")
                initial_print = False
                
            print(flux_source)
            
            
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
         print("All the observation sources are present in 'source_fluxes.txt'")