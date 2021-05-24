#!/usr/bin/python

import numpy as np
import itertools
import os
import sys
import time

#### Updated 24/05/2021
#### Script for parameter testing on a difx installation. 
#### usage: ./difxParamTest <v2d_file> <parameter_input_file> scan_start_integer scan_end_integer 
#### WARNING: Currently has no in-built file clean up or time-out functionality - so please make sensible guesses for parameter values.

cwd = os.getcwd()

def paramInput(input_filename):
    with open(input_filename, 'r') as input_file:
        input_lines = input_file.readlines()
    # Determine what sections exist in parameter file
    section = []
    section_index = []
    for i in range(0,len(input_lines)):
        if '--' in input_lines[i]:
            section.append(input_lines[i].rstrip().strip('-'))
            section_index.append(i)
    # Extract parameters for testing
    parameter_list = []
    for j in range(0,len(section)):
        parameter_temp = []
        index = section_index[j] + 1 # starting position
        while (index not in section_index) and (index < len(input_lines)): # check index is valid for our input file and is not a section title
            parameter_temp.append(input_lines[index].rstrip().split('=')) # append parameters, split about the equal signS
            index += 1
        parameter_list.append(parameter_temp)
    # Convert values to floats + extract them to determine all possible combinations
    all_values = []
    all_params = []
    section_ref = []
    for i in range(0, len(parameter_list)): # iterate through sections
        if len(parameter_list[0]) != 0:
            for j in range(0, len(parameter_list[i])): # iterate through params in one section
                values = np.array(parameter_list[i][j][1].split(',')).astype(float)
                all_values.append(values) # append all possible values - needed for determining combinations
                all_params.append(parameter_list[i][j][0]) # append parameter name for reference
                section_ref.append(section[i])
    combinations = list(itertools.product(*all_values))     
    # return combinations, parameters and section reference list
    return combinations, all_params, section_ref

def removeLines(lines_list, params):
    remove_indices = []
    for line in lines_list:
        if any(param in line for param in params):
            remove_indices.append(lines_list.index(line))
    for index in sorted(remove_indices, reverse=True):
        del lines_list[index]    
    return lines_list
    
def generateV2D(reference_v2d, combin, params, sec_refs): # reference file name, combination list, parameter names, section references
    print("Reading in reference v2d file: " + str(reference_v2d))
    with open(reference_v2d, 'r') as v2d_file:
        v2d_lines = v2d_file.readlines()
    print("\nCreating " + str(len(combin)) + " individual v2d files.")
    with open(reference_v2d.split('.')[0] + '_reference_table.out', 'w') as ref_file: # open reference table file to append filenames and combinations
        print('# ' + " ".join(params), file=ref_file)
        for combination in combin:
            print(reference_v2d.split('.')[0] + 'Test' + str(combin.index(combination)+1) + ' ' + str(combination), file=ref_file) 
            v2d_lines = removeLines(v2d_lines, params)
            for i in range(0, len(combination)):
                param_string = params[i] + '=' + str(combination[i])
                if sec_refs[i] == 'GLOBAL':
                    v2d_lines.insert(0, param_string)
                else:
                    for line in v2d_lines:
                        if sec_refs[i] in line:
                            section_start = v2d_lines.index(line)
                            v2d_lines.insert(section_start+2, param_string) # +2 should allways work
            with open(reference_v2d.split('.')[0] + 'Test' + str(combin.index(combination)+1) + '.v2d', 'w') as out_file: # write out the v2d file
                for line in v2d_lines:
                    print(line.rstrip(), file=out_file)
    print('\nDone!')
    
def corrScan(scan_tag, num_processes):
    # Setup start time
    start_time = time.time()
    # Execute correlation
    execution_string = "mpirun -machinefile "+ scan_tag + ".machines -mca rmaps seq -np "+str(num_processes)+" mpifxcorr " + scan_tag +".input"
    os.system(execution_string)
    # Calculate time delta
    end_time= time.time()
    time_delta = end_time - start_time
    return time_delta
    
def padCheck(file_tag, cwd):
    dir_list = os.listdir(cwd)
    pad_length = 0
    for file in dir_list:
        if file_tag + '_' in file:
            pad = file.split('_')[1].split('.')[0] # splits just the padded scan number string out
            pad_length = pad_length + len(pad)
            break
    return pad_length
    
def scanSuffix(scan_num_start, scan_num_end, padlength):
    scan_list = list(range(int(scan_num_start), int(scan_num_end)+1))
    scan_suffix = []
    for suffix in scan_list:
        scan_suffix.append(str(suffix).zfill(padlength))
    return scan_suffix
    
def runDifx(exp_tag, scan_start, scan_end, cwd):
    with open(exp_tag + '_reference_table.out', 'r') as ref_file:
        # read reference table for param information and v2d file name references
        reference_lines = ref_file.readlines()
        # create data table - ready to have data appended to it
        with open(exp_tag + '_data_table.out', 'w') as data_table: # write out the time data into file.
            print(reference_lines[0].rstrip() + ' corr_time', file=data_table)        
        # start running difx loop on each v2d file
        for i in range(1,len(reference_lines)):
            file = reference_lines[i].split()[0]
            try:
                print('vex2difx ' + file)
                os.system('vex2difx ' + file)
            except:
                print('vex2difx on ' + file + '.v2d failed, moving on...')
                continue
            pad_length = padCheck(file, cwd)
            total_time = 0
            if pad_length > 0:
                scan_suffix_list = scanSuffix(scan_start, scan_end, pad_length)
                for suffix in scan_suffix_list:
                    scan_file_name = file.split('.')[0] + '_' + suffix
                    try:
                        print('difxcalc ' + scan_file_name + '.calc')
                        #os.system('difxcalc ' + scan_file_name + '.calc')
                        print('correlating ' + scan_file_name + '.input')
                        time_delta = corrScan(scan_file_name, 172) # correlate the scan and determine time taken
                        total_time = total_time + time_delta
                    except:
                        continue
            with open(exp_tag + '_data_table.out', 'a') as data_table:
                print(reference_lines[i].rstrip() + ' ' + str(total_time), file=data_table) # append time to file + combination string
            
def main(v2d_file_name, name_input_file, start_scan_num, end_scan_num):
    exp_tag = v2d_file_name.split('.')[0]
    # determine possible combinations of parameters
    combinations, all_params, section_ref = paramInput(name_input_file)
    # generate v2d files for all parameter combos
    generateV2D(v2d_file_name, combinations, all_params, section_ref)
    # run difx on desired range of scans for all v2d files.
    runDifx(exp_tag, start_scan_num, end_scan_num, cwd)
            
if __name__ == '__main__':
    # difxParmaTest.py executed as a script
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])            

