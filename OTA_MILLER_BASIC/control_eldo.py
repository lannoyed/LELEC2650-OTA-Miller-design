# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 14:32:52 2022

@author: dlannoye


to run :


csh
python3 -m venv env
source env/bin/activate.csh

pip install --upgrade pip
pip install numpy
pip install chardet
pip install matplotlib
pip install scipy


"""

import numpy as np


import os


import Miller_sizing



def run_sim():
    os.system('eldo to_test.cir')


def modification(line, File_new, Cf, CL, W, L, VIN, IBIAS):
    
    
    if line == '.param Cl = 10p\n':
        File_new.write('.param Cl = {}p\n'.format(CL))
    elif line == '.param Cf = 25p\n':
        File_new.write('.param Cf = {}p\n'.format(Cf))
    elif line.split('=')[0] == '.param W1 ':
        File_new.write('.param W1 = {}u\n'.format(W[0]))
    elif line.split('=')[0] == '.param L1 ':
        File_new.write('.param L1 = {}u\n'.format(L[0]))
    elif line.split('=')[0] == '.param W2 ':
        File_new.write('.param W2 = {}u\n'.format(W[1]))
    elif line.split('=')[0] == '.param L2 ':
        File_new.write('.param L2 = {}u\n'.format(L[1]))
    elif line.split('=')[0] == '.param W3 ':
        File_new.write('.param W3 = {}u\n'.format(W[2]))
    elif line.split('=')[0] == '.param L3 ':
        File_new.write('.param L3 = {}u\n'.format(L[2]))
    elif line.split('=')[0] == '.param W4 ':
        File_new.write('.param W4 = {}u\n'.format(W[3]))
    elif line.split('=')[0] == '.param L4 ':
        File_new.write('.param L4 = {}u\n'.format(L[3]))
    elif line.split('=')[0] == '.param W5 ':
        File_new.write('.param W5 = {}u\n'.format(W[4]))
    elif line.split('=')[0] == '.param L5 ':
        File_new.write('.param L5 = {}u\n'.format(L[4]))
    elif line.split('=')[0] == '.param W6 ':
        File_new.write('.param W6 = {}u\n'.format(W[5]))
    elif line.split('=')[0] == '.param L6 ':
        File_new.write('.param L6 = {}u\n'.format(L[5]))
    elif line.split('=')[0] == '.param W7 ':
        File_new.write('.param W7 = {}u\n'.format(W[6]))
    elif line.split('=')[0] == '.param L7 ':
        File_new.write('.param L7 = {}u\n'.format(L[6]))
    elif line.split('=')[0] == '.param W8 ':
        File_new.write('.param W8 = {}u\n'.format(W[7]))
    elif line.split('=')[0] == '.param L8 ':
        File_new.write('.param L8 = {}u\n'.format(L[7]))
    elif line.split('=')[0] == '.param IBIAS ':
        File_new.write('.param IBIAS = {}u\n'.format(IBIAS))
    elif line.split('=')[0] == '.param VIN ':
        File_new.write('.param VIN = {}\n'.format(VIN))
    else:
        File_new.write(line)

def safe_solution(W_list, VIN, IBIAS):
    if VIN > 1.1 or VIN < 0.0:
        return False
    
    if IBIAS > 1000 or IBIAS < 0.0:
        return False
    
    for i in range(len(W_list)):
        if W_list[i] > 1000 or W_list[i] < 0.0:
            return False
    return True

def all_transistor_satured (list_vds, list_vdsat):
    for i in range(len(list_vds)):
        if np.abs(list_vds[i]) < np.abs(list_vdsat[i]):
            return 0
    return 1
        
import chardet
import codecs

import subprocess


if __name__ == '__main__':
    p = 1e-12
    u = 1e-6
    
    L = 1e-6 #um
    
    fT = 13.7e06
    CL = 10e-12 
    Cf = 3e-9
    
    gmid1_list = np.arange(14, 16,1, dtype=float)
    #do the same for the other gmid
    gmid3_list = np.arange(12, 14, dtype=float)
    
    gmid5_list = np.arange(9, 11, dtype=float)
    gmid6_list = np.arange(9, 11, dtype=float)
    gmid7_list = np.arange(12, 14, dtype=float)
    gmid8_list = np.arange(9, 11, dtype=float)
    
    output_file = open('output_file_optimisation.txt', 'w', encoding='utf-8')
    output_file.write('gmid1 gmid3 gmid5 gmid6 gmid7 gmid8 all_satured AV0 FT PHASE_MARGINE\n')
      
    for gmid1 in gmid1_list: 
        for gmid3 in gmid3_list: 
            for gmid5 in gmid5_list: 
                for gmid6 in gmid6_list: 
                    for gmid7 in gmid7_list: 
                        for gmid8 in gmid8_list: 
    
                            Cf_result, CL_result, W_list, L_list, IBIAS, VIN, vdsat_list = Miller_sizing.Miller_sizing(Cf, fT, CL, L, gmid1, gmid3, gmid5, gmid6, gmid7, gmid8)  
                                                      
                            if safe_solution(W_list, VIN, IBIAS):
                                # Detect the encoding of the file
                                
                                with open('Miller_OTA_Basic.cir', 'rb') as f:
                                    encoding = chardet.detect(f.read())['encoding']
                                
                                with codecs.open('Miller_OTA_Basic.cir', 'r', encoding=encoding) as File:
                                    with codecs.open('to_test.cir', 'w', encoding=encoding) as File_new:
                                        for line in File:
                                            modification(line, File_new, Cf_result, CL_result, W_list, L_list, VIN, IBIAS)
            
                                      
                                # Start the program in a separate process
                                process = subprocess.Popen(['eldo', 'to_test.cir'], stdout=subprocess.PIPE)
                                
                                # Create a list to store the last 3 lines of output
                                vds = []
                                Av0 = 0
                                FT = 0
                                phase_margin = 0
                                
                                # Read 3 lines from the process's standard output stream
                                for i in range(1000):
                                    line = str( process.stdout.readline()) 
                                    if not line:
                                        break
                                    
                                    line = line.split(' ')
                                    
                                    if len(line) > 4:
                                        value = line[1]    
                                        vds_str = line[4]
                                    
                                        if value == 'VSD1':
                                            vds.append(float(vds_str))
                                        elif value == 'VSD2':
                                            vds.append(float(vds_str))
                                        elif value == 'VDS3':
                                            vds.append(float(vds_str))
                                        elif value == 'VDS4':
                                            vds.append(float(vds_str))
                                        elif value == 'VSD5':
                                            vds.append(float(vds_str))
                                        elif value == 'VSD6':
                                            vds.append(float(vds_str))
                                        elif value == 'VDS7':
                                            vds.append(float(vds_str))
                                        elif value == 'VSD8':
                                            vds.append(float(vds_str))
                                        elif value == 'AV0':
                                            Av0 = float(vds_str)
                                        elif value == 'FT':
                                            FT = float(vds_str)
                                        elif value == 'PHASE_MARGIN':
                                            phase_margin = float(vds_str)
                                            break
                                        

                                # Print the last 3 lines of output
                                print(vds[:])
                                
                                output_file.write('{} {} {} {} {} {} {} {} {} {}\n'.format(gmid1, gmid3, gmid5, gmid6, gmid7, gmid8, all_transistor_satured(vds, vdsat_list), Av0, FT, phase_margin))
                                                                                  
                                
                                
                            else :
                                print('\n\n\n\nSolution not safe \n\n\n\n')

                            
                            
