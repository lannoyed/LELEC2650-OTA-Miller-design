# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 14:32:52 2022

@author: dlannoye
"""

import numpy as np


import os
os.system('cmd /k "Your Command Prompt Command"')

def run_sim():
    os.system('.csh')
    os.system('eldo to_test.cir')


def modification (line, File_new, Cf, CL, W, L, VIN, IBIAS):
    if line == '.param Cl = 10p\n' : 
        File_new.write(f'.param Cl = {CL}p\n')
    elif line == '.param Cf = 25p\n' : 
        File_new.write(f'.param Cf = {Cf}p\n')
    
    elif line.split('=')[0] == '.param W1 ' : 
        File_new.write(f'.param W1 = {W[0]}u\n')
    elif line.split('=')[0] == '.param L1 ' : 
        File_new.write(f'.param L1 = {L[0]}u\n')
        
    elif line.split('=')[0] == '.param W2 ' : 
        File_new.write(f'.param W2 = {W[1]}u\n')
    elif line.split('=')[0] == '.param L2 ' : 
        File_new.write(f'.param L2 = {L[1]}u\n')

    elif line.split('=')[0] == '.param W3 ' : 
        File_new.write(f'.param W3 = {W[2]}u\n')
    elif line.split('=')[0] == '.param L3 ' : 
        File_new.write(f'.param L3 = {L[2]}u\n')

    elif line.split('=')[0] == '.param W4 ' : 
        File_new.write(f'.param W4 = {W[3]}u\n')
    elif line.split('=')[0] == '.param L4 ' : 
        File_new.write(f'.param L4 = {L[3]}u\n')

    elif line.split('=')[0] == '.param W5 ' : 
        File_new.write(f'.param W5 = {W[4]}u\n')
    elif line.split('=')[0] == '.param L5 ' : 
        File_new.write(f'.param L5 = {L[4]}u\n')

    elif line.split('=')[0] == '.param W6 ' : 
        File_new.write(f'.param W6 = {W[5]}u\n')
    elif line.split('=')[0] == '.param L6 ' : 
        File_new.write(f'.param L6 = {L[5]}u\n')

    elif line.split('=')[0] == '.param W7 ' : 
        File_new.write(f'.param W7 = {W[6]}u\n')
    elif line.split('=')[0] == '.param L7 ' : 
        File_new.write(f'.param L7 = {L[6]}u\n')

    elif line.split('=')[0] == '.param W8 ' : 
        File_new.write(f'.param W8 = {W[7]}u\n')
    elif line.split('=')[0] == '.param L8 ' : 
        File_new.write(f'.param L8 = {L[7]}u\n')     
        
        
    elif line.split('=')[0] == '.param IBIAS ':
        File_new.write(f'.param VIN = {IBIAS}u\n')
    elif line.split('=')[0] == '.param VIN ':
        File_new.write(f'.param VIN = {VIN}\n')        
    else :
        File_new.write(line)




if __name__ == '__main__':
    p = 1e-12
    u = 1e-6
    
    W = np.array([1,2,3,4,5,6,7,8]) #µm
    L = np.array([8,7,6,5,4,3,2,1]) #µm
    CL = 20 #pF
    Cf = 25 #pF
    
    VIN = 2.5 
    IBIAS = 69 #µA
    
    File = open('Miller_OTA_Basic.cir', 'r')
    File_new =  open('to_test.cir', 'w')
    
    for line in File : 
        modification (line, File_new, Cf, CL, W, L, VIN, IBIAS)
        
    
    File_new.close()
    
    run_sim()    
    