# -*- coding: utf-8 -*-
"""
Created on Mon May 15 11:07:44 2023

@author: Admin
"""

import os
import shutil

#src_str = 'C:/Users/Admin/Desktop/pop6.0'
def create_pointer(set_pointer, model_index):
    pointer = set_pointer+'/testset'+str(model_index)
    os.mkdir(pointer)
    return pointer

def get_test_set_pointer(pointer,l_pop = ['5.7','6.0','6.7','7.0','7.7','8.0','8.7','9.0']):
    L_popfile = [pointer+'/'+'pop'+pop for pop in l_pop]
    for popfile in L_popfile:
        os.mkdir(popfile)
    return L_popfile
    
def move_cor_files_to(popfile,gensis = 'D:/test_gensis_A/test_file_generator'):
    L = popfile.split('/')
    #print(L)
    L1 = L[-1].split('p')
    #print(L1)
    pophead = L1[-1]
    #print(pophead)
    orgtail = '/orgs'
    popcountail = '/popcounts'
    rawtail = '/raw'
    
    L_cor = ['cor'+str(round(0.1*(i+1),1)) for i in range(9)]
    
    for cor in L_cor:
        
        popcorfile = popfile+'/'+cor
        os.mkdir(popcorfile)
        orgfile = popcorfile + orgtail
        os.mkdir(orgfile)
        popcountfile = popcorfile + popcountail
        os.mkdir(popcountfile)
        rawfile = popcorfile + rawtail
        os.mkdir(rawfile)
        
        L_file = os.listdir(gensis)
        for file in L_file:
            if pophead in file and cor in file:
                #print(pophead, cor, file)s
                if 'org' in file:
                    shutil.move(file,orgfile)
                elif 'raw' in file:
                    shutil.move(file,rawfile)
                else:
                    shutil.move(file,popcountfile)
            else:
                pass
    return
