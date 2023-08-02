# -*- coding: utf-8 -*-
"""
Created on Wed May 31 19:34:59 2023

@author: Admin
"""
import os
import shutil
import test_generator as gen
import set_files_in_order as order

def main():
    '''set up the file to store all test sets (3 of)'''
    set_pointer = 'D:/test_gensis_A/set_3'
    os.mkdir(set_pointer)
    '''for each set, with model_index (int) stat from 0, form file (pointer)'''
    for model_index in range(3):
        model_index = model_index + 1
        gen.test_generator(model_index = model_index)
        '''put simulated results into pointer files
           1. get pointer and pop files
           2. transplant them
        '''
        pointer = order.create_pointer(set_pointer, model_index)
        L_popfiles = order.get_test_set_pointer(pointer)
        
        for popfile in L_popfiles:
            order.move_cor_files_to(popfile)
        
    return

main()