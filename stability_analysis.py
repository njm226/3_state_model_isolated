#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 14:17:18 2022

@author: fabio
"""


import numpy as np 
#import pandas as pd
import matplotlib.pyplot as plt
#from smaller_model_function_AtoU import simple_small as ss
from all_local_model_40 import simple as ss 
import multiprocessing
#import time
import pickle
pool = multiprocessing.Pool(multiprocessing.cpu_count())


#noise = [0,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1, 0.5, 1, 5, 10]
noise = [0,0.0001,0.0002,0.0003,0.0004,0.0005,0.0006,0.0007,0.0008,0.0009,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10]
fraction_Citrine_ON_day_5_NS = []
fraction_mCherry_ON_day_5_NS = []

fraction_Citrine_ON_day_5_5kb = []
fraction_mCherry_ON_day_5_5kb = []

for n in noise: 
    parameters_NS=[[20,100,100,n]]
    parameters_5kb=[[40,100,100,n]]

    reps=1000
    
    repeat_NS=reps*parameters_NS
    repeat_5kb=reps*parameters_5kb

    duration=121*3#201
    
    if __name__ == '__main__':
        status_NS = pool.map(ss, repeat_NS) 
        status_5kb = pool.map(ss, repeat_5kb) 


    Citrine_list_NS = np.zeros([len(repeat_NS),duration]) 
    mCherry_list_NS = np.zeros([len(repeat_NS),duration])
    
    Citrine_list_5kb = np.zeros([len(repeat_5kb),duration]) 
    mCherry_list_5kb = np.zeros([len(repeat_5kb),duration])


    for elt in range(len(repeat_NS)):
        
        Citrine_NS = np.array(status_NS[elt][0])
        mCherry_NS = np.array(status_NS[elt][1])
        
        #switch the values of the list (1 stands now for timepoint when reporter is on)
        Citrine_NS=1-Citrine_NS
        mCherry_NS=1-mCherry_NS
        
        Citrine_list_NS[elt]=Citrine_NS
        mCherry_list_NS[elt]=mCherry_NS
        
        
        
        Citrine_5kb = np.array(status_5kb[elt][0])
        mCherry_5kb = np.array(status_5kb[elt][1])
        
        #switch the values of the list (1 stands now for timepoint when reporter is on)
        Citrine_5kb=1-Citrine_5kb
        mCherry_5kb=1-mCherry_5kb
        
        Citrine_list_5kb[elt]=Citrine_5kb
        mCherry_list_5kb[elt]=mCherry_5kb
    
    
    
    
    #output
    Citrine_list_NS = (sum(Citrine_list_NS))/reps
    #
    mCherry_list_NS = (sum(mCherry_list_NS))/reps
    
    #output
    Citrine_list_5kb = (sum(Citrine_list_5kb))/reps
    #
    mCherry_list_5kb = (sum(mCherry_list_5kb))/reps
    
    
    fraction_Citrine_ON_day_5_NS.append(Citrine_list_NS[120*3])
    fraction_mCherry_ON_day_5_NS.append(mCherry_list_NS[120*3])
    
    fraction_Citrine_ON_day_5_5kb.append(Citrine_list_5kb[120*3])
    fraction_mCherry_ON_day_5_5kb.append(mCherry_list_5kb[120*3])



with open('fraction_Citrine_ON_day_5_NS_TRE_0_silencing_long_with_rep.txt', 'wb') as F:
        pickle.dump(fraction_Citrine_ON_day_5_NS, F)
        
with open('fraction_mCherry_ON_day_5_NS_TRE_0_silencing_long_with_rep.txt', 'wb') as F:
        pickle.dump(fraction_mCherry_ON_day_5_NS, F)
        
with open('fraction_Citrine_ON_day_5_5kb_TRE_0_silencing_long_with_rep.txt', 'wb') as F:
        pickle.dump(fraction_Citrine_ON_day_5_5kb, F)
        
with open('fraction_mCherry_ON_day_5_5kb_TRE_0_silencing_long_with_rep.txt', 'wb') as F:
        pickle.dump(fraction_mCherry_ON_day_5_5kb, F)