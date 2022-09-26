#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 13:16:19 2019
histogramm and timecourse data
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


# noise = [0.0001,0.001,0.01,0.1,1]
# for n in noise:

parameters_silencing_NS = [[20,0,200,0.5,0]]  # N, starting state, duration, noise, TRE
parameters_silencing_5kb =[[40,0,200,0.5,0]]

parameters_silencing_TRE_NS = [[20,0,200,0.5,1000]]  # N, starting state, duration, noise, TRE
parameters_silencing_TRE_5kb =[[40,0,200,0.5,1000]]

parameters_reactivation_NS = [[20,2,400,0.5,0]]
parameters_reactivation_5kb =[[40,2,400,0.5,0]]

parameters_reactivation_SC = [[21,2,800,0.5,0]]



reps=1000

repeat_silencing_NS =reps*parameters_silencing_NS
repeat_silencing_5kb =reps*parameters_silencing_5kb

repeat_silencing_TRE_NS =reps*parameters_silencing_TRE_NS
repeat_silencing_TRE_5kb =reps*parameters_silencing_TRE_5kb

repeat_reactivation_NS =reps*parameters_reactivation_NS
repeat_reactivation_5kb =reps*parameters_reactivation_5kb


repeat_reactivation_SC =reps*parameters_reactivation_SC


if __name__ == '__main__':
    status_silencing_NS = pool.map(ss, repeat_silencing_NS) 
    status_silencing_5kb = pool.map(ss, repeat_silencing_5kb) 
    
    status_silencing_TRE_NS = pool.map(ss, repeat_silencing_TRE_NS) 
    status_silencing_TRE_5kb = pool.map(ss, repeat_silencing_TRE_5kb) 
    
    status_reactivation_NS = pool.map(ss, repeat_reactivation_NS) 
    status_reactivation_5kb = pool.map(ss, repeat_reactivation_5kb) 
    
    status_reactivation_SC = pool.map(ss, repeat_reactivation_SC) 


Citrine_list_silencing_NS = np.zeros([len(repeat_silencing_NS),200*3]) 
mCherry_list_silencing_NS = np.zeros([len(repeat_silencing_NS),200*3])

Citrine_list_silencing_5kb = np.zeros([len(repeat_silencing_5kb),200*3]) 
mCherry_list_silencing_5kb = np.zeros([len(repeat_silencing_5kb),200*3])


Citrine_list_silencing_TRE_NS = np.zeros([len(repeat_silencing_TRE_NS),200*3]) 
mCherry_list_silencing_TRE_NS = np.zeros([len(repeat_silencing_TRE_NS),200*3])

Citrine_list_silencing_TRE_5kb = np.zeros([len(repeat_silencing_TRE_5kb),200*3]) 
mCherry_list_silencing_TRE_5kb = np.zeros([len(repeat_silencing_TRE_5kb),200*3])


Citrine_list_reactivation_NS = np.zeros([len(repeat_reactivation_NS),400*3]) 
mCherry_list_reactivation_NS = np.zeros([len(repeat_reactivation_NS),400*3])

Citrine_list_reactivation_5kb = np.zeros([len(repeat_reactivation_5kb),400*3]) 
mCherry_list_reactivation_5kb = np.zeros([len(repeat_reactivation_5kb),400*3])


Citrine_list_reactivation_SC = np.zeros([len(repeat_reactivation_SC),800*3])
mCherry_list_reactivation_SC = np.zeros([len(repeat_reactivation_SC),800*3])


delay_list_NS = np.zeros([len(repeat_silencing_NS)])
delay_list_5kb = np.zeros([len(repeat_silencing_5kb)])



for elt in range(len(repeat_silencing_NS)):
    
    Citrine_silencing_TRE_NS = np.array(status_silencing_TRE_NS[elt][0])
    mCherry_silencing_TRE_NS = np.array(status_silencing_TRE_NS[elt][1])
    

    delay_list_NS[elt] = list(mCherry_silencing_TRE_NS).index(1) - list(Citrine_silencing_TRE_NS).index(1)

    
    #switch the values of the list (1 stands now for timepoint when reporter is on)
    Citrine_silencing_TRE_NS=1-Citrine_silencing_TRE_NS
    mCherry_silencing_TRE_NS=1-mCherry_silencing_TRE_NS
    
    Citrine_list_silencing_TRE_NS[elt]=Citrine_silencing_TRE_NS
    mCherry_list_silencing_TRE_NS[elt]=mCherry_silencing_TRE_NS
    
    
    
    
    Citrine_silencing_TRE_5kb = np.array(status_silencing_TRE_5kb[elt][0])
    mCherry_silencing_TRE_5kb = np.array(status_silencing_TRE_5kb[elt][1])
    

    delay_list_5kb[elt] = list(mCherry_silencing_TRE_5kb).index(1) - list(Citrine_silencing_TRE_5kb).index(1)

    
    #switch the values of the list (1 stands now for timepoint when reporter is on)
    Citrine_silencing_TRE_5kb=1-Citrine_silencing_TRE_5kb
    mCherry_silencing_TRE_5kb=1-mCherry_silencing_TRE_5kb
    
    Citrine_list_silencing_TRE_5kb[elt]=Citrine_silencing_TRE_5kb
    mCherry_list_silencing_TRE_5kb[elt]=mCherry_silencing_TRE_5kb
    
    
    
    Citrine_silencing_NS = np.array(status_silencing_NS[elt][0])
    mCherry_silencing_NS = np.array(status_silencing_NS[elt][1])
    
    
    #switch the values of the list (1 stands now for timepoint when reporter is on)
    Citrine_silencing_NS=1-Citrine_silencing_NS
    mCherry_silencing_NS=1-mCherry_silencing_NS
    
    Citrine_list_silencing_NS[elt]=Citrine_silencing_NS
    mCherry_list_silencing_NS[elt]=mCherry_silencing_NS
    
    
    
    
    Citrine_silencing_5kb = np.array(status_silencing_5kb[elt][0])
    mCherry_silencing_5kb = np.array(status_silencing_5kb[elt][1])
    
    
    #switch the values of the list (1 stands now for timepoint when reporter is on)
    Citrine_silencing_5kb=1-Citrine_silencing_5kb
    mCherry_silencing_5kb=1-mCherry_silencing_5kb
    
    Citrine_list_silencing_5kb[elt]=Citrine_silencing_5kb
    mCherry_list_silencing_5kb[elt]=mCherry_silencing_5kb
    
    
    
    
    
    Citrine_reactivation_NS = np.array(status_reactivation_NS[elt][0])
    mCherry_reactivation_NS = np.array(status_reactivation_NS[elt][1])
    
    #switch the values of the list (1 stands now for timepoint when reporter is on)
    Citrine_reactivation_NS=1-Citrine_reactivation_NS
    mCherry_reactivation_NS=1-mCherry_reactivation_NS
    
    Citrine_list_reactivation_NS[elt]=Citrine_reactivation_NS
    mCherry_list_reactivation_NS[elt]=mCherry_reactivation_NS
    
    

    
    Citrine_reactivation_5kb = np.array(status_reactivation_5kb[elt][0])
    mCherry_reactivation_5kb = np.array(status_reactivation_5kb[elt][1])
    
    #switch the values of the list (1 stands now for timepoint when reporter is on)
    Citrine_reactivation_5kb=1-Citrine_reactivation_5kb
    mCherry_reactivation_5kb=1-mCherry_reactivation_5kb
    
    Citrine_list_reactivation_5kb[elt]=Citrine_reactivation_5kb
    mCherry_list_reactivation_5kb[elt]=mCherry_reactivation_5kb
    
    
    
    
    
        
    Citrine_reactivation_SC = np.array(status_reactivation_SC[elt][0])
    mCherry_reactivation_SC = np.array(status_reactivation_SC[elt][1])
    
    #switch the values of the list (1 stands now for timepoint when reporter is on)
    Citrine_reactivation_SC=1-Citrine_reactivation_SC
    mCherry_reactivation_SC=1-mCherry_reactivation_SC
    
    Citrine_list_reactivation_SC[elt]=Citrine_reactivation_SC
    mCherry_list_reactivation_SC[elt]=mCherry_reactivation_SC

   



#output
Citrine_list_silencing_NS = (sum(Citrine_list_silencing_NS))/reps
#
mCherry_list_silencing_NS = (sum(mCherry_list_silencing_NS))/reps


Citrine_list_silencing_5kb = (sum(Citrine_list_silencing_5kb))/reps
#
mCherry_list_silencing_5kb = (sum(mCherry_list_silencing_5kb))/reps


Citrine_list_silencing_TRE_NS = (sum(Citrine_list_silencing_TRE_NS))/reps
#
mCherry_list_silencing_TRE_NS = (sum(mCherry_list_silencing_TRE_NS))/reps


Citrine_list_silencing_TRE_5kb = (sum(Citrine_list_silencing_TRE_5kb))/reps
#
mCherry_list_silencing_TRE_5kb = (sum(mCherry_list_silencing_TRE_5kb))/reps


Citrine_list_reactivation_NS = (sum(Citrine_list_reactivation_NS))/reps
#
mCherry_list_reactivation_NS = (sum(mCherry_list_reactivation_NS))/reps


Citrine_list_reactivation_5kb = (sum(Citrine_list_reactivation_5kb))/reps
#
mCherry_list_reactivation_5kb = (sum(mCherry_list_reactivation_5kb))/reps


Citrine_list_reactivation_SC = (sum(Citrine_list_reactivation_SC))/reps
#
mCherry_list_reactivation_SC = (sum(mCherry_list_reactivation_SC))/reps



with open('Citrine_list_silencing_TRE_NS.txt', 'wb') as F:
    pickle.dump(Citrine_list_silencing_TRE_NS, F)
     
with open('mCherry_list_silencing_TRE_NS.txt', 'wb') as F:
    pickle.dump(mCherry_list_silencing_TRE_NS, F)


with open('Citrine_list_silencing_TRE_5kb.txt', 'wb') as F:
    pickle.dump(Citrine_list_silencing_TRE_5kb, F)
     
with open('mCherry_list_silencing_TRE_5kb.txt', 'wb') as F:
    pickle.dump(mCherry_list_silencing_TRE_5kb, F)


with open('Citrine_list_silencing_NS.txt', 'wb') as F:
    pickle.dump(Citrine_list_silencing_NS, F)
     
with open('mCherry_list_silencing_NS.txt', 'wb') as F:
    pickle.dump(mCherry_list_silencing_NS, F)


with open('Citrine_list_silencing_5kb.txt', 'wb') as F:
    pickle.dump(Citrine_list_silencing_5kb, F)
     
with open('mCherry_list_silencing_5kb.txt', 'wb') as F:
    pickle.dump(mCherry_list_silencing_5kb, F)
    
    
with open('Citrine_list_reactivation_NS.txt', 'wb') as F:
    pickle.dump(Citrine_list_reactivation_NS, F)
     
with open('mCherry_list_reactivation_NS.txt', 'wb') as F:
    pickle.dump(mCherry_list_reactivation_NS, F)

    
with open('Citrine_list_reactivation_5kb.txt', 'wb') as F:
    pickle.dump(Citrine_list_reactivation_5kb, F)
     
with open('mCherry_list_reactivation_5kb.txt', 'wb') as F:
    pickle.dump(mCherry_list_reactivation_5kb, F)
    
    
with open('Citrine_list_reactivation_SC.txt', 'wb') as F:
    pickle.dump(Citrine_list_reactivation_SC, F)
     
with open('mCherry_list_reactivation_SC.txt', 'wb') as F:
    pickle.dump(mCherry_list_reactivation_SC, F)
    
    
with open('delay_list_NS.txt', 'wb') as F:
    pickle.dump(delay_list_NS, F)
    
with open('delay_list_5kb.txt', 'wb') as F:
    pickle.dump(delay_list_5kb, F)





