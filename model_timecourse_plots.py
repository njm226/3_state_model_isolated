#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 13:16:19 2019
histogramm and timecourse data
@author: fabio
"""


import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from pylab import *
from scipy import stats
import uncertainties as unc
import uncertainties.unumpy as unp

import pickle


duration1 = 200*3#200*3#40
time1 = np.array(range(duration1))


with open('Citrine_list_silencing_TRE_NS.txt', 'rb') as F:
    Citrine_list_silencing_TRE_NS = pickle.load(F)
     
with open('mCherry_list_silencing_TRE_NS.txt', 'rb') as F:
    mCherry_list_silencing_TRE_NS = pickle.load(F)


with open('Citrine_list_silencing_TRE_5kb.txt', 'rb') as F:
    Citrine_list_silencing_TRE_5kb = pickle.load(F)
     
with open('mCherry_list_silencing_TRE_5kb.txt', 'rb') as F:
    mCherry_list_silencing_TRE_5kb = pickle.load(F)
    

#plot the fitted function
plt.figure(figsize=[15,10])
plt.plot(time1, Citrine_list_silencing_TRE_NS, 'b', label='Citrine NS', lw=3.)#: a=%5.1f, b=%5.1f' % tuple(popt1), lw=3.)
plt.plot(time1, mCherry_list_silencing_TRE_NS, 'g', label='mCherry NS', lw=3.)#: a=%5.1f, b=%5.1f' % tuple(popt3), lw=3.)
plt.plot(time1, Citrine_list_silencing_TRE_5kb, 'r', label='Citrine 5kb', lw=3.)#: a=%5.1f, b=%5.1f' % tuple(popt1), lw=3.)
plt.plot(time1, mCherry_list_silencing_TRE_5kb, 'k', label='mCherry 5kb', lw=3.)#: a=%5.1f, b=%5.1f' % tuple(popt3), lw=3.)

#plt.ylim([-0.5,1])
plt.yticks(fontsize=40)
plt.xticks(fontsize=40)
plt.xlim([0,120])
#plt.yscale('log')
#plt.xticks([5, 10, 15, 20, 25, 30, 35, 40])
#plt.xticks(ticks=[20*3, 40*3, 60*3, 80*3, 100*3, 120*3, 140*3, 160*3, 180*3, 200*3], labels=[20, 40, 60, 80, 100, 120, 140, 160, 180, 200])
plt.xticks(ticks=[24*3, 48*3, 72*3, 96*3, 120*3], labels=[1,2,3,4,5])
#plt.xticks(ticks=[0, 120*3, 240*3, 360*3], labels=[0,5,10,15])
#plt.xticks(ticks=[120*3, 240*3, 360*3, 480*3, 600*3, 720*3], labels=[5, 10, 15, 20, 25,30])
plt.yticks([1,0.5,0],[1, 0.5, 0])
#plt.yticks([0,0.1,1],[0, 0.1, 1])
plt.ylim([0,1])
#plt.ylabel('fraction of "ON" cells', fontsize=30)
#plt.xlabel('time (generations)', fontsize=30)
plt.legend(fontsize=30)#, loc='upper right')
plt.tick_params(width=4,length=4)

plt.savefig('silencing_TRE.pdf', format='pdf')


with open('Citrine_list_silencing_NS.txt', 'rb') as F:
    Citrine_list_silencing_NS = pickle.load(F)
     
with open('mCherry_list_silencing_NS.txt', 'rb') as F:
    mCherry_list_silencing_NS = pickle.load(F)


with open('Citrine_list_silencing_5kb.txt', 'rb') as F:
    Citrine_list_silencing_5kb = pickle.load(F)
     
with open('mCherry_list_silencing_5kb.txt', 'rb') as F:
    mCherry_list_silencing_5kb = pickle.load(F)
    

#plot the fitted function
plt.figure(figsize=[15,10])
plt.plot(time1, Citrine_list_silencing_NS, 'b', label='Citrine NS', lw=3.)#: a=%5.1f, b=%5.1f' % tuple(popt1), lw=3.)
plt.plot(time1, mCherry_list_silencing_NS, 'g', label='mCherry NS', lw=3.)#: a=%5.1f, b=%5.1f' % tuple(popt3), lw=3.)
plt.plot(time1, Citrine_list_silencing_5kb, 'r', label='Citrine 5kb', lw=3.)#: a=%5.1f, b=%5.1f' % tuple(popt1), lw=3.)
plt.plot(time1, mCherry_list_silencing_5kb, 'k', label='mCherry 5kb', lw=3.)#: a=%5.1f, b=%5.1f' % tuple(popt3), lw=3.)

#plt.ylim([-0.5,1])
plt.yticks(fontsize=40)
plt.xticks(fontsize=40)
plt.xlim([0,120])
#plt.yscale('log')
#plt.xticks([5, 10, 15, 20, 25, 30, 35, 40])
#plt.xticks(ticks=[20*3, 40*3, 60*3, 80*3, 100*3, 120*3, 140*3, 160*3, 180*3, 200*3], labels=[20, 40, 60, 80, 100, 120, 140, 160, 180, 200])
plt.xticks(ticks=[24*3, 48*3, 72*3, 96*3, 120*3], labels=[1,2,3,4,5])
#plt.xticks(ticks=[0, 120*3, 240*3, 360*3], labels=[0,5,10,15])
#plt.xticks(ticks=[120*3, 240*3, 360*3, 480*3, 600*3, 720*3], labels=[5, 10, 15, 20, 25,30])
plt.yticks([1,0.5,0],[1, 0.5, 0])
#plt.yticks([0,0.1,1],[0, 0.1, 1])
plt.ylim([0,1])
#plt.ylabel('fraction of "ON" cells', fontsize=30)
#plt.xlabel('time (generations)', fontsize=30)
plt.legend(fontsize=30)#, loc='upper right')
plt.tick_params(width=4,length=4)

plt.savefig('silencing.pdf', format='pdf')


## save state_list
with open('delay_list_NS.txt', 'rb') as F:
    delay_NS = pickle.load(F)
    
## save state_list
with open('delay_list_5kb.txt', 'rb') as F:
    delay_5kb= pickle.load(F)


fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(4,8))

# Plot violin plot on axes 1
ax1.violinplot([delay_NS, delay_5kb],showmedians=True, showextrema=False, bw_method='scott')
ax1.set_ylim([-15,150])
ax1.set_xticks([1,2])
ax1.set_xticklabels(['NS','5kb'], fontsize=30)
ax1.set_yticks([-15,0,30,60,90,120,150])
ax1.set_yticklabels([-5,0,10,20,30,40,50], fontsize=20)


plt.show()
plt.savefig('delay.pdf')



duration2 = 400*3#200*3#40
time2 = np.array(range(duration2))


with open('Citrine_list_reactivation_NS.txt', 'rb') as F:
    Citrine_list_reactivation_NS = pickle.load(F)
     
with open('mCherry_list_reactivation_NS.txt', 'rb') as F:
    mCherry_list_reactivation_NS = pickle.load(F)


with open('Citrine_list_reactivation_5kb.txt', 'rb') as F:
    Citrine_list_reactivation_5kb = pickle.load(F)
     
with open('mCherry_list_reactivation_5kb.txt', 'rb') as F:
    mCherry_list_reactivation_5kb = pickle.load(F)




#plot the fitted function
plt.figure(figsize=[15,10])
plt.plot(time2, Citrine_list_reactivation_NS, 'b', label='Citrine NS', lw=3.)#: a=%5.1f, b=%5.1f' % tuple(popt1), lw=3.)
plt.plot(time2, mCherry_list_reactivation_NS, 'g', label='mCherry NS', lw=3.)#: a=%5.1f, b=%5.1f' % tuple(popt3), lw=3.)
plt.plot(time2, Citrine_list_reactivation_5kb, 'r', label='Citrine 5kb', lw=3.)#: a=%5.1f, b=%5.1f' % tuple(popt1), lw=3.)
plt.plot(time2, mCherry_list_reactivation_5kb, 'k', label='mCherry 5kb', lw=3.)#: a=%5.1f, b=%5.1f' % tuple(popt3), lw=3.)

#plt.ylim([-0.5,1])
plt.yticks(fontsize=40)
plt.xticks(fontsize=40)
plt.xlim([0,360])
#plt.yscale('log')
#plt.xticks([5, 10, 15, 20, 25, 30, 35, 40])
#plt.xticks(ticks=[20*3, 40*3, 60*3, 80*3, 100*3, 120*3, 140*3, 160*3, 180*3, 200*3], labels=[20, 40, 60, 80, 100, 120, 140, 160, 180, 200])
#plt.xticks(ticks=[24*3, 48*3, 72*3, 96*3, 120*3], labels=[1,2,3,4,5])
plt.xticks(ticks=[0, 120*3, 240*3, 360*3], labels=[0,5,10,15])
#plt.xticks(ticks=[120*3, 240*3, 360*3, 480*3, 600*3, 720*3], labels=[5, 10, 15, 20, 25,30])
plt.yticks([1,0.5,0],[1, 0.5, 0])
#plt.yticks([0,0.1,1],[0, 0.1, 1])
plt.ylim([0,1])
#plt.ylabel('fraction of "ON" cells', fontsize=30)
#plt.xlabel('time (generations)', fontsize=30)
plt.legend(fontsize=30)#, loc='upper right')
plt.tick_params(width=4,length=4)

plt.savefig('reactivation.pdf', format='pdf')


duration3 = 800*3#200*3#40
time3 = np.array(range(duration3))

with open('Citrine_list_reactivation_SC.txt', 'rb') as F:
    Citrine_list_reactivation_SC = pickle.load(F)
     
with open('mCherry_list_reactivation_SC.txt', 'rb') as F:
    mCherry_list_reactivation_SC = pickle.load(F)




#plot the fitted function
plt.figure(figsize=[15,10])
plt.plot(time3, Citrine_list_reactivation_SC, 'b', label='Citrine SC', lw=3.)#: a=%5.1f, b=%5.1f' % tuple(popt1), lw=3.)
plt.plot(time3, mCherry_list_reactivation_SC, 'g', label='mCherry SC', lw=3.)#: a=%5.1f, b=%5.1f' % tuple(popt3), lw=3.)

#plt.ylim([-0.5,1])
plt.yticks(fontsize=40)
plt.xticks(fontsize=40)
plt.xlim([0,720])
#plt.yscale('log')
#plt.xticks([5, 10, 15, 20, 25, 30, 35, 40])
#plt.xticks(ticks=[20*3, 40*3, 60*3, 80*3, 100*3, 120*3, 140*3, 160*3, 180*3, 200*3], labels=[20, 40, 60, 80, 100, 120, 140, 160, 180, 200])
#plt.xticks(ticks=[24*3, 48*3, 72*3, 96*3, 120*3], labels=[1,2,3,4,5])
#plt.xticks(ticks=[0, 120*3, 240*3, 360*3], labels=[0,5,10,15])
plt.xticks(ticks=[120*3, 240*3, 360*3, 480*3, 600*3, 720*3], labels=[5, 10, 15, 20, 25,30])
plt.yticks([1,0.5,0],[1, 0.5, 0])
#plt.yticks([0,0.1,1],[0, 0.1, 1])
plt.ylim([0,1])
#plt.ylabel('fraction of "ON" cells', fontsize=30)
#plt.xlabel('time (generations)', fontsize=30)
plt.legend(fontsize=30)#, loc='upper right')
plt.tick_params(width=4,length=4)

plt.savefig('reactivation_SC.pdf', format='pdf')



