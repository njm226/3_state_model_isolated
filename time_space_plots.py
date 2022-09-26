#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 11:52:01 2022

@author: fabio
"""


import numpy as np
import pyximport
pyximport.install(reload_support=True, setup_args={"include_dirs":np.get_include()})
from all_local_40 import t_loop
import all_local_40
import importlib
importlib.reload(all_local_40)
import pandas as pd
import matplotlib.pyplot as plt
import pickle

    
with open('System_state.txt', 'rb') as F:
    System_state = pickle.load(F)
#print(nucs[:100])
positions = list(range(40))
x = positions*(80)
y= [i for i in range(80) for _ in range(40)]


# show every 200th element

Nucs = System_state[0:2400:30]
#Nucs = System_state[720:960:3]

# convert 2D nucleosome list to a 1D list
States = Nucs.flatten()
print(len(States))
print(len(Nucs))

colors = []
for elt in States:
    if elt == 0:
        colors.append('blue')
    elif elt == 1:
        colors.append('lightgrey')
    elif elt==2:
        colors.append('red')


#plotting
fig, ax = plt.subplots(figsize=(6,21))
#plt.title("beta= 0.01, sigma=1", fontsize = 25)
plt.scatter(x, y, c=colors, marker="o", lw=0.6)
plt.yticks([10,20,30,40,50,60,70,80], [100,200,300,400,500,600,700,800],fontsize=20)
plt.ylabel("Generations", fontsize =25)
plt.xticks(fontsize=20)
#plt.savefig("timecourse_203_new.pdf")






#print(nucs[:100])
positions = list(range(40))
x = positions*(80)
y= [i for i in range(80) for _ in range(40)]


# show every 200th element

Nucs = System_state[1440:1680:3]
#Nucs = System_state[720:960:3]

# convert 2D nucleosome list to a 1D list
States = Nucs.flatten()
print(len(States))
print(len(Nucs))

colors = []
for elt in States:
    if elt == 0:
        colors.append('blue')
    elif elt == 1:
        colors.append('lightgrey')
    elif elt==2:
        colors.append('red')


#plotting
fig, ax = plt.subplots(figsize=(6,21))
#plt.title("beta= 0.01, sigma=1", fontsize = 25)
plt.scatter(x, y, c=colors, marker="o", lw=0.6)
plt.yticks([10,20,30,40,50,60,70,80], [510,520,530,540,550,560,570,580],fontsize=20)
plt.ylabel("Generations", fontsize =25)
plt.xticks(fontsize=20)
#plt.savefig("timecourse_203_new.pdf")









#print(nucs[:100])
positions = list(range(40))
x = positions*(80)
y= [i for i in range(80) for _ in range(40)]


# show every 200th element

Nucs = System_state[960:1200:3]
#Nucs = System_state[720:960:3]

# convert 2D nucleosome list to a 1D list
States = Nucs.flatten()
print(len(States))
print(len(Nucs))

colors = []
for elt in States:
    if elt == 0:
        colors.append('blue')
    elif elt == 1:
        colors.append('lightgrey')
    elif elt==2:
        colors.append('red')


#plotting
fig, ax = plt.subplots(figsize=(6,21))
#plt.title("beta= 0.01, sigma=1", fontsize = 25)
plt.scatter(x, y, c=colors, marker="o", lw=0.6)
plt.yticks([10,20,30,40,50,60,70,80], [410,420,430,440,450,460,470,480],fontsize=20)
plt.ylabel("Generations", fontsize =25)
plt.xticks(fontsize=20)
#plt.savefig("timecourse_203_new.pdf")