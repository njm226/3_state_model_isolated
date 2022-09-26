
import numpy as np
import pyximport
pyximport.install(reload_support=True, setup_args={"include_dirs":np.get_include()})
from all_local_40 import t_loop
import all_local_40
import importlib
importlib.reload(all_local_40)
import pandas as pd
import matplotlib.pyplot as plt

def simple(X_Y):  
    
    
    N=X_Y[0] 
    
    #mating type region (array of 140 nucleosomes) starting all with silent nucs (2)
    mt_region = np.ones(N, dtype=np.int32)*X_Y[1]
    #indices of the mt_region corresponding to positions of nucleosomes
    positions = np.arange(len(mt_region), dtype=np.int32)
    

    duration = X_Y[2]


     
    
    direct = X_Y[3]
    #direct = X_Y[3] 
    
   
   
    
    # local recruitment-rate M-catalysed change of U to M (recruited conversion)
    alpha1 = 250*len(mt_region)#*0.05 
    # local recruitment-rate A-catalysed change of M to U (recruited conversion)
    alpha2 = 250*len(mt_region)#*97
    # local recruitment-rate A-catalysed change of U to A (recruited conversion)
    alpha3 = 250*len(mt_region)#*97
    # global recruitment-rate (recruited conversion of A (0) to U (1))
    alpha4 = 250*len(mt_region)#*0.05
    # spontaneous conversion-rate (direct conversion of A to U)
    beta1 = direct*len(mt_region)#*0.05 
    # spontaneous conversion-rate (direct conversion)
    beta2 = direct*len(mt_region)#*0.05
    # spontaneous conversion-rate (direct conversion)
    beta3 = direct*len(mt_region)#*0.05 
    # spontaneous conversion-rate (direct conversion)
    beta4 = direct*len(mt_region)#*0.05
    # spontaneous conversion-rate in cenH region (only A to U)
    beta5 = X_Y[4]

    if N == 21:
        beta6 = 25#50#*len(mt_region)#   
        beta5 = 0#1000
   
    if N == 21:
        rates = np.array([beta1, beta2, beta3, beta4, beta5, alpha1, alpha2, alpha3, alpha4, beta6], dtype=np.double)
    else: 
        rates = np.array([beta1, beta2, beta3, beta4, beta5, alpha1, alpha2, alpha3, alpha4], dtype=np.double)
    #print(cenH_status_list)
        
    print(X_Y)
    
    Citrine_status_list, mCherry_status_list, states, System_state= t_loop(duration, mt_region, positions, rates)#Citrine_status_list, mCherry_status_list 
    
    
    
    # #print(nucs[:100])
    # positions = list(range(len(mt_region)))
    # x = positions*(80)
    # y= [i for i in range(80) for _ in range(len(mt_region))]
    
    
    # # show every 200th element
    
    # Nucs = System_state[0:240:3]
    # #Nucs = System_state[720:960:3]
    
    # # convert 2D nucleosome list to a 1D list
    # States = Nucs.flatten()
    # print(len(States))
    # print(len(Nucs))
    
    # colors = []
    # for elt in States:
    #     if elt == 0:
    #         colors.append('blue')
    #     elif elt == 1:
    #         colors.append('lightgrey')
    #     elif elt==2:
    #         colors.append('red')
    
    
    # #plotting
    # fig, ax = plt.subplots(figsize=(6,21))
    # #plt.title("beta= 0.01, sigma=1", fontsize = 25)
    # plt.scatter(x, y, c=colors, marker="o", lw=0.6)
    # plt.yticks(fontsize=20)
    # plt.ylabel("Generations", fontsize =25)
    # plt.xticks(fontsize=20)
    # #plt.savefig("timecourse_203_new.pdf")
        
    

    return list(Citrine_status_list), list(mCherry_status_list) 

# # #values for full(special set to 300)
# if __name__ == '__main__':
#     import time
#     #import cProfile
#     t1 = time.time()
#     simple([40, 0, 80, 1, 1000])
#     print(time.time() - t1)


