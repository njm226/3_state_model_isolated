
import numpy as np
cimport numpy as np
import random
import cython
from libc.math cimport log

@cython.boundscheck(False)
@cython.wraparound(False) 
@cython.cdivision(True)  
def t_loop(int duration, int[:] mt_region, int[:] positions, double[:] rates):
    
    # Very important function to create new random numbers for this cell to avoid dependencies!!!!
    np.random.seed()
    
    # total time (in units of generations) 
    cdef double T = 0 
    # time (set to 0 again after cell devides) 
    cdef double t1 = 0
    # time (set to 0 again after cell devides) 
    cdef double t2 = 0
    # duration of a simulation
    cdef double p1   
    cdef double p2 
     
    #sum of all rates
    cdef double total_rate = np.sum(rates)
    # array of the cumulative sums of all rates
    cdef double[:] cumsum_rates = np.cumsum(rates)
    #number of rates   
    cdef int n_rates = len(rates) 
   
    
    
    cdef int nn = len(mt_region) 
  
    cdef int TRE = 3
    

    cdef int Citrinel = 10
    cdef int Citriner = 11#13
    
    cdef int insulator = 14
    
    if nn == 20:
        
        mCherryl = 17#16
        mCherryr = 18#20
    
    elif nn == 21:
        
        mCherryl = 18#17
        mCherryr = 19#21 
        
    elif nn == 40:
        
        mCherryl = 37#36
        mCherryr = 38#40
        
        
    
    # silencing_threshold
    cdef int threshold1 = 1#2
    cdef int threshold2 = 1#2

    cdef int low_t_index
    cdef int pos_conv
    cdef int pos_rec
    cdef int i
    cdef int nuc_conv
    cdef int long_nn = 40#203
    cdef int nuc_rec 
    cdef double ran
    cdef int x

    
    cdef int n = 100000 # large number!
    
    # generates n random integers from 0 to size of positions
    cdef int[:] random_integers = np.random.randint(len(positions), size=n, dtype=np.int32)
    cdef int[:] random_integers2 = np.random.randint(len(positions), size=n, dtype=np.int32)
    # generates 4 * n different random numbers uniformly distributed between 0 and 1
    cdef double[:] random_doubles1 = np.random.random(n)
    cdef double[:] random_doubles2 = np.random.random(n)
    cdef double[:] random_doubles3 = np.random.random(n)
    cdef double[:] random_doubles4 = np.random.random(n)
    cdef double[:] random_doubles5 = np.random.random(n)
    
    #index for generating new random uniformly distributed numbers if running out of them

    cdef int j = 0
    cdef int k
    cdef int l = 0
    cdef double rand
    
        # cenH and EcoRV are not silent (0) at the beginning
    cdef int Citrine_silent = 0
    cdef int mCherry_silent = 0
    
    # cenH and EcoRV status at each time point
    cdef int[:] Citrine_status_list = np.zeros(duration*3, dtype=np.int32)
    cdef int[:] mCherry_status_list = np.zeros(duration*3, dtype=np.int32)

    # list to store the colorcoded nucleosome states of current mt_region
    states = []
    
    System_state = np.zeros([int(duration*3), nn])
    for i in xrange(nn):
                System_state[0,i] = mt_region[i]
    cdef int elt = 0
    
    S_nucleosomes = []
    A_nucleosomes = []
    U_nucleosomes = []
    S_nucleosomes_cenH = []
    

    current_states = []
    region_state = 0
    
    cdef int[:] mod_prob = np.zeros(240, dtype=np.int32)
    
    while T <= duration:
        j += 1
        
        
        #generating new random uniformly distributed numbers if running out of them
        if j >= n:
            # used for chosing a position
            random_integers = np.random.randint(len(positions), size=n, dtype=np.int32)
            random_integers2 = np.random.randint(len(positions), size=n, dtype=np.int32)
            #used for generating time_increase
            random_doubles1 = np.random.random(n)
            random_doubles2 = np.random.random(n)
            random_doubles3 = np.random.random(n)
            random_doubles4 = np.random.random(n)
            random_doubles5 = np.random.random(n)
            j = 0 

        # choses the time of the fastest reaction (equivalent to generating 8 different numbers)
        time_increase = -log(random_doubles1[j]) / total_rate
        
        # generates a random number within the range of total_rate (sum of all rates)
        rand = total_rate * random_doubles2[j] 

        low_t_index = n_rates - 1
        
        # goes through all rates 
        for k in range(n_rates):
            # if rand is smaller then the cumsum at pos of kth rate
            if rand < cumsum_rates[k]:
                # the index of the fastest rate is chosen
                low_t_index = k
                break
        
        # increases time by the length of the 
        T += time_increase
        t1 += time_increase
        t2 += time_increase
        

        # if the spontaneous conversion-rate (direct conversion of A to U) is chosen
        if low_t_index == 0:
            # a position of a nucleosome to be converted is chosen
            pos_conv = random_integers[j]
            # the nucleosome at that posion is selected
            nuc_conv = mt_region[pos_conv]
            
            if nuc_conv == 0:
                mt_region[pos_conv]=1
                
                
        # if the spontaneous conversion-rate (direct conversion of U to A) is chosen
        elif low_t_index == 1:            
           # a position of a nucleosome to be converted is chosen
            pos_conv = random_integers[j]
            # the nucleosome at that posion is selected
            nuc_conv = mt_region[pos_conv]
            
            if nuc_conv == 1:
                mt_region[pos_conv]=0
                
                
        # if the spontaneous conversion-rate (direct conversion of U to M) is chosen
        elif low_t_index == 2:   
            # a position of a nucleosome to be converted is chosen
            pos_conv = random_integers[j]
            # the nucleosome at that posion is selected  
            nuc_conv = mt_region[pos_conv]
    
            if nuc_conv == 1:
                mt_region[pos_conv]= 2
    
    
        # if the spontaneous conversion-rate (direct conversion of M to U) is chosen
        elif low_t_index == 3:      
            # a position of a nucleosome to be converted is chosen
            pos_conv = random_integers[j]
            # the nucleosome at that posion is selected
            nuc_conv = mt_region[pos_conv]
            
            if nuc_conv == 2:
                mt_region[pos_conv]=1 
                
                
         # if the spontaneous conversion-rate (direct conversion in special region of A to U) is chosen
        elif low_t_index == 4:
            
            #if mt_region[TRE] == 1:
            
            mt_region[TRE]=2 
          
  
            # else, nothing happens
                
            
        # if the global recruitment-rate M-catalysed change of U to M (recruited conversion) is chosen
        elif low_t_index == 5:                            
            
            # a position of a nucleosome to be converted is chosen
            pos_rec = random_integers[j]
            # the nucleosome at that posion is selected
            nuc_rec = mt_region[pos_rec]
            
            if random_doubles3[j] > 0.5:
                pos_conv = pos_rec - 1
            else:
                pos_conv = pos_rec + 1
                
                    
            if pos_conv >= 0 and pos_conv <= nn-1:
                # the nucleosome at that posion is selected
                nuc_conv = mt_region[pos_conv]
            else:
                nuc_conv = -1
            
            
            # if the recruiting nucleosome is in state M
            if nuc_rec == 2:
                # and the nucleosome to be converted is in state U
                if nuc_conv == 1:
                    # then it is changed to an A
                   mt_region[pos_conv]= 2
            
            
                   
        
        # if the local recruitment-rate A-catalysed change of M to U (recruited conversion) is chosen
        elif low_t_index == 6:
          # a position of a recruiting nucleosome is chosen
            pos_rec = random_integers[j]
            # the nucleosome at that posion is selected
            nuc_rec = mt_region[pos_rec]
            
            
            if random_doubles3[j] > 0.5:
                pos_conv = pos_rec - 1
            else:
                pos_conv = pos_rec + 1
                
                    
            if pos_conv >= 0 and pos_conv <= nn-1:
                # the nucleosome at that posion is selected
                nuc_conv = mt_region[pos_conv]
            else:
                nuc_conv = -10
            
            # if the recruiting nucleosome is in state A
            if nuc_rec == 0:
                # and the nucleosome to be converted is in state M
                if nuc_conv == 2:
                    # then it is changed to an U
                   mt_region[pos_conv]= 1
                   
                   
        # if the local recruitment-rate A-catalysed change of U to A (recruited conversion) is chosen         
        elif low_t_index == 7:            
            # a position of a recruiting nucleosome is chosen
            pos_rec = random_integers[j]
            # the nucleosome at that posion is selected
            nuc_rec = mt_region[pos_rec]
            
            
            if random_doubles3[j] > 0.5:
                pos_conv = pos_rec - 1
            else:
                pos_conv = pos_rec + 1
                
                    
            if pos_conv >= 0 and pos_conv <= nn-1:
                # the nucleosome at that posion is selected
                nuc_conv = mt_region[pos_conv]
            else:
                nuc_conv = -10 
            
            # if the recruiting nucleosome is in state A
            if nuc_rec == 0:
                # and the nucleosome to be converted is in state U
                if nuc_conv == 1:
                    # then it is changed to an A
                   mt_region[pos_conv]= 0 
                   
                    
        # if the local recruitment-rate (recruited conversion of A (0) to U (1)) is chosen         
        elif low_t_index == 8:  
            
            # a position of a recruiting nucleosome is chosen
            pos_rec = random_integers[j]
            # the nucleosome at that posion is selected
            nuc_rec = mt_region[pos_rec]
            
            if random_doubles3[j] > 0.5:
                pos_conv = pos_rec - 1
            else:
                pos_conv = pos_rec + 1
                
                    
            if pos_conv >= 0 and pos_conv <= nn-1:
                # the nucleosome at that posion is selected
                nuc_conv = mt_region[pos_conv]
            else:
                nuc_conv = -10
            
                
            # and if the recruiting nucleosome is in state M (2)
            if nuc_rec == 2:
                # and the nucleosme to be converted is in state A
                if nuc_conv == 0:
                    # then the nucleosome to be converted is changed to a U
                    mt_region[pos_conv] = 1
                    
                    
        elif low_t_index == 9:
            
            if mt_region[insulator] == 1:
                mt_region[insulator]=0
        
        if t1 >= 0.333333333333333333333333333333333:   
            t1=0
            
            # if t2 >= 20:
            #     #the state of the mt_region at this time point is stored 
            #     #mt_matrix[m]=mt_region
            #     for i in positions:
            #         rand = random.choice([0,1])
            #         if rand == 1:
            #             mt_region[i]=1
              
            #     t2=0
                
            elt += 1
            for i in xrange(nn):
                System_state[elt,i] = mt_region[i]
                
            Citrine_blue = np.count_nonzero(System_state[elt,Citrinel:Citriner+1]== 0)
            Citrine_red =  np.count_nonzero(System_state[elt,Citrinel:Citriner+1] == 2)
            
            mCherry_blue =  np.count_nonzero(System_state[elt,mCherryl:mCherryr+1] == 0)
            mCherry_red = np.count_nonzero(System_state[elt,mCherryl:mCherryr+1] == 2)
            
            if Citrine_red>= threshold1: #- Citrine_blue 
                Citrine_silent = 1
            else:
                Citrine_silent = 0
                
            if mCherry_red  >= threshold2:  #- mCherry_blue
                mCherry_silent = 1
            else:
                mCherry_silent = 0 
        
                    
            Citrine_status_list[elt]=Citrine_silent
            mCherry_status_list[elt]=mCherry_silent 
                
            # if t2 >= 20:
            #   #the state of the mt_region at this time point is stored 
            #   #mt_matrix[m]=mt_region
            #   for i in positions:
            #       rand = random.choice([0,1])
            #       if rand == 1:
            #           mt_region[i]=1
            
            #   t2=0
              
    
    
    return  Citrine_status_list, mCherry_status_list, states, System_state
    
    

