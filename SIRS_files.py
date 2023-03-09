import matplotlib
matplotlib.use('TKAgg')

import sys
import math
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd

def SIRS_animation(lx, ly, p1, p2, p3, spin):

    nstep=1100 #discard first 100 sweeps
    sweeps=0
    infected = []
    
    """
    fig = plt.figure()
    im=plt.imshow(spin, animated=True)
    """
    #update loop here - for SIRS

    for n in range(nstep):

        I = 0

        #create lists with 2500 random numbers
        i_values = np.random.choice(lx, lx*ly)
        j_values = np.random.choice(ly, lx*ly)
        
        for i in range(lx*ly):          
    #select spin randomly
            itrial = i_values[i]
            jtrial = j_values[i]

    #perform test
            if spin[itrial, jtrial] == 0:
                # check neighbours
                if ((spin[(itrial+1)%lx,jtrial]) == -1) or ((spin[(itrial-1)%lx,jtrial]) == -1) or ((spin[itrial,(jtrial+1)%lx]) == -1) or ((spin[itrial,(jtrial-1)%lx]) == -1):
                    x = random.random()
                    #print("p1"+ str(x)) #debug
                    if x <= p1:
                        spin[itrial,jtrial] = -1  #you only need to change one spin, not update the whole thing
                else:
                    spin[itrial, jtrial] = 0
            
            elif spin[itrial, jtrial] == -1:
                I += 1
                x = random.random()
                #print("p2"+ str(x)) #debug
                if x <= p2:
                    spin[itrial,jtrial] = 1

            elif spin[itrial, jtrial] == 1:
                x = random.random()
                #print("p3"+ str(x)) #debug
                if x <= p3:
                    spin[itrial,jtrial] = 0
        
            #print("i" + str(itrial) + " j" + str(jtrial)) #debug

        #print(spin) #debug
                
        #occasionally plot or update measurements, eg every 10 sweeps
        if (n%1==0): 

            #calculate and save data every sweep
            infected.append(I)
            
            """
        #       show animation
            plt.cla()
            im=plt.imshow(spin, animated=True, vmin = -1, vmax =1)
            plt.draw()
            plt.pause(0.0001)
            """

            sweeps += 1
            print(" sweeps:" + str(sweeps) + "/" + str(nstep), end="\r")

    
    #save measurements to file after 1000 sweeps (disregard first 100)
    n_del = 100
    del infected[:n_del]

    inf_array = np.asarray(infected)
    inf_avg = np.mean(infected)
    inf_var = np.var(inf_array)
    inf_frac = inf_avg/(lx*ly) 

    return inf_frac, inf_var


def main():
    
    #input

    if(len(sys.argv) != 2):
        print ("Usage python ising.animation.py N p1 p2 p3")
        sys.exit()

    lx=int(sys.argv[1]) 
    ly=lx 
    
    p1_arr = np.arange(0, 1.05, 0.05)
    p2 = 0.5
    p3_arr = np.arange(0, 1.05, 0.05)

    size = len(p1_arr)
    
    lattice = np.zeros([size, size])
    var_lattice = np.zeros([size, size])


    for i in range(size):
        for j in range(size):
            
            spin = np.random.choice([-1, 0, 1], [lx, ly])

            p1 = p1_arr[i]
            p3 = p3_arr[j]
            
            lattice[i,j] = SIRS_animation(lx, ly, p1, p2, p3, spin)[0]  
            var_lattice[i,j] = SIRS_animation(lx, ly, p1, p2, p3, spin)[1]

    df1 = pd.DataFrame(lattice) 
    df1.to_csv("infected_fraction.csv")

    df2 = pd.DataFrame(lattice) 
    df2.to_csv("infected_fraction_var.csv")

main()