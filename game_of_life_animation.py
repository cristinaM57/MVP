import matplotlib
matplotlib.use('TKAgg')

import sys
import math
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import csv
import time

def writefile(lx, initial):
    outfile = f"{initial}_Size{lx}.csv"
    header = ['Run', 'Active Sites', 'Time', 'Sweeps']
    with open(outfile, "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)

def game_of_life(lx, ly, spin, initial, run, nsims):

    start_time = time.time()

    nstep=10000 #discard first 100 sweeps
    sweeps=0
    
    """
    fig = plt.figure()
    im=plt.imshow(spin, animated=True)
    """

    #update loop here - for game of life

    for n in range(nstep):
        copy = np.copy(spin)
        active_sites = np.sum(copy)
        
        for i in range(lx):
            for j in range(ly):
                #check neighbours
                central = spin[i,j]
                N = spin[(i-1)%lx,j]
                E = spin[i,(j+1)%ly]
                S = spin[(i+1)%lx,j]
                W = spin[i,(j-1)%ly]
                NE = spin[(i-1)%lx,(j+1)%ly]
                SE = spin[(i+1)%lx,(j+1)%ly]
                SW = spin[(i+1)%lx,(j-1)%ly]
                NW = spin[(i-1)%lx,(j-1)%ly]
                sum_neighbours = N + E + S + W + NE + SE + SW + NW
                
                # if cell is alive
                if central == 1:
                    if sum_neighbours < 2:
                        copy[i,j] = 0
                    elif sum_neighbours > 3:
                        copy[i,j] = 0
                # if cell is dead
                else:
                    if sum_neighbours == 3:
                        copy[i,j] = 1
        
        spin = np.copy(copy)
        end_time = time.time()
        time_passed = end_time - start_time 
        
        if (n%1==0):
    
            outfile = f"{initial}_Size{lx}.csv"
            data = [run, active_sites, np.round(time_passed, 2), sweeps]
            with open(outfile, "a", newline='') as f:
                writer = csv.writer(f)
                writer.writerow(data)
        
            """
            #show animation
            plt.cla()
            im=plt.imshow(spin, animated=True)
            plt.draw()
            plt.pause(0.0001)
            """

            sweeps += 1
            print("simulation:" + str(run) + "/" + str(nsims) + " sweeps:" + str(sweeps) + "/" + str(nstep), end="\r")

def main():

    if(len(sys.argv) != 3):
        print ("Usage python name_of_code.py N model_type")
        sys.exit()

    lx=int(sys.argv[1]) 
    ly=lx 
    initial = str(sys.argv[2])
    nsims = 1000 #number of simulations

    writefile(lx, initial)
    
    if initial == "random":
        
        for i in range (nsims):
            run = i
            spin = np.random.choice([0,1], [lx, ly])
            game_of_life(lx, ly, spin, initial, run, nsims)

main()