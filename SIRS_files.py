import matplotlib
matplotlib.use('TKAgg')

import sys
import math
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd
import csv

def writefile(lx, calc):
    outfile = f"{calc}_Size{lx}.csv"
    header = ['Infected_Fraction', 'Immune_Percent']
    with open(outfile, "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)

def bootstrap(infected_array, lx, ly):
    L = len(infected_array)
    nsubsets = 1000

    infected_variances  = []

    for i in range(nsubsets):
        indices = np.random.randint(0,L, L)
        subset_infected = infected_array[indices]

        infected_variances.append(np.var(subset_infected)/(lx*ly))

    return np.std(infected_variances)

def SIRS_phases(lx, ly, p1, p2, p3, spin):

    nstep=1200 #discard first 200 sweeps
    sweeps=0
    infected = []
    
    """
    fig = plt.figure()
    im=plt.imshow(spin, animated=True)
    """
    #update loop here - for SIRS

    for n in range(nstep):

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
            I = np.count_nonzero(spin == -1)
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
    n_del = 200
    del infected[:n_del]

    inf_array = np.asarray(infected)
    inf_avg = np.mean(infected)
    inf_var = np.var(inf_array)/(lx*ly)
    inf_frac = inf_avg/(lx*ly) 

    return inf_frac, inf_var

def SIRS_var(lx, ly, p1, p2, p3, spin):

    nstep=10500 #discard first 500 sweeps
    sweeps=0
    infected = []
    
    """
    fig = plt.figure()
    im=plt.imshow(spin, animated=True)
    """
    #update loop here - for SIRS

    for n in range(nstep):
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
            I = np.count_nonzero(spin == -1)
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
    n_del = 500
    del infected[:n_del]

    inf_array = np.asarray(infected)
    inf_avg = np.mean(infected)
    inf_var = np.var(inf_array)/(lx*ly)
    inf_frac = inf_avg/(lx*ly) 

    var_err = bootstrap(inf_array, lx, ly)

    return inf_frac, inf_var, var_err

def immune(spin, lx, ly, per_immune):

    total_people = lx*ly
    immune_people = int(per_immune*total_people)

    #create lists with 2500 random numbers
    i_values = np.random.choice(lx, lx*ly)
    j_values = np.random.choice(ly, lx*ly)
    
    count = np.count_nonzero(spin == 2)
    while count <= immune_people:
        
        itrial = np.random.randint(0,lx)
        jtrial = np.random.randint(0,ly)

        if spin[itrial, jtrial] != 2:
            spin[itrial, jtrial] = 2
        
        count = np.count_nonzero(spin == 2)

    return spin

def SIRS_immunity(lx, ly, p1, p2, p3, spin):

    nstep=10500 #discard first 500 sweeps
    sweeps=0
    infected = []
    
    """
    fig = plt.figure()
    im=plt.imshow(spin, animated=True)
    """
    #update loop here - for SIRS

    for n in range(nstep):

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
            I = np.count_nonzero(spin == -1)
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
    n_del = 500
    del infected[:n_del]

    inf_array = np.asarray(infected)
    inf_avg = np.mean(infected)
    inf_var = np.var(inf_array)/(lx*ly)
    inf_frac = inf_avg/(lx*ly) 

    return inf_frac, inf_var

def main():
    
    #input

    lx=int(sys.argv[1]) 
    ly=lx
    calc = str(sys.argv[2]) 

    if calc == "heatmap":

        p1_arr = np.arange(0, 1, 0.05)
        p2 = 0.5
        p3_arr = np.arange(0, 1, 0.05)

        size = len(p1_arr)
        
        lattice = np.zeros([size, size])
        var_lattice = np.zeros([size, size])

        for i in range(size):
            for j in range(size):
                
                spin = np.random.choice([-1, 0, 1], [lx, ly])
                p1 = p1_arr[i]
                p3 = p3_arr[j]
                lattice[i,j] = SIRS_phases(lx, ly, p1, p2, p3, spin)[0]
                #var_lattice[i,j] = SIRS_phases(lx, ly, p1, p2, p3, spin)[1]

        df = pd.DataFrame(data=lattice.astype(float))
        df.to_csv(f'SIRS_heatmap_size{lx}', sep=' ', header = False, float_format='%.5f', index=False)
        #df1 = pd.DataFrame(data=var_lattice.astype(float))
        #df1.to_csv(f'SIRS_heatmap_var_size{lx}', sep=' ', header = False, float_format='%.5f', index=False)

    elif calc == "variance":

        p1_arr = np.arange(0.2, 0.50, 0.01)
        p2 = 0.5
        p3 = 0.5

        size = len(p1_arr)
        
        average_fraction = []
        var_1D = []
        var_err_1D = []
            
        for i in range(size):
            
            spin = np.random.choice([-1, 0, 1], [lx, ly])
            p1 = p1_arr[i]
                
            #average_fraction.append(SIRS_phases(lx, ly, p1, p2, p3, spin)[0])  
            var_1D.append(SIRS_var(lx, ly, p1, p2, p3, spin)[1])
            var_err_1D.append(SIRS_var(lx, ly, p1, p2, p3, spin)[2])
        
        
        #df = pd.DataFrame(data=average_fraction.astype(float))
        #df.to_csv(f'SIRSVAR_size{lx}', sep=' ', header = False, float_format='%.5f', index=False)
        df1 = pd.DataFrame(data=var_1D)
        df1.to_csv(f'SIRSVAR_var_size{lx}', sep=' ', header = False, float_format='%.5f', index=False)
        df2 = pd.DataFrame(data=var_err_1D)
        df2.to_csv(f'SIRSVAR_var_err_size{lx}', sep=' ', header = False, float_format='%.5f', index=False)
    
    elif calc == "immunity":

        writefile(lx, calc)

        p1 = 0.5
        p2 = 0.5
        p3 = 0.5

        probs = np.arange(0, 1, 0.01)
        spin = np.random.choice([-1, 0, 1], [lx, ly])

        for i in range(len(probs)):
            prob = probs[i]
            spin_update = immune(spin, lx, ly, prob)
            x = (SIRS_immunity(lx, ly, p1, p2, p3, spin_update)[0])
            outfile = f"{calc}_Size{lx}.csv"
            data = [x, prob]
            with open(outfile, "a", newline='') as f:
                writer = csv.writer(f)
                writer.writerow(data)

main()