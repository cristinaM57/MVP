import matplotlib
matplotlib.use('TKAgg')

import sys
import math
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import csv

def writefile(lx, model):
    outfile = f"{model}_Size{lx}.csv"
    header = ['Temperature', 'Magnetisation', 'Susceptibility', 'Energy', 'HeatCapacity', 'mag_error', 'sus_error', 'eng_error', 'shc_error']
    with open(outfile, "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)

def Glauber(lx, ly, kT, spin, model):

    J=1.0
    nstep=10100 #discard first 100 sweeps
    sweeps=0

    """
    fig = plt.figure()
    im=plt.imshow(spin, animated=True)
    """

    #calculate initial energy
    E_initial = 0
    for i in range(lx):
        for j in range(ly):
            E = ((-1)*J*spin[i,j]*(spin[(i+1)%lx,j] + spin[(i-1)%lx,j] + spin[i,(j+1)%lx] + spin[i,(j-1)%lx]))/2 #account for counting twice
            E_initial += E

    #create lists to store data
    magnetisation = []
    energies = []

    #update loop here - for Glauber dynamics

    for n in range(nstep):
        #create lists with 2500 random numbers
        i_values = np.random.choice(lx, lx*ly)
        j_values = np.random.choice(lx, lx*ly)
        
        for i in range(lx*ly):

    #select spin randomly
            itrial = i_values[i]
            jtrial = j_values[i]

    #compute delta E eg via function (account for periodic BC by using mod length)

            delta_E = 2*J*spin[itrial,jtrial]*(spin[(itrial+1)%lx,jtrial] + spin[(itrial-1)%lx,jtrial] + spin[itrial,(jtrial+1)%lx] + spin[itrial,(jtrial-1)%lx])

    #perform metropolis test

            if delta_E <= 0 or random.random() <= math.exp(-1*delta_E/kT):
                spin[itrial,jtrial] *= -1  #you only need to change one spin, not update the whole thing, but not sure how
            else:
                delta_E = 0
            
            # update initial energy
            E_initial += delta_E
    
        #occasionally plot or update measurements, eg every 10 sweeps
        if (n%10==0): 

            #calculate and save data every 10 sweeps
            
            M = np.sum(spin)
            magnetisation.append(M)
            energies.append(E_initial)



        #       update measurements
        #       dump outputs
            """
            f=open('spins.dat','w')
            for i in range(lx*ly):
                for j in range(ly):
                    f.write('%d %d %lf\n'%(i,j,spin[i,j]))
            f.close()
            """
            """
        #       show animation
            plt.cla()
            im=plt.imshow(spin, animated=True)
            plt.draw()
            plt.pause(0.0001)
            """

            sweeps += 10
            print("temperature:" + str(kT) + " sweeps:" + str(sweeps), end="\r")

    #save measurements to file after 10000 sweeps (disregard first 100)
    n_del = 10
    del magnetisation[:n_del]
    del energies[:n_del]

    mag_array = np.asarray(magnetisation)
    eng_array = np.asarray(energies)

    mag = abs(np.mean(magnetisation))
    sus = (1/((lx*ly)*kT))*np.var(mag_array)
    eng = np.mean(energies)
    shc = (1/((lx*ly)*(kT**2)))*np.var(eng_array)

    mag_err, sus_err, eng_err, shc_err = bootstrap(mag_array, eng_array, lx, ly, kT)

    outfile = f"{model}_Size{lx}.csv"
    data = [kT, mag, sus, eng, shc, mag_err, sus_err, eng_err, shc_err]
    with open(outfile, "a", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(data)

    return spin

def Kawasaki(lx, ly, kT, spin, model):
    J=1.0
    nstep=10100
    sweeps = 0

    """
    fig = plt.figure()
    im=plt.imshow(spin, animated=True)
    """

    #calculate initial energy
    E_initial = 0
    for i in range(lx):
        for j in range(ly):
            E = ((-1)*J*spin[i,j]*(spin[(i+1)%lx,j] + spin[(i-1)%lx,j] + spin[i,(j+1)%lx] + spin[i,(j-1)%lx]))/2 #account for counting twice
            E_initial += E

    #create lists to store data
    magnetisation = []
    energies = []

    #update loop here - for Kawasaki dynamics

    for n in range(nstep):
        #create lists with 2500 random numbers
        i_values1 = np.random.choice(lx, lx*ly)
        j_values1 = np.random.choice(lx, lx*ly)
        i_values2 = np.random.choice(lx, lx*ly)
        j_values2 = np.random.choice(lx, lx*ly)

        for i in range(lx*ly):
        #select two spins randomly
            itrial1 = i_values1[i]
            jtrial1 = j_values1[i]
            spin1 = spin[itrial1,jtrial1]

            itrial2 = i_values2[i]
            jtrial2 = j_values2[i]
            spin2 = spin[itrial2,jtrial2]

            #if same spin selected twice, change both spins until they are different
            while itrial1 == itrial2 & jtrial1 == jtrial2: 
                itrial2 = np.random.randint(0,lx)
                jtrial2 = np.random.randint(0,ly)

            #if spins are the same, energy change is zero
            if spin1 == spin2:
                continue
            
            else:
                #compute delta E eg via function (account for periodic BC by using mod length)
                Glauber_1 = 2*J*spin[itrial1,jtrial1]*(spin[(itrial1+1)%lx,jtrial1] + spin[(itrial1-1)%lx,jtrial1] + spin[itrial1,(jtrial1+1)%lx] + spin[itrial1,(jtrial1-1)%lx])
                Glauber_2 = 2*J*spin[itrial2,jtrial2]*(spin[(itrial2+1)%lx,jtrial2] + spin[(itrial2-1)%lx,jtrial2] + spin[itrial2,(jtrial2+1)%lx] + spin[itrial2,(jtrial2-1)%lx])
                delta_E = Glauber_1 + Glauber_2

                if (itrial1 == itrial2 & (jtrial1 == (jtrial2 - 1)%lx or jtrial1 == (jtrial2 + 1)%lx)) or (jtrial1 ==  jtrial2 & (itrial1 == (itrial2 - 1)%lx or itrial1 == (itrial2 + 1)%lx)):
                    delta_E += 4

        #perform metropolis test

            if delta_E <= 0 or random.random() <= math.exp(-1*delta_E/kT):
                #swap spins (this is the same as changing the sign)                                
                spin[itrial1,jtrial1] *= -1
                spin[itrial2,jtrial2] *= -1
            else:
                delta_E = 0

            # update initial energy
            E_initial += delta_E

        #occasionally plot or update measurements, eg every 10 sweeps
        if(n%10==0): 

            #calculate and save data every 10 sweeps
                
            M = np.sum(spin)
            magnetisation.append(M)
            energies.append(E_initial)

        #       update measurements
        #       dump output
            """
            f=open('spins.dat','w')
            for i in range(lx):
                for j in range(ly):
                    f.write('%d %d %lf\n'%(i,j,spin[i,j]))
            f.close()
            """
            """
        #   show animation
            plt.cla()
            im=plt.imshow(spin, animated=True)
            plt.draw()
            plt.pause(0.0001)
            """

            sweeps += 10
            print("temperature:" + str(kT) + " sweeps:" + str(sweeps), end="\r")

    #save measurements to file after 10000 sweeps (disregard first 100)
    n_del = 10
    del magnetisation[:n_del]
    del energies[:n_del]

    mag_array = np.asarray(magnetisation)
    eng_array = np.asarray(energies)

    mag = abs(np.mean(magnetisation))
    sus = (1/((lx*ly)*kT))*np.var(mag_array)
    eng = np.mean(energies)
    shc = (1/((lx*ly)*(kT**2)))*np.var(eng_array)

    mag_err, sus_err, eng_err, shc_err = bootstrap(mag_array, eng_array, lx, ly, kT)

    outfile = f"{model}_Size{lx}.csv"
    data = [kT, mag, sus, eng, shc, mag_err, sus_err, eng_err, shc_err]
    with open(outfile, "a", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(data)

    return spin

def bootstrap(magnetisation, energy, lx, ly, kT):
    L = len(energy)
    nsubsets = 6

    shc_subset  = []
    sus_subset = []

    for i in range(nsubsets):
        indices = np.random.choice(L, int(L/nsubsets))
        eng_subset = energy[indices]
        mag_subset = magnetisation[indices]

        sus_subset.append((1/((lx*ly)*kT))*np.var(mag_subset))
        shc_subset.append((1/((lx*ly)*(kT**2)))*np.var(eng_subset))

    return np.std(magnetisation), np.std(sus_subset), np.std(energy), np.std(shc_subset)

def main():
    
    T_range = np.arange(1,3.1, 0.1)
    temperatures = []
    for i in range (len(T_range)):
        temperatures.append(np.round(T_range[i], 2))  
    
    #input

    if(len(sys.argv) != 3):
        print ("Usage python ising.animation.py N T")
        sys.exit()

    lx=int(sys.argv[1]) 
    ly=lx 
    model = str(sys.argv[2])

    writefile(lx, model)
    
    if model == "Glauber":
    #initialise all in 1 and then discard first 100 sweeps
        spin = np.random.choice([1], [lx, ly])

        for i in range (len(T_range)):
            updated_spin = Glauber(lx, ly, temperatures[i], spin, model)
            spin = updated_spin
    
    elif model == "Kawasaki":
    #initialise spins randomly, discard first 100 sweeps
        spin = np.random.choice([-1,1], [lx, ly])
        
        for i in range (len(T_range)):
            updated_spin = Kawasaki(lx, ly, temperatures[i], spin, model)
            print("Prog")
            spin = updated_spin

main()