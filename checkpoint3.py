import numpy as np
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
import csv
import math

def charge_density(lx, ly, lz):
    rho = np.random.choice([0], (lx, ly, lz))
    lx_centre = int(lx/2)
    ly_centre = int(ly/2)
    lz_centre = int(lz/2)
    rho[lx_centre, ly_centre, lz_centre] = 1
    return rho


def magnetisation(lx, ly, lz):
    mag = np.random.choice([0], (lx, ly, lz))
    lx_centre = int(lx/2)
    ly_centre = int(ly/2)
    mag[lx_centre, ly_centre] = 1
    return mag
    

def chemical_potential(grid, a, k, dx):
    b = k*(1/(dx**2))
    laplacian_1 = (np.roll(grid, 1, axis=0) + np.roll(grid, -1, axis =0) + np.roll(grid, 1, axis =1) + np.roll(grid, -1, axis =1) - np.multiply(grid,4))
    x = (-1*a*grid) + (a)*(np.power(grid, 3)) - (b*laplacian_1)
    return x

def writefile_cahnhilliard(lx, phi):
    outfile = f"Cahn_Hilliard_Size{lx}_phi{phi}.csv"
    header = ['Time', 'Free_Energy_Density']
    with open(outfile, "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)


def free_energy_density(grid, a, k, dx):
    grad_y = (np.roll(grid, -1, axis=1) - np.roll(grid, 1, axis=1)) / (2*dx)
    grad_x = (np.roll(grid, -1, axis=0) - np.roll(grid, 1, axis=0)) / (2*dx)
    grad_sqr = np.power(grad_y,2) + np.power(grad_x,2)
    return (((-1)*(a/2)*np.power(grid, 2)) + ((a/4)*np.power(grid, 4)) + ((k/2)*grad_sqr))


def cahn_hilliard_anim(phi, lx, ly, nstep, dt, dx):
    
    offset = 0.4

    phi_zero = np.random.uniform(phi - offset, phi + offset, (lx, ly))
    a = 0.1
    k = 0.1
    M = 0.1
    c = (1/(dx**2))

    #define the first grid
    phi_zero_old = phi_zero

    fig = plt.figure()
    im=plt.imshow(phi_zero_old, interpolation='gaussian', animated=True)

    sweeps = 0

    for n in range(nstep):
        
        #obtain chemical potential
        grid = chemical_potential(phi_zero_old, a, k, dx)

        #laplacian of chemical potential
        laplacian_2 = (np.roll(grid, 1, axis=0) + np.roll(grid, -1, axis =0) + np.roll(grid, 1, axis =1) + np.roll(grid, -1, axis =1) - np.multiply(grid,4))
        
        #calculate the new phi 
        phi_zero_old += M*dt*c*laplacian_2

        if (n%100==0):
            
            #show animation
            plt.cla()
            im=plt.imshow(phi_zero_old, interpolation='gaussian', animated=True, cmap = 'PiYG')
            plt.draw()
            plt.pause(0.0001)
            

            sweeps += 100
            print(" sweeps:" + str(sweeps) + "/" + str(nstep), end="\r")


def cahn_hilliard_file(phi, lx, ly, nstep, dt, dx):
    
    outfile = f"Cahn_Hilliard_Size{lx}_phi{phi}.csv"
    offset = 0.4

    phi_zero = np.random.uniform(phi - offset, phi + offset, (lx, ly))
    a = 0.1
    k = 0.1
    M = 0.1
    c = (1/(dx**2))

    #define the first grid
    phi_zero_old = phi_zero

    sweeps = 0

    for n in range(nstep):
        
        #obtain chemical potential
        grid = chemical_potential(phi_zero_old, a, k, dx)

        #laplacian of chemical potential
        laplacian_2 = (np.roll(grid, 1, axis=0) + np.roll(grid, -1, axis =0) + np.roll(grid, 1, axis =1) + np.roll(grid, -1, axis =1) - np.multiply(grid,4))
        
        #calculate the new phi 
        phi_zero_old += M*dt*c*laplacian_2

        if (n%1==0):

            energy_density_arr = free_energy_density(phi_zero_old, a, k, dx)
            energy_density = np.sum(energy_density_arr)
            time = n*dt

            data = [time, energy_density]
            with open(outfile, "a", newline='') as f:
                writer = csv.writer(f)
                writer.writerow(data)
            

            sweeps += 1
            print(" sweeps:" + str(sweeps) + "/" + str(nstep), end="\r")


def writefile_jacobi(lx, condition):
    outfile = f"Jacobi_Size{lx}_{condition}.csv"
    header = ['Convergence_Sweeps']
    with open(outfile, "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)

# working
def jacobi_pad(lx, ly, lz, dt, dx, condition):
    
    outfile = f"Jacobi_Size{lx}_{condition}.csv"

    phi_zero = np.zeros((lx, ly, lz))
    
    # plot
    size_added = lx + 2
    size_small = lx // 2

    
    if condition == "electric":
        #define the point charge
        charge = charge_density(lx, ly, lz)

        #pad
        phi_zero[[0,-1], :, :]=0
        phi_zero[:, [0,-1], :]=0
        phi_zero[:, :, [0,-1]]=0
        phi_zero_old = phi_zero

        sweeps = 0

        while True:
        
            #calculate the new phi
            phi_zero_update = (1/6)*(np.roll(phi_zero_old, 1, axis=0) + np.roll(phi_zero_old, -1, axis=0) + np.roll(phi_zero_old, 1, axis=1) + np.roll(phi_zero_old, -1, axis=1) + np.roll(phi_zero_old, 1, axis=2) + np.roll(phi_zero_old, -1, axis=2) + (dx**2)*charge)
            phi_zero_update[[0,-1], :, :]=0
            phi_zero_update[:, [0,-1], :]=0
            phi_zero_update[:, :, [0,-1]]=0
        
            #difference between new and old     
            #takes x sweeps #need to reset this
            dist = np.mean(np.abs(phi_zero_update-phi_zero_old))
            print(dist)
            if np.isclose(dist,0,atol=5e-8):
                break   
        

            #update phi_zero
            phi_zero_old = np.copy(phi_zero_update)  

            sweeps += 1
            print(" sweeps:" + str(sweeps), end="\r")

        data = [sweeps]
        with open(outfile, "a", newline='') as f:
            writer = csv.writer(f)
            writer.writerow(data)

        # reshaping the array from 3D
        # matrice to 2D matrice.
        arr_reshaped = phi_zero_old.reshape(phi_zero_old.shape[0], -1)
  
        # saving reshaped array to file.
        np.savetxt(f"Jacobi_Size{lx}_{condition}.txt", arr_reshaped)

        """
        #plt.imshow(phi_zero_old[0:size_added,size_small,0:size_added],interpolation='gaussian')
        plt.imshow(phi_zero_old[0:size_added,size_small,0:size_added])
        plt.colorbar()
        plt.show()
        """

    elif condition == "magnetic":
        #define the magnetic plane
        charge = magnetisation(lx, ly, lz)

        #pad
        phi_zero[[0,-1], :, :]=0
        phi_zero[:, [0,-1], :]=0
        phi_zero_old = phi_zero

        sweeps = 0

        while True:
        
            #calculate the new phi
            phi_zero_update = (1/6)*(np.roll(phi_zero_old, 1, axis=0) + np.roll(phi_zero_old, -1, axis=0) + np.roll(phi_zero_old, 1, axis=1) + np.roll(phi_zero_old, -1, axis=1) + np.roll(phi_zero_old, 1, axis=2) + np.roll(phi_zero_old, -1, axis=2) + (dx**2)*charge)
            phi_zero_update[[0,-1], :, :]=0
            phi_zero_update[:, [0,-1], :]=0
        
            #difference between new and old     
            #takes x sweeps #need to reset this
            dist = np.mean(np.abs(phi_zero_update-phi_zero_old))
            print(dist)
            if np.isclose(dist,0,atol=5e-8):
                break   
        

            #update phi_zero
            phi_zero_old = np.copy(phi_zero_update)  

            sweeps += 1
            print(" sweeps:" + str(sweeps), end="\r")

        data = [sweeps]
        with open(outfile, "a", newline='') as f:
            writer = csv.writer(f)
            writer.writerow(data)

        # reshaping the array from 3D
        # matrice to 2D matrice.
        arr_reshaped = phi_zero_old.reshape(phi_zero_old.shape[0], -1)
  
        # saving reshaped array to file.
        np.savetxt(f"Jacobi_Size{lx}_{condition}.txt", arr_reshaped)

        """
        #plt.imshow(phi_zero_old[0:size_added,size_small,0:size_added],interpolation='gaussian')
        plt.imshow(phi_zero_old[0:size_added,size_small,0:size_added])
        plt.colorbar()
        plt.show()
        """
          
    

def writefile_gaussseidel(lx, condition):
    outfile = f"Gauss_Seidel_Size{lx}_{condition}.csv"
    header = ['Convergence_Sweeps']
    with open(outfile, "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)

# working
def gauss_seidel_pad(lx, ly, lz, dt, dx, condition):
    
    outfile = f"Gauss_Seidel_Size{lx}_{condition}.csv"
    
    phi_zero = np.zeros((lx, ly, lz))

    # plot
    size_added = lx + 2
    size_small = lx // 2

    if condition == 'electric':

        #pad
        phi_zero[[0,-1], :, :]=0
        phi_zero[:, [0,-1], :]=0
        phi_zero[:, :, [0,-1]]=0
        phi = phi_zero

        #define the first charge density
        charge = charge_density(lx, ly, lz) 
    
        sweeps = 0
        while True:
            previous_phi = np.copy(phi)
            for i in range(1, lx-1):
                for j in range(1, ly-1):
                    for k in range(1, lz-1):
                            phi[i,j,k] = (1/6)*(phi[i+1, j, k] + phi[i, j+1, k] + phi[i, j, k+1] + phi[i-1, j, k] + phi[i, j-1, k] + phi[i, j, k-1] + charge[i, j, k])

        #difference between new and old     
        #takes x sweeps #need to reset this
            
            dist = np.mean(np.abs(phi-previous_phi))
            print(dist)
            if np.isclose(dist,0,atol=5e-8): break  

            sweeps += 1
            print(" sweeps:" + str(sweeps), end="\r")


        data = [sweeps]
        with open(outfile, "a", newline='') as f:
            writer = csv.writer(f)
            writer.writerow(data)

        # reshaping the array from 3D
        # matrice to 2D matrice.
        arr_reshaped = phi.reshape(phi.shape[0], -1)
  
        # saving reshaped array to file.
        np.savetxt(f"Gauss_Seidel_Size{lx}_{condition}.txt", arr_reshaped)
        
        """
        plt.imshow(phi[0:size_added,size_small,0:size_added])
        plt.colorbar()
        plt.show()
        """

    elif condition == 'magnetic':

        #pad
        phi_zero[[0,-1], :, :]=0
        phi_zero[:, [0,-1], :]=0
        phi = phi_zero

        #define magnetic plane
        charge = magnetisation(lx, ly, lz)

        sweeps = 0

        while True:
            previous_phi = np.copy(phi)
            for i in range(1, lx-1):
                for j in range(1, ly-1):
                    for k in range(0, lz):
                            phi[i,j,k] = (1/6)*(phi[i+1, j, k] + phi[i, j+1, k] + phi[i, j, (k+1)%lz] + phi[i-1, j, k] + phi[i, j-1, k] + phi[i, j, (k-1)%lz] + charge[i, j, k])

        #difference between new and old     
        #takes x sweeps #need to reset this
            
            dist = np.mean(np.abs(phi-previous_phi))
            print(dist)
            if np.isclose(dist,0,atol=5e-8): break  

            sweeps += 1
            print(" sweeps:" + str(sweeps), end="\r")

        data = [sweeps]
        with open(outfile, "a", newline='') as f:
            writer = csv.writer(f)
            writer.writerow(data)

        # reshaping the array from 3D
        # matrice to 2D matrice.
        arr_reshaped = phi.reshape(phi.shape[0], -1)

        # saving reshaped array to file.
        np.savetxt(f"Gauss_Seidel_Size{lx}_{condition}.txt", arr_reshaped)
        
        """
        plt.imshow(phi[0:size_added,size_small,0:size_added])
        plt.colorbar()
        plt.show()
        """


def writefile_gaussseidel_sor(lx, condition,w):
    outfile = f"Gauss_Seidel_SOR_Size{lx}_{condition}_{w}.csv"
    header = ['Convergence_Sweeps']
    with open(outfile, "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)

def gauss_seidel_sor(lx, ly, lz, dt, dx, condition, w):
    
    outfile = f"Gauss_Seidel_SOR_Size{lx}_{condition}_{w}.csv"
    
    phi_zero = np.zeros((lx, ly, lz))

    # plot
    size_added = lx + 2
    size_small = lx // 2

    if condition == 'electric':

        #pad
        phi_zero[[0,-1], :, :]=0
        phi_zero[:, [0,-1], :]=0
        phi_zero[:, :, [0,-1]]=0
        phi = phi_zero

        #define the first charge density
        charge = charge_density(lx, ly, lz) 
    
        sweeps = 0
        while True:
            previous_phi = np.copy(phi)
            for i in range(1, lx-1):
                for j in range(1, ly-1):
                    for k in range(1, lz-1):
                            phi_gs = (1/6)*(phi[i+1, j, k] + phi[i, j+1, k] + phi[i, j, k+1] + phi[i-1, j, k] + phi[i, j-1, k] + phi[i, j, k-1] + charge[i, j, k])
                            phi[i,j,k] = w*phi_gs + (1-w)*previous_phi[i, j, k]

        #difference between new and old     
        #takes x sweeps #need to reset this
            
            dist = np.mean(np.abs(phi-previous_phi))
            print(dist)
            if np.isclose(dist,0,atol=5e-8): break  

            sweeps += 1
            print(" sweeps:" + str(sweeps), end="\r")


        data = [sweeps]
        with open(outfile, "a", newline='') as f:
            writer = csv.writer(f)
            writer.writerow(data)

        # reshaping the array from 3D
        # matrice to 2D matrice.
        arr_reshaped = phi.reshape(phi.shape[0], -1)
  
        # saving reshaped array to file.
        np.savetxt(f"Gauss_Seidel_SOR_Size{lx}_{condition}_{w}.txt", arr_reshaped)
        
        """
        plt.imshow(phi[0:size_added,size_small,0:size_added])
        plt.colorbar()
        plt.show()
        """

    elif condition == 'magnetic':

        #pad
        phi_zero[[0,-1], :, :]=0
        phi_zero[:, [0,-1], :]=0
        phi = phi_zero

        #define magnetic plane
        charge = magnetisation(lx, ly, lz)

        sweeps = 0

        while True:
            previous_phi = np.copy(phi)
            for i in range(1, lx-1):
                for j in range(1, ly-1):
                    for k in range(0, lz):
                            phi_gs = (1/6)*(phi[i+1, j, k] + phi[i, j+1, k] + phi[i, j, (k+1)%lz] + phi[i-1, j, k] + phi[i, j-1, k] + phi[i, j, (k-1)%lz] + charge[i, j, k])
                            phi[i,j,k] = w*phi_gs + (1-w)*previous_phi[i, j, k]
        
        #difference between new and old     
        #takes x sweeps #need to reset this
            
            dist = np.mean(np.abs(phi-previous_phi))
            print(dist)
            if np.isclose(dist,0,atol=5e-8): break  

            sweeps += 1
            print(" sweeps:" + str(sweeps), end="\r")

        data = [sweeps]
        with open(outfile, "a", newline='') as f:
            writer = csv.writer(f)
            writer.writerow(data)

        # reshaping the array from 3D
        # matrice to 2D matrice.
        arr_reshaped = phi.reshape(phi.shape[0], -1)

        # saving reshaped array to file.
        np.savetxt(f"Gauss_Seidel_SOR_Size{lx}_{condition}_{w}.txt", arr_reshaped)
        
        """
        plt.imshow(phi[0:size_added,size_small,0:size_added])
        plt.colorbar()
        plt.show()
        """


def main():

    dt = 1
    dx = 1

    initial = str(sys.argv[1])
    lx = int(sys.argv[2])
    ly = lx
    lz = lx


    if initial == 'cahn_hilliard_anim':
        phi = float(sys.argv[3])
        cahn_hilliard_anim(phi, lx, ly, 20000, dt, dx)
        

    elif initial == 'cahn_hilliard_file':
        phi = float(sys.argv[3])
        writefile_cahnhilliard(lx, phi)
        cahn_hilliard_file(phi, lx, ly, 1000000, dt, dx)
    
    elif initial == 'gauss_seidel_pad':
        condition = str(sys.argv[3])
        writefile_gaussseidel(lx, condition)
        gauss_seidel_pad(lx, ly, lz, dt, dx, condition)

    elif initial == 'jacobi_pad':
        condition = str(sys.argv[3])
        writefile_jacobi(lx, condition)
        jacobi_pad(lx, ly, lz, dt, dx, condition)

    elif initial == 'gauss_seidel_sor':
        
        condition = str(sys.argv[3])
        w = float(sys.argv[4])
        writefile_gaussseidel_sor(lx, condition, w) #write function
        gauss_seidel_sor(lx, ly, lz, dt, dx, condition, w)
        
main()