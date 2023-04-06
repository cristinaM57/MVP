import os
import numpy as np


def main():

    lx = 50

    w_range = np.arange(1.7, 1.9, 0.01)

    electric = "electric"

    gauss_seidel_sor = "gauss_seidel_sor"
    
    print('Batch Simulation Starting \n')

    for i in range(len(w_range)):
        cmd = f'python checkpoint3.py gauss_seidel_sor 50 electric {np.round(w_range[i], 2)}'
        os.system(cmd)

        print(f'Succesful Simulation @ Parameters: w={np.round(w_range[i], 2)}')

    print('\nBatch Simulations Complete')


main()