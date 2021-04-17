#!/usr/bin/env python3
import numpy as np
import numpy.linalg as LA
import os, sys
import pymatgen.io.vasp.outputs as vaspout
import scipy.constants as constants

import hse_band_edge

def final_energy(filename='OSZICAR'):
    # static_eletronic_potential
    oszicar = vaspout.Oszicar(filename)
    final_energy = oszicar.final_energy
    return final_energy


def polyfit(xx, yy, rank=2):
    result = np.polyfit(xx, yy, rank, full=True)
    num_fit = len(xx)
    fit = np.zeros((num_fit))
    for i in range(num_fit):
        for j in range(rank+1):
            fit[i] += result[0][j] * xx[i] ** (rank-j)

    R_square = 1.0 - ((fit - yy)**2).sum() / ((yy - yy.mean())**2).sum()
    return (result, R_square)


def get_lattice(dir_name='.'):
    current_dir = os.getcwd()
    os.chdir(dir_name)
    data_band = hse_band_edge.band_calculation()
    os.chdir(current_dir)

    return data_band.lat
    #return vasprun.final_structure.lattice._matrix
    
    
def get_c2d():
    i_direction = [0, 1]
    s_direction = ['0', '1']
    f_strain = [0.990, 0.995, 1.000, 1.005, 1.010]
    s_strain = ['0.990', '0.995', '1.000', '1.005', '1.010']
    total_energy = np.zeros((2,5), dtype=float)

    for i in range(2):
        for j in range(5):
            dir_name = s_direction[i] + '_' + s_strain[j]
            total_energy[i][j] = final_energy(filename=dir_name+r'/OSZICAR')
    
    lat = get_lattice(dir_name='0_1.000')
    area = LA.norm(np.cross(lat[0], lat[1]))

    direction1, r_square1 = polyfit(f_strain, total_energy[0])
    direction2, r_square2 = polyfit(f_strain, total_energy[1])

    d1_p2 = direction1[0][0]
    d1_p1 = direction1[0][1]
    d1_p0 = direction1[0][2]

    d2_p2 = direction2[0][0]
    d2_p1 = direction2[0][1]
    d2_p0 = direction2[0][2]

    c2d = np.zeros((2))
    c2d[0] = d1_p2 * 2.0 / area
    c2d[1] = d2_p2 * 2.0 / area
    
    c2d_SI = c2d * constants.electron_volt / (constants.angstrom**2)
    return (c2d_SI[0], c2d_SI[1], r_square1, r_square2)

   #print('d1_p2, r_square1: %12.6f    %12.6f' %(d1_p2, r_square1))
   #print('d1_p2, r_square1: %12.6f    %12.6f' %(d1_p2, r_square1))

if __name__ == '__main__':
    result = get_c2d()
    c2d_SI_x = result[0]
    c2d_SI_y = result[1]
    r_square_x = result[2]
    r_square_y = result[3]
    print('c2d_x: %12.6f    %12.6f' % (c2d_SI_x, r_square_x))
    print('c2d_y: %12.6f    %12.6f' % (c2d_SI_y, r_square_y))
