#!/usr/bin/env python3
import numpy as np
import numpy.linalg as LA
import os, sys
import pymatgen.io.vasp.outputs as vaspout

import hse_band_edge


def get_vacuum(filename='LOCPOT'):
    # static_eletronic_potential
    sep = vaspout.Locpot.from_file(filename)
    data_sep = sep.data['total']
    pot_mesh = data_sep.shape
    result_z = np.einsum('abc->c', data_sep)
    return result_z.max() / pot_mesh[0] / pot_mesh[1]


def get_fermi_cbm_minus_fermi_vbm_minus_fermi():
    data_band = hse_band_edge.band_calculation()
    e_fermi = data_band.efermi
    e_cbm = data_band.e_cbm
    e_vbm = data_band.e_vbm
    return e_fermi, e_cbm, e_vbm


def polyfit(xx, yy, rank=2):
    result = np.polyfit(xx, yy, rank, full=True)
    num_fit = len(xx)
    fit = np.zeros((num_fit))
    for i in range(num_fit):
        for j in range(rank+1):
            fit[i] += result[0][j] * xx[i] ** (rank-j)

    R_square = 1.0 - ((fit - yy)**2).sum() / ((yy - yy.mean())**2).sum()
    return (result, R_square)


def get_deformation_potential_constant():
    i_direction = [0, 1]
    s_direction = ['0', '1']
    f_strain = np.array([0.990, 0.995, 1.000, 1.005, 1.010])
    s_strain = ['0.990', '0.995', '1.000', '1.005', '1.010']

    fermi = np.zeros((2,5))
    cbm_minus_fermi = np.zeros((2,5))
    vbm_minus_fermi = np.zeros((2,5))
    vacuum = np.zeros((2,5))

    for i in range(2):
        for j in range(5):
            print(s_direction[i]+'_'+s_strain[j], end='  ', flush=True)
            currentdir = os.getcwd()
            os.chdir(s_direction[i]+'_'+s_strain[j])
            vacuum[i][j] = get_vacuum()
            t1, t2, t3 =  get_fermi_cbm_minus_fermi_vbm_minus_fermi()
            fermi[i][j] = t1
            cbm_minus_fermi[i][j] = t2
            vbm_minus_fermi[i][j] = t3
            os.chdir(currentdir)
    print('')
   #print('data')
   #print(vacuum[0])
   #print(fermi[0])
   #print(cbm_minus_fermi[0])
   #print(vbm_minus_fermi[0])
    cbm_minus_vacuum = cbm_minus_fermi + fermi - vacuum
    vbm_minus_vacuum = vbm_minus_fermi + fermi - vacuum
    fit_cbm_x, r_square_cbm_x = polyfit(f_strain, cbm_minus_vacuum[0], rank=1)
    fit_cbm_y, r_square_cbm_y = polyfit(f_strain, cbm_minus_vacuum[1], rank=1)
    fit_vbm_x, r_square_vbm_x = polyfit(f_strain, vbm_minus_vacuum[0], rank=1)
    fit_vbm_y, r_square_vbm_y = polyfit(f_strain, vbm_minus_vacuum[1], rank=1)
   #print('fit_result')
   #print(fit_cbm_x[0][0], r_square_cbm_x)
   #print(fit_cbm_y[0][0], r_square_cbm_y)
   #print(fit_vbm_x[0][0], r_square_vbm_x)
   #print(fit_vbm_y[0][0], r_square_vbm_y)
    result_fit = [fit_cbm_x[0][0], fit_cbm_y[0][0], fit_vbm_x[0][0], fit_vbm_y[0][0]]
    result_r_square = [r_square_cbm_x, r_square_cbm_y, r_square_vbm_x, r_square_vbm_y]
    return (result_fit, result_r_square)

    
if __name__ == '__main__':
    result_fit, result_r = get_deformation_potential_constant()
    for i in range(4):
        print(result_fit[i], result_r[i])

