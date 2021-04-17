#!/usr/bin/env python3
import numpy as np
import numpy.linalg as LA
import os, sys
import matplotlib
import scipy.constants as constants
matplotlib.use('Agg')

# homemade script
import hse_band_edge
import hse_polyfit_emass_cbm
import hse_polyfit_emass_vbm
import get_c2d
import get_deformation_potential_constant

# homemade documents
# hse_polyfit_emass_cbm, unit: electron_mass
#    relative_emass1 = emass1 / constants.electron_mass
#    relative_emass2 = emass2 / constants.electron_mass
#    return (relative_emass1, relative_emass2, R_square1, R_square2)
#
# hse_polyfit_emass_vbm, unit: electron_mass
#    relative_emass1 = emass1 / constants.electron_mass
#    relative_emass2 = emass2 / constants.electron_mass
#    return (relative_emass1, relative_emass2, R_square1, R_square2)
#
# get_c2d, SI_unit: J/m2
#    return (c2d_SI[0], c2d_SI[1], r_square1, r_square2)
#
# get_deformation_potential_constant, eV
#    result_fit = [fit_cbm_x, fit_cbm_y, fit_vbm_x, fit_vbm_y]
#    result_r_square = [r_square_cbm_x, r_square_cbm_y, r_square_vbm_x, r_square_vbm_y]
#    return (result_fit, result_r_square)
#
# the directory architechture is at the end of this file


def calculate_mobility():
    ######################## initialization #######################
    print('initialization')
    # cbm_x, cbm_y, vbm_x, vbm_y_
    c2d = np.zeros((4), dtype=float)
    dp_constant = np.zeros((4), dtype=float)
    emass = np.zeros((4), dtype=float)
    emass_d = np.zeros((4), dtype=float)
    # r_squares
    rs_c2d = np.zeros((4), dtype=float)
    rs_dp = np.zeros((4), dtype=float)
    rs_emass = np.zeros((4), dtype=float)
    
    # record the root direcitory of our system
    system_dir = os.getcwd()
    cal_dir = [r'/c2d', r'/deformation_potential_constant', r'/emass/cbm', r'/emass/vbm']

    ######################## c2d, SI unit, J/m^2  #######################
    print('calculating c2d')
    os.chdir(system_dir + cal_dir[0])
    c2d1, c2d2, r_s1, r_s2 = get_c2d.get_c2d()
    c2d = np.array([c2d1, c2d2, c2d1, c2d2])
    # c2d_x, c2d_y, c2d_x, c2d_y
    rs_c2d = np.array([r_s1, r_s2, r_s1, r_s2])

    ####################### emass#####################################
    print('calculating emass')
    os.chdir(system_dir)
    os.chdir(system_dir + cal_dir[2])
    # cbm-x-emass, cbm-y-emass
    result_cbm = hse_polyfit_emass_cbm.hse_polyfit_emass_cbm()
    emass[0:2] = np.array(result_cbm[0:2])
    rs_emass[0:2] = np.array(result_cbm[2:4])
    os.chdir(system_dir)
    os.chdir(system_dir + cal_dir[3])
    # vbm-x-emass, vbm-y-emass
    result_vbm = hse_polyfit_emass_vbm.hse_polyfit_emass_vbm()
    emass[2:4] = np.array(result_vbm[0:2])
    rs_emass[2:4] = np.array(result_vbm[2:4])

    emass_d[0] = emass[1]
    emass_d[1] = emass[0]
    emass_d[2] = emass[3]
    emass_d[3] = emass[2]

    ######################## deformation potential constant, dp, unit: eV#################
    print('calculating deformation potential: ', end=' ', flush=True)
    os.chdir(system_dir)
    os.chdir(system_dir + cal_dir[1])
    # cbm-x, cbm-y, vbm-x, vbm-y
    result = get_deformation_potential_constant.get_deformation_potential_constant()
    dp_constant = np.array(result[0])
    rs_dp = result[1]
    
    ###################### 3 main variable end ##########################
    # formula
    # miu_2d = e * h_bar**3 * C2d / (kB * T * me * md * dp**2)
    Temp = 300 # room temperature

    up = 1.0 * constants.electron_volt
    up = up * constants.hbar**3
    up = up * c2d
    down = constants.Boltzmann * Temp
    down = down * emass
    down = down * np.sqrt(emass*emass_d)
    down = down * constants.electron_mass**2
    down = down * (dp_constant*constants.electron_volt)**2

    print('######################################################')
    print('mobility test:')
    print('  ', end='')
    print(up/down)
    print('######################################################')

    print('calculation details')
    print('c2d and R-square, unit: J*m^-2')
    print('  ', end='')
    print(c2d)
    print('  ', end='')
    print(rs_c2d)
    print('deformation_potential_constant and R-square, unit: eV')
    print('  ', end='')
    print(dp_constant)
    print('  ', end='')
    print(rs_dp)
    print('relative emass and R-square, unit: me_0')
    print('  ', end='')
    print(emass)
    print('  ', end='')
    print(rs_emass)

    hbar = constants.hbar
    ev = constants.electron_volt
    em = constants.electron_mass
    kb = constants.Boltzmann
    # mobility = constants.electron_volt * constants.hbar**3 * c2d / (constants.Boltzmann * Temp * emass * np.sqrt(emass*emass_d) * constants.electron_mass**2 * (dp_constant*constants.electron_volt)**2)
    # unit = hbar ** 3 / (kb * e**2 * m**2)
    unit = hbar / kb * hbar / ev * hbar / em / em
    print('Temperature: %d' % Temp)
    mobility_float = c2d / (Temp * emass * np.sqrt(emass*emass_d) * dp_constant**2) 
    mobility = mobility_float * unit

    print('mobility in SI unit: m^2*V^-1*s^-1')
    print('cbm_x, cbm_y, vbm_x, vbm_y:')
    print('  ', end='')
    print(mobility)

    print('mobility in conventional unit: 10^3 cm^2*V^-1*s^-1')
    print('cbm_x, cbm_y, vbm_x, vbm_y:')
    print('  ', end='')
    print(mobility*10.0)
    #test = constants.electron_volt * constants.hbar**3 * c2d / (constants.Boltzmann * Temp * emass * np.sqrt(emass*emass_d) * constants.electron_mass**2 * (dp_constant*constants.electron_volt)**2)
    #print('testmobility: ')
    #print(test)
    

if __name__ == '__main__':
    result = calculate_mobility()


# calculate_mobility
#    directory architechture
#       
#           |         |-- 0
#           |         |-- 0_0.990
#           |         |-- 0_0.995
#           |         |-- 0_1.000
#           |         |-- 0_1.005
#           |         |-- 0_1.010
#           |-- c2d --|-- 1_0.990
#           |         |-- 1_0.995
#           |         |-- 1_1.000
#           |         |-- 1_1.005
#           |         `-- 1_1.010
#           |
#           |                                    |-- 0
#           |                                    |-- 0_0.990
#           |                                    |-- 0_0.995
#           |                                    |-- 0_1.000
#           |                                    |-- 0_1.005
#           |-- deformation_potential_constant --|-- 0_1.010 
#           |                                    |-- 1_0.990
#           |                                    |-- 1_0.995
#           |                                    |-- 1_1.000
#           |                                    |-- 1_1.005
#           |                                    `-- 1_1.010
#           |
#           |
#           |           |-- cbm
#           |           |      
#  system --|-- emass --|      
#           |           |      
#           |           `-- vbm
#           |
#           |
#           |-- relax
#           |-- scf
#           `-- strain
#    
#
#
#
