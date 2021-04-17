#!/usr/bin/env python3
import numpy as np
import hse_band_edge
import os
import scipy.constants as constants

def hse_polyfit_emass_vbm():
    x = np.zeros((20))
    y = np.zeros((20))
    
    data_band = hse_band_edge.band_calculation()
    
    #print(len(data_band.klsum_cart))
    x1 = data_band.klsum_cart[0:10]
    x2 = data_band.klsum_cart[10:20]
    
    y1 = data_band.eigen[0:10,data_band.nb_vbm]
    y2 = data_band.eigen[10:20,data_band.nb_vbm]
    
    direction1 = np.polyfit(x1, y1, 2, full=True)
    direction2 = np.polyfit(x2, y2, 2, full=True)
    
    d1_p2 = direction1[0][0]
    d1_p1 = direction1[0][1]
    d1_p0 = direction1[0][2]
    
    d2_p2 = direction2[0][0]
    d2_p1 = direction2[0][1]
    d2_p0 = direction2[0][2]
    
    y_fit1 = np.array([ d1_p2*t1*t1 + d1_p1*t1 + d1_p0 for t1 in x1])
    R_square1 = 1.0 - ((y_fit1 - y1)**2).sum() / ((y1 - y1.mean())**2).sum()
    y_fit2 = np.array([ d2_p2*t2*t2 + d2_p1*t2 + d2_p0 for t2 in x2])
    R_square2 = 1.0 - ((y_fit2 - y2)**2).sum() / ((y2 - y2.mean())**2).sum()
    
    residual1 = direction1[1]
    residual2 = direction2[1]
    
   #print('residual1: %f  residual2: %f' % (residual1, residual2))
   #print('R_square1: %f R_square2: %f' % (R_square1, R_square2))
   #print('directionx: %f %f %f' %(direction1[0][0], direction1[0][1], direction1[0][2]))
   #print('directiony: %f %f %f' %(direction2[0][0], direction2[0][1], direction2[0][2]))
    
    d1_p2_SI = d1_p2 * constants.electron_volt * constants.angstrom ** 2
    d2_p2_SI = d2_p2 * constants.electron_volt * constants.angstrom ** 2
    emass1 = constants.hbar ** 2 / 2.0 / d1_p2_SI
    emass2 = constants.hbar ** 2 / 2.0 / d2_p2_SI
    
    relative_emass1 = emass1 / constants.electron_mass
    relative_emass2 = emass2 / constants.electron_mass
    return (relative_emass1, relative_emass2, R_square1, R_square2)
    
if __name__ == '__main__':
    r_em1, r_em2, rs1, rs2 = hse_polyfit_emass_vbm()
    print('hole emass1: %16.8f r_square1: %16.8f' % (r_em1, rs1))
    print('hole emass2: %16.8f r_square2: %16.8f' % (r_em2, rs2))
