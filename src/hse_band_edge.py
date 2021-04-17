#!/usr/bin/env python3
# reclat here is reclat_outcar * 2pi, with unit: A^-1
# k_cart = self.k_frac x reclat
#                        [rec00 rec01 rec02]
# [c0 c1 c2] = [f0 f1 f2][rec10 rec11 rec12]
#                        [rec20 rec21 rec22]
import numpy as np
import numpy.linalg as LA
import matplotlib
matplotlib.use('Agg')
import os, sys


class band_calculation():
    def __init__(self):
        self.lat = np.zeros((3,3))
        self.reclat = np.zeros((3,3)) # reclat = np.transpose(np.inv(lat)), this is in vasp OUTCAR
        self.k_cbm = np.zeros((3))
        self.k_vbm = np.zeros((3))
        self.e_cbm = 0.0 
        self.e_vbm = 0.0

        self.read_band_cal()
        self.cal_kl()
        self.write_band()
        self.find_band_edge()
        #self.output()


    def read_band_cal(self):
        os.system('grep E-fermi OUTCAR | tail -1 > tmp')
        os.system('grep -A3 "reciprocal lattice vectors" OUTCAR | sed -n 2,4p >> tmp')
        
        ftmp = open('tmp', 'r')
        txt = ftmp.readlines()
        self.efermi = float(txt[0].split()[2])
        for i in range(3):
            txt_tmp = [float(i) for i in txt[i+1].split()]
            self.lat[i][0:3] = np.array(txt_tmp[0:3])
            self.reclat[i][0:3] = np.pi * 2 * np.array(txt_tmp[3:6])
        
        ftmp.close()
        f_kpt = open('KPOINTS', 'r')
        tmp_kpt = f_kpt.readlines()
        f_kpt.close()

        self.num_tot_hsekpt = int(tmp_kpt[1].split()[0])
        self.num_scf_kpt = 0
        for i in range(3,self.num_tot_hsekpt+3):
            if float(tmp_kpt[i].split()[3]) > 0.0001:
                self.num_scf_kpt += 1


        self.num_band_kpt = self.num_tot_hsekpt - self.num_scf_kpt
        

        f_eigen = open('EIGENVAL', 'r')
        for i in range(6):
            tmp = f_eigen.readline()
        
        # self.nk = int(tmp.split()[1])
        self.nk = self.num_band_kpt
        self.nband = int(tmp.split()[2])
        
        self.k_frac = np.zeros((self.nk,3))
        self.eigen = np.zeros((self.nk,self.nband))
        self.kl_frac = np.zeros((self.nk))
        self.kl_cart = np.zeros((self.nk))
        self.klsum_frac = np.zeros((self.nk))
        self.klsum_cart = np.zeros((self.nk))

        self.k_scf = np.zeros((self.num_scf_kpt,3), dtype=float)
        self.k_scf_weight = np.zeros((self.num_scf_kpt), dtype=float)

        # read scf kpts whose weights ne 0.0
        for i in range(self.num_scf_kpt):
            tmp = f_eigen.readline()
            tmp = f_eigen.readline()
            tk = tmp.split()
            self.k_scf[i][0] = float(tk[0])
            self.k_scf[i][1] = float(tk[1])
            self.k_scf[i][2] = float(tk[2])
            self.k_scf_weight[i] = float(tk[3])
            for j in range(self.nband):
                tmp_band = f_eigen.readline()
            
        
        # read k_frac and eigen
        for i in range(self.num_band_kpt):
            tmp = f_eigen.readline()
            tmp = f_eigen.readline()
            tk = tmp.split()
            self.k_frac[i][0] = float(tk[0])
            self.k_frac[i][1] = float(tk[1])
            self.k_frac[i][2] = float(tk[2])
            for j in range(self.nband):
                tmp_band = f_eigen.readline()
                self.eigen[i][j] = float(tmp_band.split()[1])
        
        f_eigen.close()


    def show_eigen(self, ki, bi):
        print('eigenvalue at (kpt_number, band_number). eigen_shifted(%d,%d)=%12.6f' % (ki, bi, self.eigen_shifted[ki, bi]))
        

    # calcuself.late self.kl_cart and klsum
    # self.kl_cart[0] = klsum[0] = 0
    def cal_kl(self):
        for i in range(1,self.nk):
            k_frac_delta = self.k_frac[i] - self.k_frac[i-1]
            k_cart_delta = np.dot(k_frac_delta, self.reclat)
        
            self.kl_frac[i] = LA.norm(k_frac_delta)
            self.kl_cart[i] = LA.norm(k_cart_delta)
            self.klsum_frac[i] = 0.0;
            self.klsum_cart[i] = 0.0;
            for j in range(i+1):
                self.klsum_frac[i] += self.kl_frac[j]
                self.klsum_cart[i] += self.kl_cart[j]

    # write band.dat
    def write_band(self):
        fout = open('band_frac.dat', 'w')
        fout.write('self.kl_frac\tenergy\n')
        for i in range(self.nband):
            for j in range(self.nk):
                fout.write('%8.4f\t%8.4f\n' % (self.klsum_frac[j], self.eigen[j][i] - self.efermi))
        
            fout.write('\n')
        fout.close()
        
        fout = open('band_cart.dat', 'w')
        fout.write('self.kl_cart\tenergy\n')
        for i in range(self.nband):
            for j in range(self.nk):
                fout.write('%8.4f\t%8.4f\n' % (self.klsum_cart[j], self.eigen[j][i] - self.efermi))
        
            fout.write('\n')
        fout.close()
    # find band edge
    def find_band_edge(self):
        self.e_cbm = 10.0
        self.e_vbm = -10.0
        self.k_cbm = np.zeros((3))
        self.k_vbm = np.zeros((3))
        self.klc = 0.0
        self.klv = 0.0
        self.klcf = 0.0
        self.klvf = 0.0
        self.nk_cbm = 0
        self.nk_vbm = 0
        self.nb_cbm = 0
        self.nb_vbm = 0
        
        self.eigen_shifted = self.eigen - self.efermi
        for i in range(self.nk):
            for j in range(self.nband):
                if self.eigen_shifted[i][j] > 0 and self.eigen_shifted[i][j] < self.e_cbm:
                    self.e_cbm = self.eigen_shifted[i][j]
                    self.nk_cbm = i
                    self.nb_cbm = j
                    self.k_cbm = self.k_frac[i]
                    self.klc = self.klsum_cart[i]
                    self.klcf = self.klsum_frac[i]
                    
                    
                if self.eigen_shifted[i][j] < 0 and self.eigen_shifted[i][j] > self.e_vbm:
                    self.e_vbm = self.eigen_shifted[i][j]
                    self.nk_vbm = i
                    self.nb_vbm = j
                    self.k_vbm = self.k_frac[i]
                    self.klv = self.klsum_cart[i]
                    self.klvf = self.klsum_frac[i]



    def output(self):
        if self.klc == self.klv:
            print('direct gap')
        else:
            print('indirect gap')

        print('Efermi = %16.10f' % self.efermi)
        print('(self.kl_cart, CBM-Fermi)= %8.4f %8.4f' % (self.klc, self.e_cbm))
        print('(self.kl_cart, VBM-Fermi)= %8.4f %8.4f' % (self.klv, self.e_vbm))
        print('(self.kl_frac, CBM-Fermi)= %8.4f %8.4f' % (self.klcf, self.e_cbm))
        print('(self.kl_frac, VBM-Fermi)= %8.4f %8.4f' % (self.klvf, self.e_vbm))
        print('gap = %8.4f' % (self.e_cbm-self.e_vbm))
        print('kc_frac= %8.4f %8.4f %8.4f' % (self.k_cbm[0], self.k_cbm[1], self.k_cbm[2]))
        print('kv_frac= %8.4f %8.4f %8.4f' % (self.k_vbm[0], self.k_vbm[1], self.k_vbm[2]))
        print('band_egde. nk_cbm, nk_vbm, nb_cbm, nb_vbm: %d, %d, %d, %d' % (self.nk_cbm, self.nk_vbm, self.nb_cbm, self.nb_vbm))

    def write_effect_mass_kpt(self):
        xx = np.array([0.1, 0.0, 0.0])
        yy = np.array([0.0, 0.1, 0.0])
        t2 = 2.0/3.0*self.reclat[1]
        t1 = 1.0/3.0*self.reclat[0]

        xc_0 = np.zeros((3))
        xc_1 = np.zeros((3))
        yc_0 = np.zeros((3))
        yc_1 = np.zeros((3))
        xv_0 = np.zeros((3))
        xv_1 = np.zeros((3))
        yv_0 = np.zeros((3))
        yv_1 = np.zeros((3))
        
        xc_0 = self.k_cbm - xx
        xc_1 = self.k_cbm + xx
        xv_0 = self.k_vbm - xx
        xv_1 = self.k_vbm + xx
       #
       #fsym = open('sym', 'r')
       #txt = fsym.readline()
       #sym = txt.split()[0]

       #if sym == 'square' or sym == 'rectangle':
       #    yc_0 = self.k_cbm - yy
       #    yc_1 = self.k_cbm + yy
       #    yv_0 = self.k_vbm - yy
       #    yv_1 = self.k_vbm + yy
       #elif sym == 'hex60':  # reciprocal 120 degree
       #    yc_0 = self.k_cbm - 0.1 * (t2 + t1)
       #    yc_1 = self.k_cbm + 0.1 * (t2 + t1)
       #    yv_0 = self.k_vbm - 0.1 * (t2 + t1)
       #    yv_1 = self.k_vbm + 0.1 * (t2 + t1)
       #elif sym == 'hex120':  # reciprocal 60 degree
       #    yc_0 = self.k_cbm - 0.1 * (t2 - t1)
       #    yc_1 = self.k_cbm + 0.1 * (t2 - t1)
       #    yv_0 = self.k_vbm - 0.1 * (t2 - t1)
       #    yv_1 = self.k_vbm + 0.1 * (t2 - t1)
       #else:
       #    sys.exit()
       #
        try:
            fsym = open('sym', 'r')
            txt = fsym.readline()
            sym = txt.split()[0]
         
            if sym == 'square' or sym == 'rectangle':
                yc_0 = self.k_cbm - yy
                yc_1 = self.k_cbm + yy
                yv_0 = self.k_vbm - yy
                yv_1 = self.k_vbm + yy
            elif sym == 'hex60':  # reciprocal 120 degree
                yc_0 = self.k_cbm - 0.1 * (t2 + t1)
                yc_1 = self.k_cbm + 0.1 * (t2 + t1)
                yv_0 = self.k_vbm - 0.1 * (t2 + t1)
                yv_1 = self.k_vbm + 0.1 * (t2 + t1)
            elif sym == 'hex120':  # reciprocal 60 degree
                yc_0 = self.k_cbm - 0.1 * (t2 - t1)
                yc_1 = self.k_cbm + 0.1 * (t2 - t1)
                yv_0 = self.k_vbm - 0.1 * (t2 - t1)
                yv_1 = self.k_vbm + 0.1 * (t2 - t1)
            else:
                print('unsupported sym')
                sys.exit()
        except IOError:
            print('file "sym" not found, use hex120')
            yc_0 = self.k_cbm - 0.1 * (t2 - t1)
            yc_1 = self.k_cbm + 0.1 * (t2 - t1)
            yv_0 = self.k_vbm - 0.1 * (t2 - t1)
            yv_1 = self.k_vbm + 0.1 * (t2 - t1)


        xc_array = np.zeros((10,3), dtype=float)
        yc_array = np.zeros((10,3), dtype=float)
        xv_array = np.zeros((10,3), dtype=float)
        yv_array = np.zeros((10,3), dtype=float)

        for i in range(10):
            xc_array[i] = xc_0 * float(10-i) / 10.0 + xc_1 * float(i) / 10.0
            yc_array[i] = yc_0 * float(10-i) / 10.0 + yc_1 * float(i) / 10.0
            xv_array[i] = xv_0 * float(10-i) / 10.0 + xv_1 * float(i) / 10.0
            yv_array[i] = yv_0 * float(10-i) / 10.0 + yv_1 * float(i) / 10.0


        fhsecbm = open('kpt_hse_cbm', 'w')
        fhsecbm.write('k-points along high symmetry lines\n')
        fhsecbm.write('%d\n' % (20+self.num_scf_kpt))
        fhsecbm.write('Reciprocal\n')
        for i in range(self.num_scf_kpt):
            fhsecbm.write('%20.14f %20.14f %20.14f    %20.14f\n' % (self.k_scf[i][0], self.k_scf[i][1], self.k_scf[i][2], self.k_scf_weight[i]))

        for i in range(10):
            fhsecbm.write('%20.14f %20.14f %20.14f    0.0\n' % (xc_array[i][0], xc_array[i][1], xc_array[i][2])) 

        for i in range(10):
            fhsecbm.write('%20.14f %20.14f %20.14f    0.0\n' % (yc_array[i][0], yc_array[i][1], yc_array[i][2])) 
        fhsecbm.close()


        fhsevbm = open('kpt_hse_vbm', 'w')
        fhsevbm.write('k-points along high symmetry lines\n')
        fhsevbm.write('%d\n' % (20+self.num_scf_kpt))
        fhsevbm.write('Reciprocal\n')
        for i in range(self.num_scf_kpt):
            fhsevbm.write('%20.14f %20.14f %20.14f    %20.14f\n' % (self.k_scf[i][0], self.k_scf[i][1], self.k_scf[i][2], self.k_scf_weight[i]))

        for i in range(10):
            fhsevbm.write('%20.14f %20.14f %20.14f    0.0\n' % (xv_array[i][0], xv_array[i][1], xv_array[i][2])) 

        for i in range(10):
            fhsevbm.write('%20.14f %20.14f %20.14f    0.0\n' % (yv_array[i][0], yv_array[i][1], yv_array[i][2])) 
        fhsevbm.close()


        fcbm = open('kpt_cbm', 'w')
        fcbm.write('kpts for effect_mass -0.1 ~ 0.1 of cbm along x or zigzag direction\n')
        fcbm.write('20\n')
        fcbm.write('line mode\n')
        fcbm.write('reciprocal coordinates\n')
        fcbm.write('%8.5f %8.5f %8.5f    0\n' % (xc_0[0], xc_0[1], xc_0[2]))
        fcbm.write('%8.5f %8.5f %8.5f    0\n' % (xc_1[0], xc_1[1], xc_1[2]))
        fcbm.write('                      \n')
        fcbm.write('%8.5f %8.5f %8.5f    0\n' % (yc_0[0], yc_0[1], yc_0[2]))
        fcbm.write('%8.5f %8.5f %8.5f    0\n' % (yc_1[0], yc_1[1], yc_1[2]))
        fcbm.close()
        
        fvbm = open('kpt_vbm', 'w')
        fvbm.write('kpts for effect_mass -0.1 ~ 0.1 of vbm along y or armchair direction\n')
        fvbm.write('20\n')
        fvbm.write('line mode\n')
        fvbm.write('reciprocal coordinates\n')
        fvbm.write('%8.5f %8.5f %8.5f    0\n' % (xv_0[0], xv_0[1], xv_0[2]))
        fvbm.write('%8.5f %8.5f %8.5f    0\n' % (xv_1[0], xv_1[1], xv_1[2]))
        fvbm.write('                      \n')
        fvbm.write('%8.5f %8.5f %8.5f    0\n' % (yv_0[0], yv_0[1], yv_0[2]))
        fvbm.write('%8.5f %8.5f %8.5f    0\n' % (yv_1[0], yv_1[1], yv_1[2]))
        fvbm.close()
        
        

if __name__ == '__main__':
    b1 = band_calculation()
    b1.write_effect_mass_kpt()
    if len(sys.argv) > 2:
        ki = int(sys.argv[1])
        bi = int(sys.argv[2])
        b1.show_eigen(ki, bi)
    b1.output()
#b1. __init__()
#b1.cal_kl()
#b1.write_band()
#b1.find_band_edge()
#b1.output()
