#!/usr/bin/python

from sys import argv;
from MDAnalysis import Universe;
from glob import glob;
import numpy as np
import sys;
import os;

ifile = file('conserHolo_check.txt','r')
pro_marks = {};
for iline in ifile:
    iwords = iline.split();
    pro_marks[ iwords[0] ] = iwords[1:]

pro_names = pro_marks.keys()

#pro_ = os.getcwd().split('/')[-1]
pro_name_ = sys.argv[1];
pro_name = pro_name_.split('.')[0];

prmfile = sys.argv[1];
trjfile = sys.argv[2];
univ = Universe( prmfile, trjfile );

# define the dihedral
i_res = int(pro_marks[pro_name][4]) + 1
pdY_phi = univ.selectAtoms( 'resnum ' + str(i_res-1) + ' and name C',
                            'resnum ' + str(i_res) + ' and name N',
                            'resnum ' + str(i_res) + ' and name CA',
                            'resnum ' + str(i_res) + ' and name C');  #C-1-N-CA-C
pdY_psi = univ.selectAtoms( 'resnum ' + str(i_res) + ' and name N',
                            'resnum ' + str(i_res) + ' and name CA',
                            'resnum ' + str(i_res) + ' and name C',
                            'resnum ' + str(i_res+1) + ' and name N'); #N-CA-C-N+1
pdY_O = univ.selectAtoms( 'resnum ' + str(i_res) + ' and name CA',
                             'resnum ' + str(i_res) + ' and name CB',
                             'resnum ' + str(i_res) + ' and name OH',
                             'resnum ' + str(i_res) + ' and name HH' ); #CA-CB-OH-HH
#pdY_chi2 = univ.selectAtoms( 'resnum ' + str(i_res) + ' and name CA',
#                             'resnum ' + str(i_res) + ' and name CB',
#                             'resnum ' + str(i_res) + ' and name CG',
#                             'resnum ' + str(i_res) + ' and name CD1'); #CA-CB-CG-CD1
print i_res, pdY_psi[0].resname;

j_res = int(pro_marks[pro_name][9])
N1_ND2 = univ.selectAtoms( 'resnum ' + str(j_res) + ' and name CA',
                            'resnum ' + str(j_res) + ' and name CB',
                            'resnum ' + str(j_res) + ' and name CG',
                            'resnum ' + str(j_res) + ' and name OD1') # CA-CB-CG-OD1
#print j_res, N1_ND2[0].resname,N1_ND2[1].resname;

k_res = int(pro_marks[pro_name][10])
#Yn_chi = univ.selectAtoms( 'resnum ' + str(k_res) + ' and name CA',
#                            'resnum ' + str(k_res) + ' and name CB',
#                            'resnum ' + str(k_res) + ' and name OH',
#                            'resnum ' + str(k_res) + ' and name HH' ); # CA-CB-OH-HH
Yn_chi = univ.selectAtoms( 'resnum ' + str(k_res) + ' and name CA',
                            'resnum ' + str(k_res) + ' and name CB',
                            'resnum ' + str(k_res) + ' and name CG',
                            'resnum ' + str(k_res) + ' and name CD1' ); # CA-CB-CZ-HZ for BAZ2B, SMARCA2/4,TIF1A Fn instead of Yn
#print k_res, Yn_chi1[0].resname;

l_res = int(pro_marks[pro_name][0])
wPf_phi = univ.selectAtoms( 'resnum ' + str(l_res-1) + ' and name C',
                            'resnum ' + str(l_res) + ' and name N',
                            'resnum ' + str(l_res) + ' and name CA',
                            'resnum ' + str(l_res) + ' and name C') # C-1-N-NA-C
wPf_psi = univ.selectAtoms( 'resnum ' + str(l_res) + ' and name N',
                            'resnum ' + str(l_res) + ' and name CA',
                            'resnum ' + str(l_res) + ' and name C',
                            'resnum ' + str(l_res+1) + ' and name N') # N-CA-C-N+1
print l_res, wPf_phi[0].resname;

m_res = int(pro_marks[pro_name][5])
pMd_phi = univ.selectAtoms( 'resnum ' + str(m_res-1) + ' and name C',
                            'resnum ' + str(m_res) + ' and name N',
                            'resnum ' + str(m_res) + ' and name CA',
                            'resnum ' + str(m_res) + ' and name C') # C-1-N-NA-C
pMd_psi = univ.selectAtoms( 'resnum ' + str(m_res) + ' and name N',
                            'resnum ' + str(m_res) + ' and name CA',
                            'resnum ' + str(m_res) + ' and name C',
                            'resnum ' + str(m_res+1) + ' and name N'); #N-CA-C-N+1
n_res = int(pro_marks[pro_name][1])-2
preP_phi = univ.selectAtoms( 'resnum ' + str(n_res-1) + ' and name C',
                            'resnum ' + str(n_res) + ' and name N',
                            'resnum ' + str(n_res) + ' and name CA',
                            'resnum ' + str(n_res) + ' and name C') # C-1-N-NA-C
preP_psi = univ.selectAtoms( 'resnum ' + str(n_res) + ' and name N',
                            'resnum ' + str(n_res) + ' and name CA',
                            'resnum ' + str(n_res) + ' and name C',
                            'resnum ' + str(n_res+1) + ' and name N') # N-CA-C-N+1

ofile = file( pro_name_ + '_waterHB_Dih.dat', 'w');
print >>ofile, 'frame'.ljust(10) + 'pdY_phi pdY_psi pdY_O N1_ND2 Yn_chi wPf_phi wPf_psi pMd_phi pMd_psi preP_phi preP_psi';
for ts in univ.trajectory:
    print >>ofile, str(ts.frame).ljust(10) + ' ' + str(pdY_phi.dihedral()) + ' ' + \
          str(pdY_psi.dihedral()) + ' ' + str(pdY_O.dihedral()) + ' ' + str(N1_ND2.dihedral()) + \
          ' ' + str(Yn_chi.dihedral()) + ' ' + str(wPf_phi.dihedral()) + ' ' + str(wPf_psi.dihedral())+ \
	' '+str(pMd_phi.dihedral()) + ' ' +str(pMd_psi.dihedral()) + ' ' + str(preP_phi.dihedral()) + ' ' +str(preP_psi.dihedral());
ofile.close();

