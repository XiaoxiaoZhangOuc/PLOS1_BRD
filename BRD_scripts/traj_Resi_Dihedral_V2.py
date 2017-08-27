#!/usr/bin/python

from sys import argv;
from MDAnalysis import Universe;
from glob import glob;
import numpy as np
import sys;
import os;
sys.path.append('/share/udata1/xiaoxiao/BRD_Ana_15.4.22/scripts/');
from pro_conser_marks import pro_conser_marks;

pro_name_ = os.getcwd().split('/')[-1];
pro_name = pro_name_.split('_run')[0];
print 'pro_name', pro_name;

prmfile = sys.argv[1];
trjfile = sys.argv[2];
univ = Universe( prmfile, trjfile );

# define the dihedral
i_res = int(pro_conser_marks(pro_name)['pDy']) + 1
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
pdY_chi = univ.selectAtoms( 'resnum ' + str(i_res) + ' and name N',
                             'resnum ' + str(i_res) + ' and name CA',
                             'resnum ' + str(i_res) + ' and name CB',
                             'resnum ' + str(i_res) + ' and name CG'); #N-CA-CB-CG
print i_res, pdY_psi[0].resname;

j_res = int(pro_conser_marks(pro_name)['N1'])
N1_ND2 = univ.selectAtoms( 'resnum ' + str(j_res) + ' and name CA',
                            'resnum ' + str(j_res) + ' and name CB',
                            'resnum ' + str(j_res) + ' and name CG',
                            'resnum ' + str(j_res) + ' and name ND2') # CA-CB-CG-ND2
print j_res, N1_ND2[0].resname,N1_ND2[1].resname;

k_res = int(pro_conser_marks(pro_name)['Yn'])
Yn_O = univ.selectAtoms( 'resnum ' + str(k_res) + ' and name CA',
                            'resnum ' + str(k_res) + ' and name CB',
                            'resnum ' + str(k_res) + ' and name OH',
                            'resnum ' + str(k_res) + ' and name HH' ); # CA-CB-OH-HH
Yn_chi = univ.selectAtoms( 'resnum ' + str(k_res) + ' and name N',
                            'resnum ' + str(k_res) + ' and name CA',
                            'resnum ' + str(k_res) + ' and name CB',
                            'resnum ' + str(k_res) + ' and name CG' ); # CA-CB-CZ-HZ for BAZ2B, SMARCA2/4,TIF1A Fn instead of Yn
print k_res, Yn_chi[0].resname;

l_res = int(pro_conser_marks(pro_name)['wPf'])
wPf_phi = univ.selectAtoms( 'resnum ' + str(l_res-1) + ' and name C',
                            'resnum ' + str(l_res) + ' and name N',
                            'resnum ' + str(l_res) + ' and name CA',
                            'resnum ' + str(l_res) + ' and name C') # C-1-N-NA-C
wPf_psi = univ.selectAtoms( 'resnum ' + str(l_res) + ' and name N',
                            'resnum ' + str(l_res) + ' and name CA',
                            'resnum ' + str(l_res) + ' and name C',
                            'resnum ' + str(l_res+1) + ' and name N') # N-CA-C-N+1
print l_res, wPf_phi[0].resname;

m_res = int(pro_conser_marks(pro_name)['pMd'])
pMd_phi = univ.selectAtoms( 'resnum ' + str(m_res-1) + ' and name C',
                            'resnum ' + str(m_res) + ' and name N',
                            'resnum ' + str(m_res) + ' and name CA',
                            'resnum ' + str(m_res) + ' and name C') # C-1-N-NA-C
pMd_psi = univ.selectAtoms( 'resnum ' + str(m_res) + ' and name N',
                            'resnum ' + str(m_res) + ' and name CA',
                            'resnum ' + str(m_res) + ' and name C',
                            'resnum ' + str(m_res+1) + ' and name N'); #N-CA-C-N+1
n_res = int(pro_conser_marks(pro_name)['pVd'])-2
preP_phi = univ.selectAtoms( 'resnum ' + str(n_res-1) + ' and name C',
                            'resnum ' + str(n_res) + ' and name N',
                            'resnum ' + str(n_res) + ' and name CA',
                            'resnum ' + str(n_res) + ' and name C') # C-1-N-NA-C
preP_psi = univ.selectAtoms( 'resnum ' + str(n_res) + ' and name N',
                            'resnum ' + str(n_res) + ' and name CA',
                            'resnum ' + str(n_res) + ' and name C',
                            'resnum ' + str(n_res+1) + ' and name N') # N-CA-C-N+1

ofile = file( pro_name_ + '_waterHB_Dih_V2.dat', 'w');
print >>ofile, 'frame'.ljust(10) + 'pdY_phi pdY_psi pdY_O pdY_chi1  N1_ND2 Yn_O Yn_chi1  wPf_phi wPf_psi pMd_phi pMd_psi preP_phi preP_psi';
for ts in univ.trajectory:
    print >>ofile, str(ts.frame).ljust(10) + ' ' + str(pdY_phi.dihedral()) + ' ' + \
          str(pdY_psi.dihedral()) + ' ' + str(pdY_O.dihedral())+ ' ' + str(pdY_chi.dihedral()) + ' ' + str(N1_ND2.dihedral()) + ' ' + str(Yn_O.dihedral()) + \
          ' ' + str(Yn_chi.dihedral()) + ' ' + str(wPf_phi.dihedral()) + ' ' + str(wPf_psi.dihedral())+ \
	' '+str(pMd_phi.dihedral()) + ' ' +str(pMd_psi.dihedral()) + ' ' + str(preP_phi.dihedral()) + ' ' +str(preP_psi.dihedral());
ofile.close();

