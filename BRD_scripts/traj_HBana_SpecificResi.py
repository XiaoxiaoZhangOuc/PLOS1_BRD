# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 14:28:32 2015

@author: xiaox
"""
import sys;
import MDAnalysis.analysis.hbonds
import os;
import MDAnalysis as mda;
import numpy as np;
import numpy.linalg as lg;
from sys import argv;
#%%
sys.path.append('/share/udata1/xiaoxiao/BRD_Ana_15.4.22/scripts/');
from pro_conser_marks import pro_conser_marks;
def BRD_sel(pdb,xtc):
    from MDAnalysis.analysis.waterdynamics import HydrogenBondLifetimes as HBL
    import MDAnalysis.analysis.hbonds
    '''
    Input: top file and coord file needed by MDAnalysis.Universe
    Output: pocket water analysis results;
    '''
    #start = 1800;
    #end = 1900;
    pro_name_ = str(pdb)[:-4];
    pro_name = pro_name_.split('_run')[0];
    marks = pro_conser_marks(pro_name);
    # res id of conserve residue at protein pocket
    pdY_ = str(int(marks['pDy'])+2);
    #Pdy = str(int(marks['pDy'])-1);
    #pDy = str(int(marks['pDy']));
    pvD_ = str(int(marks['pVd'])-2);
    sel_all = 'resnum '+ pvD_ + '-' + pdY_ + ' and backbone'
    return sel_all 
pdb=argv[1]
xtc=argv[2]
u = mda.Universe(pdb, xtc)
# Define selections
sel_all = BRD_sel(pdb, xtc)
h_all = mda.analysis.hbonds.HydrogenBondAnalysis(u, sel_all, sel_all, update_selection1=False,update_selection2=False)
h_all.run()
np.save(pdb[:-4]+'_HB_PVDPDY_bb.npy', h_all.timeseries)
