import sys;
from sys import argv;
import os;
import MDAnalysis;
import numpy as np;
import numpy.linalg as lg;
from sys import argv;

sys.path.append('/share/udata1/xiaoxiao/BRD_Ana_15.4.22/scripts/');
from pro_conser_marks import pro_conser_marks;
def BRD_sel(pdb,xtc):
    '''
    Input: top file and coord file needed by MDAnalysis.Universe
    Output: pocket water analysis results;
    '''
    #start = 1800;
    #end = 1900;
    pro_name = str(pdb)[:-4];
    traj = MDAnalysis.Universe(pdb,xtc);
    marks = pro_conser_marks(pro_name);
    wPf_O= '(( resnum '+ str(marks['wPf']) + ' and name O) or '
    preP_O = '(resnum '+str(int(marks['pVd'])-2)+' and name O) or '
    pdY_OH = '(resnum '+str(int(marks['pDy'])+1)+' and name OH) or '
    pMd_N = '(resnum '+str(marks['pMd'])+' and name N) or '
    pMd_O = '(resnum '+str(marks['pMd'])+' and name O) or '
    V1_O = '(resnum '+str(marks['V1'])+' and name O ))' 
    key_res= wPf_O + preP_O + pdY_OH + pMd_N + pMd_O + V1_O
    print key_res
    
    wats = 'resname SOL and around 7 ( (resnum '\
               + marks['wPf'] + ' or resnum ' + str(int(marks['pVd'])-2) + ' or resnum '\
               + marks['pVd'] + ' or resnum ' + str(int(marks['pDy'])+1) + ' or resnum '\
               + marks['pMd'] + ' or resnum ' + str(int(marks['pMd'])+1) + ' or resnum ' \
               + marks['V1']  + ' or resnum ' + marks['N1'] + ' or resnum ' + \
               str(int(marks['N1'])+1) + ' or resnum ' + str(int(marks['Yn'])+1) + ') and \
               (name O* or name N*) )'
    print wats
    layer1 = 'resname SOL and around 3.5 '+key_res
    layer2 = 'resname SOL and around 3.5 protein and not around 3.5 '+ key_res 
    #layer3 = 'name OW and sphlayer 3.0 6.0 protein'
    return key_res, layer1, layer2, wats

import MDAnalysis
from MDAnalysis.analysis.waterdynamics import HydrogenBondLifetimes as HBL
import MDAnalysis.analysis.hbonds

pdb=argv[1]
trj=argv[2]
import MDAnalysis
from MDAnalysis.analysis.align import *
from MDAnalysis.analysis.rms import rmsd
from MDAnalysis.analysis.density import density_from_Universe
ref = MDAnalysis.Universe(pdb,pdb)
u = MDAnalysis.Universe(pdb,trj)
#rms_fit_trj(u, ref, select='protein')

#trj="rmsfit_"+argv[2]
u = MDAnalysis.Universe(pdb, trj)
# Define selections
key_res, layer1, layer2,wats = BRD_sel(pdb, trj)

HBL_layer1 = HBL(u, layer1, key_res, 0, 200, 10)
#HBL_layer1.run()
#np.save("HBLifetime_layer1_pro.npy",HBL_layer1.timeseries)
ofile=file("HBLifetime_layers.txt",'a')
#i=0
#ofile.write("HBL_layer1_pro")
#for HBLc, HBLi in HBL_layer1.timeseries:
#      ofile.write(str(i) +" "+ str(HBLc) +" "+ str(i) +" "+ str(HBLi)+'\n')
#      i+=1
#
#HBL_layer1_wat = HBL(u, layer1, layer1, 0, 200, 10)
#HBL_layer1_wat.run()
#np.save("HBLifetime_layer1_wat.npy",HBL_layer1_wat.timeseries)
#ofile=file("HBLifetime_layers.txt",'a')
#i=0
#ofile.write("HBL_layer1_wat")
#for HBLc, HBLi in HBL_layer1_wat.timeseries:
#      ofile.write(str(i) +" "+ str(HBLc) +" "+ str(i) +" "+ str(HBLi)+'\n')
#      i+=1

#HBL_layer1_all = HBL(u, layer1, "("+layer1 + ") or (" + key_res +")", 0, 200, 10)
#HBL_layer1_all.run()
#np.save("HBLifetime_layer1_all.npy",HBL_layer1_all.timeseries)
#ofile=file("HBLifetime_layers.txt",'a')
#i=0
#ofile.write("HBL_layer1_all")
#for HBLc, HBLi in HBL_layer1_all.timeseries:
#      ofile.write(str(i) +" "+ str(HBLc) +" "+ str(i) +" "+ str(HBLi)+'\n')
#      i+=1

#HBL_layer2 = HBL(u, layer2, 'protein', 0, 200, 10)
#HBL_layer2.run()
#np.save("HBLifetime_layer2_pro.npy",HBL_layer2.timeseries)
#ofile=file("HBLifetime_layers.txt",'a')
#i=0
#ofile.write("HBL_layer2_pro")
#for HBLc, HBLi in HBL_layer2.timeseries:
#      ofile.write(str(i) +" "+ str(HBLc) +" "+ str(i) +" "+ str(HBLi)+'\n')
#      i+=1

#HBL_layer2_wat = HBL(u, layer2, layer2, 0, 200, 10)
#HBL_layer2_wat.run()
#np.save("HBLifetime_layer2_wat.npy",HBL_layer2_wat.timeseries)
#ofile=file("HBLifetime_layers.txt",'a')
#i=0
#ofile.write("HBL_layer2_wat")
#for HBLc, HBLi in HBL_layer2_wat.timeseries:
#      ofile.write(str(i) +" "+ str(HBLc) +" "+ str(i) +" "+ str(HBLi)+'\n')
#      i+=1

#HBL_layer2_all = HBL(u, layer2, "("+layer2 + ") or protein" , 0, 200, 10)
#HBL_layer2_all.run()
#np.save("HBLifetime_layer2_all.npy",HBL_layer2_all.timeseries)
#ofile=file("HBLifetime_layers.txt",'a')
#i=0
#ofile.write("HBL_layer2_all")
#for HBLc, HBLi in HBL_layer2_wat.timeseries:
#      ofile.write(str(i) +" "+ str(HBLc) +" "+ str(i) +" "+ str(HBLi)+'\n')
#      i+=1

#
#Water Orientational Relaxation
#from MDAnalysis.analysis.waterdynamics import WaterOrientationalRelaxation as WOR

#WOR_analysis = WOR(u, wats, 0, 200, 10)
#WOR_analysis.run()
#i=0
#ofile.write('WOR_OH WOR_HH WOR_dip vs t'+'\n')
#print len(WOR_analysis.timeseries)
#now we print the data ready to graph. The first two columns are WOR_OH vs t graph,
#the second two columns are WOR_HH vs t graph and the third two columns are WOR_dip vs t graph
#for WOR_OH, WOR_HH, WOR_dip in WOR_analysis.timeseries:
#      ofile.write( str(i) +" "+ str(WOR_OH) +" "+ str(i) +" "+ str(WOR_HH) +" "+ str(i) +" "+ str(WOR_dip)+'\n')
#      i+=1

# Angular Distribution
#from MDAnalysis.analysis.waterdynamics import AngularDistribution as AD

#bins = 20
#AD_analysis = AD(u,wats,bins)
#AD_analysis.run()
#now we print data ready to graph. The first two columns are P(cos(theta)) vs cos(theta) for OH vector ,
#the seconds two columns are P(cos(theta)) vs cos(theta) for HH vector and thirds two columns
#are P(cos(theta)) vs cos(theta) for dipole vector
#i=0
#ofile.write('P(cos(theta)) vs cos(theta) for OH vector, HH vector, dipole vector'+'\n')
#for i in range(bins):
#    ofile.write(str(AD_analysis.graph[0][i]) +" "+ str(AD_analysis.graph[1][i]) +" "+ str(AD_analysis.graph[2][i])+'\n')
#    i+=1


# MeanSquareDisplacement

from MDAnalysis.analysis.waterdynamics import MeanSquareDisplacement as MSD

MSD_analysis = MSD(u, layer2, 0, 200, 20)
MSD_analysis.run()
#now we print data ready to graph. The graph
#represents MSD vs t
ofile.write('MSD vs t'+'\n')
i=0
for msd in MSD_analysis.timeseries:
    ofile.write(str(i) +" "+ str(msd)+'\n')
    i += 1

#
# Survival Probability
from MDAnalysis.analysis.waterdynamics import SurvivalProbability as SP

SP_analysis = SP(u, wats, 0, 200, 20)
SP_analysis.run()
#now we print data ready to graph. The graph
#represents SP vs t
ofile.write('SP vs t'+'\n')
i=0
for sp in SP_analysis.timeseries:
    ofile.write(str(i) +" "+ str(sp)+'\n')
    i += 1


