#!/usr/bin/python
# Author: Kevin Chan; Data: 2014/2/20;
# get the pocket water of BRD MD;
# Version 2: 2/21/2014
# Modified Data: 2014/3/11
# Modified Data: 2014/6/6

import sys;
import os;
import MDAnalysis;
import numpy as np;
import numpy.linalg as lg;
from sys import argv;

sys.path.append('/share/udata1/xiaoxiao/BRD_Ana_15.4.22/scripts/');
from pro_conser_marks import pro_conser_marks;

def traj_wat_ana(gro,xtc):
    '''
    Input: top file and coord file needed by MDAnalysis.Universe
    Output: pocket water analysis results;
    '''
    #start = 1800;
    #end = 1900;
    pro_name_ = os.getcwd().split('/')[-1];
    pro_name = pro_name_.split('_run')[0];
    traj = MDAnalysis.Universe(gro,xtc);
    marks = pro_conser_marks(pro_name);

    # res id of conserve residue at protein pocket 
    wPf_O = traj.selectAtoms('resnum '+ marks['wPf'] + ' and name O' )[0];
    preP_O = traj.selectAtoms('resnum '+str(int(marks['pVd'])-2)+' and name O')[0];

    pdY_OH = traj.selectAtoms('resnum '+str(int(marks['pDy'])+1)+' and name OH')[0];
   # pdY_OH = traj.selectAtoms('resnum '+str(int(marks['pDy'])+1)+' and name HZ')[0]; # for modified brds

    pMd_N = traj.selectAtoms('resnum '+marks['pMd']+' and name N' )[0];
    pMd_O = traj.selectAtoms('resnum '+marks['pMd']+' and name O' )[0];
    pmD_CA = traj.selectAtoms('resnum '+str(int(marks['pMd'])+1)+' and name CA')[0];

    V1_O = traj.selectAtoms('resnum '+marks['V1']+' and name O')[0]

    N1_ND2 = traj.selectAtoms('resnum '+marks['N1']+' and name ND2')[0]
    N1_CB = traj.selectAtoms('resnum '+marks['N1']+' and name CB')[0]
    N1_O = traj.selectAtoms('resnum '+marks['N1']+' and name O')[0];
    N1p1_CA = traj.selectAtoms('resnum '+str(int(marks['N1'])+1)+' and name CA')[0];

    coords1, coords2, coords3, coords4, coords5 = [], [], [], [], [];
    

    wat1_id = file("wat1_B.txt",'w');
    wat2_id = file("wat2_B.txt",'w');
    wat3_id = file("wat3_B.txt",'w');
    wat4_id = file("wat4_B.txt",'w');
    wat5_id = file("wat5_B.txt",'w');
    
    #print pro_name,
    for ts in traj.trajectory: #[start:end]:
        wats1 = wats2 = wats3 = wats4 = wats5 = traj.selectAtoms( 'name XYZ');
        wats = traj.selectAtoms('name O* and resname SOL and around 7 ( (resnum '\
               + marks['wPf'] + ' or resnum ' + str(int(marks['pVd'])-2) + ' or resnum '\
               + marks['pVd'] + ' or resnum ' + str(int(marks['pDy'])+1) + ' or resnum '\
               + marks['pMd'] + ' or resnum ' + str(int(marks['pMd'])+1) + ' or resnum ' \
               + marks['V1']  + ' or resnum ' + marks['N1'] + ' or resnum ' + \
               str(int(marks['N1'])+1) + ' or resnum ' + str(int(marks['Yn'])+1) + ') and \
               (name O* or name N*) )' );
        wat1_id.write( str(ts.frame) + '   ' );
        wat2_id.write( str(ts.frame) + '   ' );
        wat3_id.write( str(ts.frame) + '   ' );
        wat4_id.write( str(ts.frame) + '   ' );
        wat5_id.write( str(ts.frame) + '   ' );
	#if len(wats)print len(wats)
        for wat in wats:

            if (wat+V1_O).bond()<3.5 and (wat+N1_CB).bond()<5.0:
                if len(wats3)!=0 and (wat+N1_CB).bond() < (wats3[0]+N1_CB).bond():
                    wats3 = traj.selectAtoms( 'name XYZ');
                    wats3 += wat;
                else:
                    wats3 += wat;

            if (wat+pMd_N).bond() < 4.0 and (wat+N1_ND2).bond() < 4.0:
                wats1 += wat;
            if len(wats1)==0 and len(wats3) !=0 and 2.0<(wat+wats3[0]).bond()<5.0 and (wat+N1_ND2).bond() < 4.0\
                and (wat+pMd_N).bond() < 4.0:
                wats1 += wat

            if  (wat+N1p1_CA).bond() < 5.0 and 3.5 < (wat+V1_O).bond() < 6.0 and \
               (wat+wPf_O).bond() <= (pdY_OH+wPf_O).bond() and (wat+pMd_O).bond()>3.5:
                if len(wats3)!=0:
                    if 2.0<(wat+wats3[0]).bond()<4.0:
                        wats2 += wat;
                    continue;
                if len(wats1)!=0:
                    if 2.0<(wat+wats1[0]).bond()<5.0:
                        wats2 += wat;
                    continue;
                wats2 += wat;

            if (wat+pMd_O).bond()<3.5 and (wat+pmD_CA).bond()<5.0 and 3.5 < (wat+V1_O).bond() < 6.0:
                if len(wats3)!=0:
                    if 2.0<(wat+wats3[0]).bond()<4.0:
                        if len(wats4)!=0:
                            wats4 = traj.selectAtoms( 'name XYZ');
                        wats4 += wat
                    continue;
                wats4 += wat;

            if (wat+wPf_O).bond()<4.5 and (wat+preP_O).bond()< 5.0  and (wat+pMd_O).bond() < 6.0:
                if len(wats4)!=0:
                    if 2.0<(wat+wats4[0]).bond()<4.0:
                        wats5 += wat
                    continue;
                wats5 += wat
	    

        if len(wats1) != 0:
            coords1.append( list(wats1[0].pos) );
            wat1_id.write( str(wats1[0].resid ) + ' ' );
        if len(wats2) != 0:
            coords2.append( list(wats2[0].pos) );
            wat2_id.write( str(wats2[0].resid ) + ' ' );
        if len(wats3) != 0:
            coords3.append( list(wats3[0].pos) );
            wat3_id.write( str(wats3[0].resid ) + ' ' );
        if len(wats4) != 0:
            coords4.append( list(wats4[0].pos) );
            wat4_id.write( str(wats4[0].resid ) + ' ' );
        if len(wats5) != 0:
            coords5.append( list(wats5[0].pos) );
            wat5_id.write( str(wats5[0].resid ) + ' ' );
        wat1_id.write( '\n' );
        wat2_id.write( '\n' );
        wat3_id.write( '\n' );
        wat4_id.write( '\n' );
        wat5_id.write( '\n' );
    
#traj_wat_ana( '../em.gro', '../md_fit.xtc')
#pro_name_ = os.getcwd().split('/')[-2];
#pro_name = pro_name_.split('_run')[0];
prmfile, trjfile = argv[1], argv[2];
traj_wat_ana(prmfile, trjfile)
