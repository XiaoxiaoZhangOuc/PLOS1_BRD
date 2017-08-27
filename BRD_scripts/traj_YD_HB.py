#!/usr/bin/python
from sys import argv;
from MDAnalysis import Universe;
import numpy as np;
import numpy.linalg as lg;
from glob import glob;
import sys;
import os;

sys.path.append('/share/udata1/xiaoxiao/BRD_Ana_15.4.22/scripts/');
from pro_conser_marks import pro_conser_marks;

def YD_HB( Tyr, Asp, traj ):
    '''
    return the hydrogen bond between Y-HH and D-OD
    '''
    Y_HH = Tyr.selectAtoms('name HH')[0]
    D_OD1 = Asp.selectAtoms('name OD1')[0]
    D_OD2 = Asp.selectAtoms('name OD2')[0]
    HB = np.empty( (traj.trajectory.numframes,3) )
    for ts in traj.trajectory:
        HB[ ts.frame-1, 0 ] = ts.frame;
        HB[ ts.frame-1, 1 ] = lg.norm( Y_HH.position - D_OD1.position );
        HB[ ts.frame-1, 2 ] = lg.norm( Y_HH.position - D_OD2.position );
    return HB;

if __name__ == '__main__':
    pro_name_ = os.getcwd().split('/')[-1];
    pro_name = pro_name_.split('_run')[0];
    print 'pro_name', pro_name,
    prmfile = argv[1];
    trjfile = argv[2];
    univ = Universe( prmfile, trjfile );
    print univ.residues[int(pro_conser_marks(pro_name)['Yab'])-univ.residues[0].resnum].name, \
          univ.residues[int(pro_conser_marks(pro_name)['Db'])-univ.residues[0].resnum].name;
    if univ.residues[int(pro_conser_marks(pro_name)['Yab'])-univ.residues[0].resnum ].name == 'TYR' and \
       univ.residues[int(pro_conser_marks(pro_name)['Db'])-univ.residues[0].resnum ].name == 'ASP':
        Tyr = univ.selectAtoms( 'resnum ' + str( int(pro_conser_marks(pro_name)['Yab'])) )
        Asp = univ.selectAtoms( 'resnum ' + str( int(pro_conser_marks(pro_name)['Db'])) )
        HB1 = YD_HB( Tyr, Asp, univ );
        ofile = file( pro_name_ + '_Y119_D128_HB.dat', 'w');
        ofile.write( 'frame        Y-HH_D-OD1     Y-HH_D-OD2\n' );
        for each in HB1:
            ofile.write( str(each[0]).ljust(12)+' '+'%8.2f' %each[1]+' '+'%8.2f' %each[2] + '\n' );
        ofile.close();

