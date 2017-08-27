#!/usr/bin/python
# Get the RMSF of the trajectory
# Author: Kevin Chan
# Date: 2014-1-24

import sys;
from sys import argv;
import os;
import numpy as np;
from MDAnalysis import Universe;
from MDAnalysis import analysis;

freq = 10

def RMSF_traj( reffile, trajfile, startframe=0, endframe=-1):
    '''
    Get the RMSF from a trajectory
    Snapshots are aligned to the average structure
    '''
    traj = Universe( reffile, trajfile );
    ref  = Universe( reffile );
    ave  = ref.atoms.CA;

    cycle, count = 20, 0;
    rmsd, rmsd_ = 0.0, 0.0;
    if endframe == -1 or endframe > traj.trajectory.numframes:
        endframe = traj.trajectory.numframes;
    
    rmsds  = np.empty( ( (endframe-startframe)/freq+1, ) );
    coords = np.empty( ( (endframe-startframe)/freq+1,traj.atoms.CA.numberOfAtoms(),3) );
    sub = np.empty( ( (endframe-startframe)/freq+1,traj.atoms.CA.numberOfAtoms(),3) );

    while count < cycle:
        for ts in traj.trajectory[startframe:endframe:freq]:
            rmsds[ (ts.frame-startframe)/freq ] =  analysis.align.alignto( traj.atoms.CA, ave )[1];
            coords[ (ts.frame-startframe)/freq ] = traj.atoms.CA.positions
            sub[ (ts.frame-startframe)/freq ] = traj.atoms.CA.positions - ave.positions;
        rmsd = np.mean( rmsds );
        print rmsd;
        ave.set_positions( np.mean( coords, axis=0 ) );
        
        if np.abs( rmsd - rmsd_) < 0.00001:
            break;
        else:
            rmsd_ = rmsd;
        count += 1;
    squ = np.multiply(sub, sub);
    dis = np.sum( squ, axis=2 );
    MSF = np.mean( dis, axis=0 );
    RMSF = np.sqrt( MSF )
    return RMSF;

if __name__ == '__main__':
    pro_name_ = os.getcwd().split('/')[-1]
    pro_name = pro_name_.split('_run')[0];
    if len(argv) < 3:
        print "Usage: ./traj_ave_stru.py refstru trajstru startframe endframe";
        sys.exit();
    elif len(argv) == 3:
        RMSF = RMSF_traj( argv[1], argv[2] );
    elif len(argv) == 4:
        RMSF = RMSF_traj( argv[1], argv[2], int(argv[3]) );
    elif len(argv) == 5:
        RMSF = RMSF_traj( argv[1], argv[2], int(argv[3]), int(argv[4]) );
    ofile = file( pro_name_ + '_Ca_RMSF.dat', 'w');
    print >>ofile, 'resi'.ljust(8) + 'RMSF';
    for i in range( len(RMSF) ):
        print >>ofile, str(i+1).ljust(8) + str(RMSF[i]);
    ofile.close();

