#!/usr/bin/python
# Get the average structure of the trajectory
# Author: Kevin Chan
# Date: 2014-1-24

import sys;
from sys import argv;
import os;
import numpy as np;
from MDAnalysis import Universe;
from MDAnalysis import analysis;

freq = 20

def ave_stru( reffile, trajfile, startframe=0, endframe=-1):
    '''
    Get the average structure from a trajectory
    '''
    traj = Universe( reffile, trajfile );
    ref  = Universe( reffile );
    ave  = ref.atoms;

    cycle, count = 20, 0;
    rmsd, rmsd_ = 0.0, 0.0;
    if endframe == -1 or endframe > traj.trajectory.numframes:
        endframe = traj.trajectory.numframes;
    
    rmsds  = np.empty( ( (endframe-startframe)/freq+1, ) );
    coords = np.empty( ( (endframe-startframe)/freq+1,traj.atoms.numberOfAtoms(),3) );

    while count < cycle:
        for ts in traj.trajectory[startframe:endframe:freq]:
            rmsds[ (ts.frame-startframe)/freq ] =  analysis.align.alignto( traj, ave, select='name CA')[1];
            coords[ (ts.frame-startframe)/freq ] = traj.atoms.positions
    
        rmsd = np.mean( rmsds );
        print rmsd;
        ave.set_positions( np.mean( coords, axis=0 ) );
        
        if np.abs( rmsd - rmsd_) < 0.00001:
            break;
        else:
            rmsd_ = rmsd;
    return ( rmsds, ave );



if __name__ == '__main__':
    pro_name_ = os.getcwd().split('/')[-1]
    pro_name = pro_name_.split('_run')[0];
    if len(argv) < 3:
        print "Usage: ./traj_ave_stru.py refstru trajstru startframe endframe";
        sys.exit();
    elif len(argv) == 3:
        ( rmsds, avestru ) = ave_stru( argv[1], argv[2] );
    elif len(argv) == 4:
        ( rmsds, avestru ) = ave_stru( argv[1], argv[2], int(argv[3]) );
    elif len(argv) == 5:
        ( rmsds, avestru ) = ave_stru( argv[1], argv[2], int(argv[3]), int(argv[4]) );

    avestru.atoms.write( pro_name_ + '_AveStru.pdb') 
    

