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

def Tyr_CaCbOhHh(prmfile, trjfile, pro_name):
    '''
    return the dihedral angle of Yn residue: Ca-Cb__Oh-Hh
    '''
    univ = Universe( prmfile, trjfile );
    frames = np.empty( (univ.trajectory.numframes,) )
    Dih_Tyr_CaCbOhHh = np.empty( (univ.trajectory.numframes,) )

    if univ.residues[ int(pro_conser_marks(pro_name)['pDy'])+1-univ.residues[0].resnum ].name != 'TYR':
        return ( None, None );

    i_res = int(pro_conser_marks(pro_name)['pDy']) + 1;
    CaCbOhHh = univ.selectAtoms( 'resnum ' + str(i_res) + ' and name CA',
                                 'resnum ' + str(i_res) + ' and name CB',
    				 'resnum ' + str(i_res) + ' and name OH',
    				 'resnum ' + str(i_res) + ' and name HH' ) ;
    index = 0;
    for ts in univ.trajectory:
        frames[index] = ts.frame;
        Dih_Tyr_CaCbOhHh[index] = CaCbOhHh.dihedral();
        index += 1;
    return (frames, Dih_Tyr_CaCbOhHh);

def yN_chi(prmfile, trjfile, pro_name):
    '''
    return the chi1, chi2 of the conserved yN residue
    '''
    univ = Universe( prmfile, trjfile );
    frames = np.empty( (univ.trajectory.numframes,) )
    Dih_yN_chi1 = np.empty( (univ.trajectory.numframes,) );
    Dih_yN_chi2 = np.empty( (univ.trajectory.numframes,) );
    if univ.residues[ int(pro_conser_marks(pro_name)['Yn'])-univ.residues[0].resnum+1 ].name not \
    in [ 'ASN', 'TYR', 'THR']:
        return (None, None, None);
    i_res = int(pro_conser_marks(pro_name)['Yn']) + 1;
    if univ.residues[ int(pro_conser_marks(pro_name)['Yn'])-univ.residues[0].resnum+1].name == 'ASN':
        yN_chi1 = univ.selectAtoms( 'resnum ' + str(i_res) + ' and name N',
                                    'resnum ' + str(i_res) + ' and name CA',
                                    'resnum ' + str(i_res) + ' and name CB',
                                    'resnum ' + str(i_res) + ' and name CG' ) ;
        yN_chi2 = univ.selectAtoms( 'resnum ' + str(i_res) + ' and name CA',
                                    'resnum ' + str(i_res) + ' and name CB',
                                    'resnum ' + str(i_res) + ' and name CG',
                                    'resnum ' + str(i_res) + ' and name OD1' ) ;
    elif univ.residues[ int(pro_conser_marks(pro_name)['Yn'])-univ.residues[0].resnum+1].name == 'TYR':
        yN_chi1 = univ.selectAtoms( 'resnum ' + str(i_res) + ' and name N',
                                    'resnum ' + str(i_res) + ' and name CA',
                                    'resnum ' + str(i_res) + ' and name CB',
                                    'resnum ' + str(i_res) + ' and name CG' ) ;
        yN_chi2 = univ.selectAtoms( 'resnum ' + str(i_res) + ' and name CA',
                                    'resnum ' + str(i_res) + ' and name CB',
                                    'resnum ' + str(i_res) + ' and name OH',
                                    'resnum ' + str(i_res) + ' and name HH');
    elif univ.residues[ int(pro_conser_marks(pro_name)['Yn'])-univ.residues[0].resnum+1].name == 'THR':
        yN_chi1 = univ.selectAtoms( 'resnum ' + str(i_res) + ' and name N',
                                    'resnum ' + str(i_res) + ' and name CA',
                                    'resnum ' + str(i_res) + ' and name CB',
                                    'resnum ' + str(i_res) + ' and name OG1');
        yN_chi2 = univ.selectAtoms( 'resnum ' + str(i_res) + ' and name CA',
                                    'resnum ' + str(i_res) + ' and name CB',
                                    'resnum ' + str(i_res) + ' and name OG1',
                                    'resnum ' + str(i_res) + ' and name HG1');

    index = 0;
    for ts in univ.trajectory:
        frames[index] = ts.frame;
        Dih_yN_chi1[index] = yN_chi1.dihedral();
        Dih_yN_chi2[index] = yN_chi2.dihedral();
        index += 1;
    return (frames, Dih_yN_chi1, Dih_yN_chi2);

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print "Usage: ./traj_resiAna_pdY_CaCbOhHh.py prmfile trjfile";
        sys.exit();
    pro_name_ = os.getcwd().split('/')[-1];
    pro_name = pro_name_.split('_run')[0];
    print 'pro_name', pro_name;
    prmfile, trjfile = argv[1], argv[2];
    (frames, Dih_Tyr_CaCbOhHh) = Tyr_CaCbOhHh( prmfile, trjfile, pro_name );
    (frames_, Dih_yN_chi1, Dih_yN_chi2) = yN_chi( prmfile, trjfile, pro_name );

    ofile = file( pro_name_ + '_ConserResi_Dihedral.dat', 'w');
    if type(frames)==type(None) and type(frames_)==type(None):
        print "pdY not Tyrosine and yN not Asparagine"
        print "check the protein structure"
        sys.exit();
    elif type(frames)==type(None) and type(frames_)!=type(None) :
        print >>ofile, 'frame'.ljust(10) + '  ' + 'a_Dih_yN_chi1    a_Dih_yN_chi2';
        for i in range( len(frames_) ):
            print >>ofile, str(frames_[i]).ljust(10) + '  ' + str(Dih_yN_chi1[i]) + '  ' + str(Dih_yN_chi2[i]);
    elif type(frames)!=type(None) and type(frames_)==type(None) :
        print >>ofile, 'frame'.ljust(10) + '  ' + 'a_Dih_Tyr_CaCbOhHh';
        for i in range( len(frames) ):
            print >>ofile, str(frames[i]).ljust(10) + '  ' + str(Dih_Tyr_CaCbOhHh[i]);
    else:
        print >>ofile, 'frame'.ljust(10) + '  ' + 'a_Dih_Tyr_CaCbOhHh' + '   ' + \
               'a_Dih_yN_chi1    a_Dih_yN_chi2';
        for i in range( len(frames) ):
            print >>ofile, str(frames[i]).ljust(10) + '  ' + str(Dih_Tyr_CaCbOhHh[i]) + '  ' + \
                  str(Dih_yN_chi1[i]) + '  ' + str(Dih_yN_chi2[i]);  
    ofile.close();




