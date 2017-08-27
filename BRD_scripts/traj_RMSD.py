#!/usr/bin/python
# print the RMSD of the trajectory

import sys;
import os;
from MDAnalysis import Universe
from sys import argv
from MDAnalysis import analysis

sys.path.append('/share/udata1/xiaoxiao/BRD_Ana_15.4.22/scripts/');

from pro_ss import pro_ss;

def RMSF2StartEndResi( RMSFfile ):
    cutoff = 1.5;
    ifile = file(RMSFfile, 'r');
    ilines = ifile.readlines();
    ifile.close();
    RMSF = [];
    for iline in ilines[1:]:
        iwords = iline.split();
        RMSF.append( float(iwords[1]) );
    start, end = 0, 0;
    for i in range(0,len(RMSF)-1,1):
        if RMSF[i] < cutoff and RMSF[i+1] < cutoff:
            start = i+1;
            break;

    for j in range(len(RMSF)-1, 0, -1):
        if RMSF[j] < cutoff and RMSF[j-1] < cutoff:
            end = j+1;
            break;
    return ( start, end );

if len(argv) < 3:
    print 'Usage: ./traj_RMSD.py reffile trjfile';
    sys.exit();

pro_name_ = os.getcwd().split('/')[-1]
pro_name = pro_name_.split('_run')[0];
(HZE, HAB, HAE, HBB, HBE, HCB, resi_len) = pro_ss( pro_name )

reffile = argv[1];
trjfile = argv[2];
(start_resi, end_resi) = RMSF2StartEndResi( pro_name_ + '_Ca_RMSF.dat' );
ref  = Universe( reffile );
start_resi = ref.residues[0].resnum-1 + start_resi;
end_resi = ref.residues[0].resnum-1 + end_resi;

#print (start_resi, end_resi);


ofile = file( pro_name_+'_RMSD.dat', 'w');
print >>ofile, 'frame  pro_Calpha  resi_'+str(start_resi)+'-'+str(end_resi) + \
      '   resi_1:'+str(start_resi-1)+','+str(end_resi+1)+':'+str(resi_len) + '   '+\
      'HZABC      loop_ZA      loop_AB      loop_BC';

univ = Universe( reffile, trjfile );
for ts in univ.trajectory:
    ofile.write( str(ts.frame).ljust(8) );
    ofile.write( '%.4f' %(analysis.align.alignto(univ, ref, select='name CA')[1]) );
    ofile.write( '       ');

    mobile1 = univ.selectAtoms( "resid "+str(start_resi)+':'+str(end_resi) + ' and name CA' )
    ref1  = ref.selectAtoms( "resid "+str(start_resi)+':'+str(end_resi) + ' and name CA' ) 
    ofile.write( '%.4f' %(analysis.align.rmsd( mobile1.coordinates(), ref1.coordinates() )) );
    ofile.write( '           ');
    
    #mobile2 = univ.selectAtoms( "(resid "+str(start_resi)+':'+str(HZE) +' or resid '+str(HAB)+':'+\
    #           str(HAE)+' or resid '+str(HBB)+':'+str(HBE) +' or resid '+str(HCB)+':'+str(end_resi) + ') and name CA')
    #ref2  = ref.selectAtoms( "(resid "+str(start_resi)+':'+str(HZE) +' or resid '+str(HAB)+':'+\
    #           str(HAE)+' or resid '+str(HBB)+':'+str(HBE) +' or resid '+str(HCB)+':'+str(end_resi) + ') and name CA')
    #ofile.write( '%.4f' %(analysis.align.rmsd( mobile2.coordinates(), ref2.coordinates() ) ) );
    #ofile.write( '           ');
    
    mobile3  = univ.selectAtoms( "(resid "+str(HZE)+':'+str(HAB) +' ) and name CA');
    ref3  = ref.selectAtoms( "(resid "+str(HZE)+':'+str(HAB) +' ) and name CA' );
    ofile.write( '%.4f' %(analysis.align.rmsd( mobile3.coordinates(), ref3.coordinates() ) ) );
    ofile.write( '           ');
    
    #mobile4 = univ.selectAtoms( "(resid "+str(HAE)+':'+str(HBB) + ') and name CA');
    #ref4 = ref.selectAtoms( "(resid "+str(HAE)+':'+str(HBB) + ') and name CA');
    #ofile.write( '%.4f' %(analysis.align.rmsd( mobile4.coordinates(), ref4.coordinates() ) ) );
    #ofile.write( '           ');    

    #mobile5 = univ.selectAtoms( "(resid "+str(HBE)+':'+str(HCB) + ') and name CA');
    #ref5 = ref.selectAtoms( "(resid "+str(HBE)+':'+str(HCB) + ') and name CA');
    #ofile.write( '%.4f' %(analysis.align.rmsd( mobile5.coordinates(), ref5.coordinates() ) ) );
    ofile.write( '\n' );

ofile.close();
