#!/usr/bin/python

from sys import argv;
from MDAnalysis import Universe;
import numpy.linalg as lg;
from glob import glob;
import sys;
import os;

# change the res_range to all the residues 2--end_resi-1
resn_short = { 'GLY':'G', 'ALA':'A', 'SER':'S', 'THR':'T', 'CYS':'C', \
               'VAL':'V', 'LEU':'L', 'ILE':'I', 'MET':'M', 'PRO':'P', \
	       'PHE':'F', 'TYR':'Y', 'TRP':'W', 'ASP':'D', 'GLU':'E', \
	       'ASN':'N', 'GLN':'Q', 'HIS':'H', 'HID':'H', 'HIE':'H', \
	       'LYS':'K', 'ARG':'R' }
Resi_CG = [ 'ARG', 'LYS', 'ASP', 'GLU', 'ASN', 'GLN', 'LEU', 'MET', 'PHE', 'TYR',\
            'TRP', 'HIS', 'HID', 'HIE' ]
Resi_CG1 = ['VAL', 'ILE']
Resi_OG = ['SER',]
Resi_OG1 = ['THR']
Resi_SG = ['CYS']

# get the conserved residues' index
ifile = file('/share/udata1/xiaoxiao/BRD_Ana_15.4.22/scripts/conserved_resi_v3.txt','r');
pro_marks = {};
for iline in ifile:
    iwords = iline.split();
    if len(iwords) != 13:
        continue;
    if len(iwords) == 13:
        pro_marks[ iwords[0] ] = iwords[1:]

#print pro_marks;

pro_names = pro_marks.keys()

pro_ = os.getcwd().split('/')[-1]
pro = pro_.split('_run')[0];

print pro.ljust(7), 
prmfile = sys.argv[1];
trjfile = sys.argv[2];
univ = Universe( prmfile, trjfile );

mark1 = univ.selectAtoms( 'resid ' + str(int(pro_marks[pro][4])+1) + ' and name CA' )[0];
# resi pDy
mark2 = univ.selectAtoms( 'resid ' + pro_marks[pro][3] + ' and name CA' )[0];
# resi I94
mark3 = univ.selectAtoms( 'resid ' + str(int(pro_marks[pro][1])-2) + ' and name CA')[0];
# resi preP
mark4 = univ.selectAtoms( 'resid ' + pro_marks[pro][1] + ' and name CA' )[0];
# resi pVd
mark5 = univ.selectAtoms( 'resid ' + pro_marks[pro][4] + ' and name CA' )[0];
# resi L92
mark6 = univ.selectAtoms( 'resid ' + str(int(pro_marks[pro][11])+1) + ' and name CA' )[0];
# resi yN
mark7 = univ.selectAtoms( 'resid ' + pro_marks[pro][0] + ' and name CA' )[0];
# resi wPf 

#print pro_marks[pro][0], str( int( pro_marks[pro][1] )-2 ), str( int( pro_marks[pro][2] )+1 ), pro_marks[pro][3] , \
#      str( int( pro_marks[pro][6] ) - 7 ), str( int( pro_marks[pro][6] ) - 4 ), str( int( pro_marks[pro][6] ) - 4 );

write1 = False;  write2 = False;  write3 = True;  write4 = False;  write5 = False;

ofile1 = file( pro_+'_distance_yN.txt', 'w');
ofile1.write( 'frame'.ljust(8) + ' wPf-yN   preP-yN  pVd-yN L92-yN   I94-yN   pdY-yN '+ '\n' );

for ts in univ.trajectory:
    ofile1.write( str( ts.frame ).ljust(8)  + '%8.2f'  %lg.norm(mark7.position-mark6.position) );
    
    ofile1.write( '  ' + '%8.2f' %lg.norm(mark3.position-mark6.position) );
    ofile1.write( '  ' + '%8.2f' %lg.norm(mark4.position-mark6.position) );
    ofile1.write( '  ' + '%8.2f' %lg.norm(mark5.position-mark6.position) );
    ofile1.write( '  ' + '%8.2f' %lg.norm(mark2.position-mark6.position) );

    ofile1.write( '  ' + '%8.2f' %lg.norm(mark1.position-mark6.position) + '\n' );
ofile1.close();



