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
ifile = file('/share/udata1/xiaoxiao/BRD_Ana_15.4.22/scripts/conserved_resi_v2.txt','r');
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

mark1 = univ.selectAtoms( 'resid ' + pro_marks[pro][4] + ' and name CA' )[0];
# resi wPf and atom name O
mark2 = univ.selectAtoms( 'resid ' + pro_marks[pro][3] + ' and name CA' )[0];
# resi pVd-2 and atom name O
mark3 = univ.selectAtoms( 'resid ' + pro_marks[pro][2] + ' and name CA')[0];
# resi pdY and name OH
mark4 = univ.selectAtoms( 'resid ' + pro_marks[pro][1] + ' and name CA' )[0];
# resi pMd and name O
mark5 = univ.selectAtoms( 'resid ' + pro_marks[pro][10] + ' and name CA' )[0];
# resi N1-3  Valine O
mark6 = univ.selectAtoms( 'resid ' + pro_marks[pro][11] + ' and name CA' )[0];
# resi N1 and name ND2
mark7 = univ.selectAtoms( 'resid ' + pro_marks[pro][0] + ' and name CA' )[0];
# resi N1 and name O

#print pro_marks[pro][0], str( int( pro_marks[pro][1] )-2 ), str( int( pro_marks[pro][2] )+1 ), pro_marks[pro][3] , \
#      str( int( pro_marks[pro][6] ) - 7 ), str( int( pro_marks[pro][6] ) - 4 ), str( int( pro_marks[pro][6] ) - 4 );

write1 = False;  write2 = False;  write3 = True;  write4 = False;  write5 = False;

ofile1 = file( pro_+'_ZACa_distance_HBC.txt', 'w');
ofile1.write( 'frame'.ljust(8) + ' pDy-Yn   I94-Yn   I94-HC1   L92-wPf   L92-HC1   pVd-wPf   pVd-Yn   wPf-HC1 '+ '\n' );

for ts in univ.trajectory:
    ofile1.write( str( ts.frame ).ljust(8)  + '%8.2f'  %lg.norm(mark1.position-mark2.position) );
    
    ofile1.write( '  ' + '%8.2f' %lg.norm(mark1.position-mark5.position) );
    ofile1.write( '  ' + '%8.2f' %lg.norm(mark2.position-mark5.position) );
    ofile1.write( '  ' + '%8.2f' %lg.norm(mark2.position-mark6.position) );

    ofile1.write( '  ' + '%8.2f' %lg.norm(mark3.position-mark7.position) );
    ofile1.write( '  ' + '%8.2f' %lg.norm(mark3.position-mark6.position) );

    ofile1.write( '  ' + '%8.2f' %lg.norm(mark4.position-mark7.position) );
    ofile1.write( '  ' + '%8.2f' %lg.norm(mark4.position-mark5.position) );

    ofile1.write( '  ' + '%8.2f' %lg.norm(mark7.position-mark6.position) + '\n' );
ofile1.close();



