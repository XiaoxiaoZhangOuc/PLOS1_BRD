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
ifile = file('conserHolo_check.txt','r')
pro_marks = {};
for iline in ifile:
    iwords = iline.split();
    pro_marks[ iwords[0] ] = iwords[1:]

#print pro_marks;

pro_names = pro_marks.keys()

#pro_ = os.getcwd().split('/')[-1]
pro_ = sys.argv[1]; 
pro = pro_.split('.')[0];

print pro.ljust(7), 
prmfile = sys.argv[1];
trjfile = sys.argv[2];
univ = Universe( prmfile, trjfile );

mark1 = univ.selectAtoms('resid ' + pro_marks[pro][0] + ' and name O' )[0];
# resi wPf and atom name O
mark2 = univ.selectAtoms( 'resid ' + str( int( pro_marks[pro][1] )-2 )  + ' and name O' )[0];
# resi preP and atom name O
mark3 = univ.selectAtoms( 'resid ' + str( int( pro_marks[pro][4] )+1 ) + ' and name OH')[0];
# resi pdY and name OH
mark4 = univ.selectAtoms( 'resid ' + pro_marks[pro][5] + ' and name O' )[0];
mark5 = univ.selectAtoms( 'resid ' + pro_marks[pro][5] + ' and name N' )[0];
# resi pMd and name O, N
mark6 = univ.selectAtoms( 'resid ' + pro_marks[pro][8] + ' and name O' )[0];
# resi V1  Valine O
mark7 = univ.selectAtoms( 'resid ' + pro_marks[pro][9] + ' and name ND2' )[0];
# resi N1 and name ND2

#print pro_marks[pro][0], str( int( pro_marks[pro][1] )-2 ), str( int( pro_marks[pro][2] )+1 ), pro_marks[pro][3] , \
#      str( int( pro_marks[pro][6] ) - 7 ), str( int( pro_marks[pro][6] ) - 4 ), str( int( pro_marks[pro][6] ) - 4 );

write1 = False;  write2 = False;  write3 = True;  write4 = False;  write5 = False;

ofile1 = file( pro_+'_waterResi_distance.txt', 'w');
ofile1.write( 'frame'.ljust(8) + ' pdY_O-wPf_O pdY_O-pMd_O pdY_O-V1_O N1_ND2-V1_O pMd_N-V1_O preP_O-V1_O pMd_N-wPf_O '+ '\n' );

for ts in univ.trajectory:
    ofile1.write( str( ts.frame ).ljust(8)  + '%8.2f'  %lg.norm(mark1.position-mark3.position) );

    ofile1.write( '  ' + '%8.2f' %lg.norm(mark1.position-mark4.position) );
    ofile1.write( '  ' + '%8.2f' %lg.norm(mark1.position-mark6.position) );
    ofile1.write( '  ' + '%8.2f' %lg.norm(mark7.position-mark6.position) );

    ofile1.write( '  ' + '%8.2f' %lg.norm(mark5.position-mark6.position) );
    ofile1.write( '  ' + '%8.2f' %lg.norm(mark2.position-mark6.position) );

    ofile1.write( '  ' + '%8.2f' %lg.norm(mark5.position-mark1.position) + '\n' );
ofile1.close();



