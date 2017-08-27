#!/usr/local/bin/python
# tackle the RCSB PDB file, get the useful information from it.
# Kevin Chan
# 2012/9/11. Version 1.1

import sys;
from glob import glob;
import os;

Str_Break = ' ' + '$' + ' '

def ReadIn( filename ):
    Date = '';
    Technique = '';
    PDB_ID = '';
    Tit = '';
    Reso = '';
    R = '';
    R_free = '';
    Ref = '';
    J_Tit = '';
    PH = '';
    Temp = '';
    Env = '';
    
    PDB_ID = filename.split('.')[0];
    ifile = file( filename, 'r');
    ilines = ifile.readlines();
    ifile.close();
    try:
        Date = ilines[0].split()[-2];
    except:
        pass;
    for iline in ilines:
        iwords = iline.split();
        if len(iwords) > 1 and iwords[0] == 'TITLE':
            Tit += iline[5:].strip() + ' ';
        if len(iwords) == 7 and iline[0:10] =='REMARK 200' and iwords[2] == 'EXPERIMENT' and iwords[3] == 'TYPE':
            Technique = iwords[-2];
        elif len(iwords) == 6 and iline[0:10] =='REMARK 210' and iwords[2] == 'EXPERIMENT' and iwords[3] == 'TYPE':
            Technique = iwords[-1];
            
        if len(Reso)==0 and len(iwords) == 5 and iwords[0] == 'REMARK' and iwords[2] == 'RESOLUTION.':
            Reso = iwords[3];
        if len(R) == 0 and len(iline.strip()) > 49 and iline[0:20] =='REMARK   3   R VALUE':
            R = iline[48:].strip();
            try:
                float(R) > 0;
            except:
                R = '';     
        if len(R_free) == 0 and len(iline.strip()) > 49 and iline[0:28]=='REMARK   3   FREE R VALUE   ':
            R_free = iline[48:].strip();
        if len(iwords) > 2 and iwords[0] =='JRNL' and iwords[1]=='TITL':
            J_Tit += iline[31:].strip() + ' ';
        if len(iwords) >= 7 and iwords[0] == 'JRNL' and iwords[1] == 'REF':
            Ref = iwords[2] + ' ' + iwords[-1] + ' ' + 'V.' + iwords[-3] + ' ' + iwords[-2];
        if len(iwords)>3 and iline[0:10] =='REMARK 280':
            Env += iline[10:].strip() + ' ';
        for eachwords in Env.split(','):
            if len( eachwords.split() ) == 2 and eachwords.split()[0] == 'PH':
                PH = eachwords.split()[1];
            if len( eachwords.split() ) == 2 and eachwords.split()[0] == 'TEMPERATURE':
                Temp = eachwords.split()[1];

    Data_Summary = PDB_ID.center(15) + Str_Break + Technique.center(15) + Str_Break + Date.center(15) + Str_Break +  Reso.center(8) + Str_Break + R.center(8) + Str_Break + \
                   R_free.center(8) + Str_Break + Ref.center(50) + Str_Break + PH.center(8) + Str_Break + Temp.center(8) + Str_Break +\
                   Tit.center(30) + Str_Break + J_Tit.center(30) + '\n';
    return Data_Summary;

    
if len(sys.argv)>1:
    Files = sys.argv[1:]
else:
    Files = glob('*/*/*.pdb')
    print Files;

dirnames = glob('subclass_*');

ofile = file( 'PDB_Data_Summary.txt', 'w');
ofile.write( 'PDB_ID'.center(10) + Str_Break + 'Technique'.center(15) + Str_Break + 'Deposition_Date'.center(15) + Str_Break +  'Reso'.center(8) + Str_Break +'R'.center(8) \
             + Str_Break + 'R_free'.center(8) + Str_Break + 'Ref'.center(50) + Str_Break + 'PH'.center(8) + Str_Break + 'Temp'.center(8) + Str_Break + \
             'Tit'.center(30) + Str_Break + 'J_Tit'.center(30) + '\n' )

for eachfile in Files:
    Data_Summary = ReadIn( eachfile );
    ofile.write( Data_Summary );
ofile.close();

