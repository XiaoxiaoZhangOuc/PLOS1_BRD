#!/usr/bin/python
# Analyse the wat in pocket during the trajectory;

import MDAnalysis;
from MDAnalysis import Universe;
from glob import glob;
from sys import argv;
import sys;
import os;
import numpy.linalg as lg;

# get the mask id;
ifile = file('/share/udata1/xiaoxiao/BRD_Ana_15.4.22/scripts/conserved_resi.txt','r');
pro_marks = {};
for iline in ifile:
    iwords = iline.split();
    if len(iwords) != 8:
        continue;
    if len(iwords) == 8:
        pro_marks[ iwords[0] ] = iwords[1:]

#pro_name = pro_marks.keys()
pro_name = os.getcwd().split('/')[-1];
pro_name = pro_name.split('_run')[0];
prmfile = sys.argv[1];
trjfile = sys.argv[2];
univ = Universe( prmfile, trjfile );

freq = 1;

# change this residues index in the same place as the follow index
# in ASH1L;
mark1, mark2, mark3, mark4, mark5 = int(pro_marks[pro_name][0]), int(pro_marks[pro_name][1])-2, \
       int(pro_marks[pro_name][2])+1, int(pro_marks[pro_name][3]), int(pro_marks[pro_name][6])-4;
print 'select conser, resi '+str(mark1)+'+'+str(mark2)+'+'+str(mark3)+'+'+str(mark4)+'+'+str(mark5);
#mark1, mark2, mark3, mark4, mark5 = 30, 33, 43, 51, 81;
# check mark4+1 and mark5+1, mark5-3;

P35_O = univ.selectAtoms( "resid " + str(mark1) +" and name O" );

N38_O = univ.selectAtoms( "resid " + str(mark2) +" and name O" );

Y48_OH = univ.selectAtoms( "resid " + str(mark3) + " and name OH");
Y48_HH = univ.selectAtoms( "resid " + str(mark3) + " and name HH");

L56_N = univ.selectAtoms( "resid " + str(mark4) + " and name N");
L56_H = univ.selectAtoms( "resid " + str(mark4) + " and name H");
L56_O = univ.selectAtoms( "resid " + str(mark4) + " and name O");

D57_CA = univ.selectAtoms( "resid " + str(mark4+1) + " and name CA");

V83_O = univ.selectAtoms( "resid " + str(mark5-3) + " and name O");

N86_ND2 = univ.selectAtoms( "resid " + str(mark5) +" and name ND2");
N86_CB = univ.selectAtoms( "resid " + str(mark5) +" and name CB");
N86_O = univ.selectAtoms( "resid " + str(mark5) + " and name O");

A87_N = univ.selectAtoms("resid " + str(mark5+1)+" and name N");

#protein = univ.selectAtoms( "protein")

#protein_writer = MDAnalysis.Writer("ASH1L_pro.pdb", len(protein));
wat1_writer = MDAnalysis.Writer("wat1.pdb");
wat2_writer = MDAnalysis.Writer("wat2.pdb");
wat3_writer = MDAnalysis.Writer("wat3.pdb");
wat4_writer = MDAnalysis.Writer("wat4.pdb");
wat5_writer = MDAnalysis.Writer("wat5.pdb");

wat1_id = file("wat1.txt",'w');
wat2_id = file("wat2.txt",'w');
wat3_id = file("wat3.txt",'w');
wat4_id = file("wat4.txt",'w');
wat5_id = file("wat5.txt",'w');

# mark1, mark2, mark3, mark4, mark5 = 35, 38, 48, 56, 86;
# P35_O, N38_O, Y48_OH, Y48_HH, L56_N, L56_H, L56_O, V83_O, N86_ND2, N86_O, 
# A87_N, Y91_OH, D57_CA
for ts in univ.trajectory:
    if ts.frame%freq == 0:
        #importWat = univ.selectAtoms( "byres (resname WAT and around 6.0 ( (resid 35 and name O) or (resid 38 and name O) or (resid 48 and name OH) or (resid 56 and name N) or (resid 56 and name O) or (resid 83 and name O) or (resid 86 and name ND2) ) )" );
	importWat = univ.selectAtoms( "byres (resname WAT and around 6.0 ( (resid " + str(mark1) + \
	            " and name O) or (resid " + str(mark2) + " and name O) or (resid " + str(mark3) \
		    + " and name OH) or (resid " + str(mark4) + " and name N) or (resid " + str(mark4) \
		    + " and name O) or (resid " + str(mark5-3) + " and name O) or (resid " + str(mark5) \
		    + " and name ND2) ) )" );
        #protein_writer.write(protein);
        wat1_id.write( str(ts.frame) + '   ' );
	wat2_id.write( str(ts.frame) + '   ' );
	wat3_id.write( str(ts.frame) + '   ' );
	wat4_id.write( str(ts.frame) + '   ' );
	wat5_id.write( str(ts.frame) + '   ' );
        wat1 = wat2 = wat3 = wat4 = wat5 = univ.selectAtoms("name ZYX");
        for each in importWat:
	    if each.name == 'O' and ( (lg.norm( Y48_OH[0].position-each.position ) < 4.0 and lg.norm( L56_N[0].position-each.position ) < 4.0 ) or \
	    ( lg.norm( Y48_OH[0].position-each.position ) < 4.0 and lg.norm( N86_ND2[0].position-each.position) < 4.0 ) or \
	    ( lg.norm( L56_N[0].position-each.position ) < 4.0 and lg.norm( N86_ND2[0].position-each.position) < 4.0 ) ):
            #if each.name == 'O' and lg.norm( Y48_OH[0].position-each.position ) < 4.0 and lg.norm( L56_N[0].position-each.position )\
	    #< 4.0 and lg.norm( N86_ND2[0].position-each.position) < 4.0:
                wat1 = each.residue;
            else:
	        if each.name == 'O' and lg.norm( V83_O[0].position-each.position ) < 3.5:
	            wat3 = each.residue;

	        elif each.name == 'O' and lg.norm( L56_O[0].position-each.position ) < 3.5 and lg.norm( D57_CA[0].position-each.position )\
		< 4.5 and lg.norm( L56_N[0].position-each.position) < 6.0:
		    if len(wat4) == 3 and lg.norm( D57_CA[0].position-each.position ) > lg.norm( D57_CA[0].position - wat4[0].position ):
		        pass;
		    else:
	                wat4 = each.residue;

	        elif each.name == 'O' and lg.norm( N38_O[0].position-each.position ) < 3.5 and lg.norm( P35_O[0].position-each.position )\
		< 3.5 and lg.norm( L56_O[0].position-each.position) < 7.0: 
		    if len(wat5) == 3 and lg.norm( L56_O[0].position-each.position) > lg.norm(L56_O[0].position - wat5[0].position ):
		        pass;
		    else:
	                wat5 = each.residue;

	        elif each.name == 'O' and ( lg.norm( Y48_OH[0].position-each.position ) < 3.5 or lg.norm( N86_CB[0].position-each.position ) < 5.0 ) \
		   and lg.norm( N86_O[0].position-each.position ) < 5.0 and lg.norm( A87_N[0].position-each.position ) < 5.0:
		    if len(wat2) == 3 and lg.norm( N86_O[0].position-each.position ) > lg.norm(N86_O[0].position - wat2[0].position ):
		        pass;
		    else:
		        wat2 = each.residue;
        if len(wat1) != 0:
            wat1_id.write( str( wat1.id ) + '  ' );
	    wat1_writer.write( wat1.atoms );
        if len(wat2) != 0:
	    wat2_id.write( str( wat2.id ) + '  ' );
	    wat2_writer.write( wat2.atoms );
        if len(wat3) != 0:
	    wat3_id.write( str( wat3.id ) + '  ' );
	    wat3_writer.write( wat3.atoms );
        if len(wat4) != 0:
	    wat4_id.write( str( wat4.id ) + '  ' );
	    wat4_writer.write( wat4.atoms );
	if len(wat5) != 0:
	    wat5_id.write( str( wat5.id ) + '  ' );
	    wat5_writer.write( wat5.atoms );
            #print "%.2f" %lg.norm( L56_O[0].position-wat5[0].position);
        wat1_id.write( '\n' );
	wat2_id.write( '\n' );
	wat3_id.write( '\n' );
	wat4_id.write( '\n' );
	wat5_id.write( '\n' );


#*X*X*X*X*X*X*X*X*X file close section
#protein_writer.close();
wat1_writer.close();
wat2_writer.close();
wat3_writer.close();
wat4_writer.close();
wat5_writer.close();
wat1_id.close();
wat2_id.close();
wat3_id.close();
wat4_id.close();
wat5_id.close();




