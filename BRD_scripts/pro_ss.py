#!/usr/bin/python

import os;
import sys;
from subprocess import Popen;

# Open the dssp analysis output file, get the beginning and ending resi of 
# Helix Z, A, B, and C;

def pro_ss( dssplog='dssp_ana.out'):
    '''
    return the beginning, ending resi number of loop ZA and BC
    return the total resi length of protein
    '''

    pro_st=[];
    
    #ifile = file('/share/udata1/ck/BRDApo_GMX_99SBildn/'+ pro_name +'/dssp_ana.out','r');
    ifile = file('dssp_ana.out', 'r');
    ilines = ifile.readlines();
    ifile.close();
    iread=False;
    for iline in ilines:
        iwords = iline.split();
        if iread == True and len(iwords)>1 and len(iline)>16:
            pro_st.append([iwords[1], iline[16]]);
        if len(iwords) > 9 and iwords[1]=='RESIDUE' and iwords[2]=='AA' and iwords[3]=='STRUCTURE':
            iread = True;
    
    HZ = [];
    HA = [];
    HB = [];
    HC = [];
    HelixL = 7;
    
    jread1 = True;
    jread2 = True;
    jread3 = True;
    jread4 = True;
    
    #print pro_st;
    for each in pro_st:
        if jread1 == True:
            if each[1] == 'H':
                HZ.append(each);
            if each[1] != 'H':
                if len(HZ) > HelixL:
                    jread1 = False;
    	        else:
    	            HZ = [];

        if len(HZ) > HelixL and jread2 == True and jread1 == False:
            if each[1] == 'H':
                HA.append(each);
            if each[1] != 'H':
    	        if len(HA) > HelixL:
    	            jread2 = False;
    	        else:
    	            HA = [];

        if len(HA) > HelixL and jread2 == False and  jread3 == True:
            if each[1] == 'H':
    	        HB.append(each);
    	if each[1] != 'H':
    	    if len(HB) > HelixL:
    	        jread3 = False;
    	    else:
    	        HB = [];
        
        if len(HB) > HelixL and jread3 == False and jread4 == True:
            if each[1] == 'H':
    	        HC.append(each);
    	if each[1] != 'H':
    	    if len(HC) > HelixL:
    	        jread4 = False;
    	    else:
    	        HC = [];
        if jread4 == False:
            break;
    
    # ZAB is the ZA loop Beginning resi;
    # ZAE is the ZA loop Ending resi;
    # BCB is the BC loop Beginning resi;
    # BCE is the BC loop Ending resi;
    # END is the End resi of the protein;
    HZB = HZ[0][0];
    HZE = HZ[-1][0];
    HAB = HA[0][0];
    HAE = HA[-1][0];
    HBB = HB[0][0];
    HBE = HB[-1][0];
    print HC
    HCB = HC[0][0];
    HCE = HC[-1][0];
    END = pro_st[-1][0];

    print ' HZB  HZE  HAB  HAE  HBB  HBE  HCB  HCE  END'; 
    print '  '+HZB+'   '+HZE+'   '+HAB+'   '+HAE+'   '+HBB+'   '+HBE+'   '+HCB+'   '+HCE+'   '+END;
    return (HZE, HAB, HAE, HBB,  HBE, HCB, END);


if __name__=='__main__':
    if len(sys.argv) > 1: 
        print pro_ss( sys.argv[1]);
    else:
        print pro_ss();
        #pro_ss();









