#!/usr/bin/python
# Author: Kevin Chan
# Data:   2014/1/25

import os;
import sys;

def pro_conser_marks( pro_name):
    '''
    open conserved_resi.txt file and get the resi index of 
    conserved residues in BRD protein
    '''
    ifile = file('/share/udata1/xiaoxiao/BRD_Ana_15.4.22/scripts/conserved_resi_v2.txt','r');
    
    iline = ifile.readline();
    iwords = iline.split();
    mark_names = iwords[1:13];

    pro_marks = {};
    for iline in ifile:
        iwords = iline.split();
        if len(iwords) != 13:
            continue;
        if len(iwords) == 13:
            pro_marks[ iwords[0] ] = iwords[1:]
    
    mark_index = {}
    for i in range( len(mark_names) ):
        mark_index[ mark_names[i] ] = pro_marks[ pro_name ][i];
    
    return mark_index;


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "Usage: ./pro_conser_marks.py pro_name";
        sys.exit();
    else:
        print pro_conser_marks( sys.argv[1] );


