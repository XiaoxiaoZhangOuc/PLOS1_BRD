
CutOff = '0.12'

Letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

from sys import argv, exit
from time import sleep
from numpy import array
from subprocess import Popen, PIPE
from MDAnalysis import Universe, Writer

if len(argv) < 3 :
    print 'Usage: python  g_cluster_outPDB.py  protein ndx '
    exit()

fileBase = argv[1]
NDX = argv[2]

T_begin = raw_input( 'Please enter the first time frame in ns (eg. D700 or 700): ' )
T_end   = raw_input( 'Please enter the  last time frame in ns (eg. 750): ' )

filename = fileBase + '.' + T_begin + '_' + T_end + '.clus'

logfile = filename + 'ZA'+'.log'
pdbfile = filename + '.pdb'

try :
    ifile = file( logfile, 'r' )
except :
    print 'No corresponding .clus.log file found.'
    if raw_input( 'If you want to run clustering using Gromacs, enter y: ') in ('y','Y','yes') :
        dt = int( ( float(T_end) - float(T_begin.strip(Letters)) ) / 2.0 )
        dt = raw_input( 'Please enter dt (~' + str(dt) + ' for 2000 structures): ' )
        p = Popen( [ 'g_cluster_mpi', '-f', str(fileBase)+'.xtc', '-s', str(fileBase) +'.pdb', '-n',  NDX, '-g', logfile, '-cl', pdbfile,
                     '-b', T_begin.strip(Letters)+'000', '-e', T_end+'000', '-dt', dt, '-cutoff', CutOff ],
                   stdin=PIPE ) # stdout=open('/dev/null', 'a'), stderr=open('/dev/null', 'a') )
        p.communicate( '3\n1\n' )
        p.wait()
        sleep(3)
        ifile = file( logfile, 'r' )
    else :
        exit()

N_str = [ ] # number of structures for each cluster
Index = { }

for line in ifile :
    if len(line) > 21 and line[4] == '|' \
    and line[6] != '#' and line[8] != ' ' :
        n = int( line[:4] )
        nst = int( line[5:10] )
        if n == len(N_str) + 1 :
            N_str.append( nst )
            try :
                Index[ nst ].append( n )
            except :
                Index[ nst ] = [ n ]
        else :
            print "Error: ", n, nst, N_str
ifile.close()

print 'number of clusters:', len(N_str)

N_str.sort()
N_str.reverse()

N_write = int( raw_input( 'How many clusters you want to look at (eg. 8): ' ) ) 

universe = Universe( pdbfile, pdbfile )

writer = Writer( filename + '_.pdb', universe.atoms.numberOfAtoms() )

k = 1
for each in N_str[ : N_write ] :
    if k > N_write :
        break
    percent = float(each) * 100 / sum(N_str)
    for i in Index[ each ] :
        print 'Writing cluster #' + str(i).ljust(3),
        universe.trajectory[i-1]
        protein = universe.atoms
        writer.write( protein )
        print '%3i %5.0f%% %6.2f%%' %( k, percent, percent )
        k += 1

writer.close()

