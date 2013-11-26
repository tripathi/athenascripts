#!/usr/bin/python

from subprocess import call
import argparse
import os

ATHENA="/Users/anjalitripathi/Atmospheric-Athena/bin/athena"
CONFIGLOG="/Users/anjalitripathi/Atmospheric-Athena/config.log"

# Inputs
parser = argparse.ArgumentParser(description="Athinput and output directory names")
parser.add_argument('-i', metavar='INFILE', type=str, default=None, help='athinput file'         )
parser.add_argument('-o', metavar='OUTDIR', type=str, default=None, help='output directory'      )
parser.add_argument('-mpi', nargs='?', metavar='PARALLEL', type=int, help='Use mpirun with specified number of procs' )

# Collect command line arguments
args = parser.parse_args()
in_file = args.i
out_dir = args.o
out_dir = out_dir.rstrip('/')
out_dir_long = out_dir+"/"

#Make directory if it doesn't already exist
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
    
# Create log-file to store output
LOG_NAME="runathena.log"
log_file = open(out_dir + "/" + LOG_NAME, "w+")

#Run athena and copy input file to the output directory    
call("cp %s %s" % (in_file, out_dir_long), shell=True )
call("cp %s %s" % (CONFIGLOG, out_dir_long), shell=True )
if args.mpi:
    mpinp = args.mpi
    call("mpirun -np %d %s -i %s -d %s" % (mpinp, ATHENA, in_file, out_dir), stdout=log_file, shell=True )
else:
    call("%s -i %s -d %s" % (ATHENA, in_file, out_dir), stdout=log_file, shell=True )




