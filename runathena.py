#!/usr/bin/python

from subprocess import call
import argparse
import os

ATHENA="/home/atripathi/Atmospheric-Athena/bin/athena"
CONFIGLOG="/home/atripathi/Atmospheric-Athena/config.log"
problemfile = "/home/atripathi/Atmospheric-Athena/src/problem.c"

# Inputs
parser = argparse.ArgumentParser(description="Athinput and output directory names")
parser.add_argument('-i', metavar='INFILE', type=str, default=None, help='athinput file'         )
parser.add_argument('-o', metavar='OUTDIR', type=str, default=None, help='output directory'      )
parser.add_argument('-mpi', nargs='?', metavar='PARALLEL', type=int, help='Use mpirun with specified number of procs' )
parser.add_argument('-r', metavar='RESTARTFILE', type=str, default=None, help='restart file'         )
parser.add_argument('-t', metavar='TLIM', type=float, default=None, help='time limit'         )

# Collect command line arguments
args = parser.parse_args()
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
if not args.r:
    in_file = args.i
    call("cp %s %s" % (in_file, out_dir_long), shell=True )
else:
    restart_file = args.r
    tlim = args.t
    call("cp %s %s" % (restart_file, out_dir_long), shell=True )
    
call("cp %s %s" % (CONFIGLOG, out_dir_long), shell=True )
call("cp %s %s" % (problemfile, out_dir_long), shell=True )


if args.mpi:
    mpinp = args.mpi
    if not args.r:
        call("/opt/openmpi/bin/mpirun -np %d %s -i %s -d %s" % (mpinp, ATHENA, in_file, out_dir), stdout=log_file, shell=True )
    else:
        call("/opt/openmpi/bin/mpirun -np %d %s -r %s -d %s time/tlim=%f" % (mpinp, ATHENA, restart_file, out_dir, tlim), stdout=log_file, shell=True )
else:
    if not args.r:
        call("%s -i %s -d %s" % (ATHENA, in_file, out_dir), stdout=log_file, shell=True )
    else:
        call("%s -r %s -d %s time/tlim=%f" % (ATHENA, restart_file, out_dir, tlim), stdout=log_file, shell=True )



