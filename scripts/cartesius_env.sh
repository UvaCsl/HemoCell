#!/bin/bash
[[ $_ != $0 ]] || (echo "Script is a subshell, this wont work, source it instead!"; exit)
module rm mpi fortran c
module load mpi/openmpi/1.10.2
module load gcc/5.2.0
module list
echo "the TMP dir on cartesius seems to be FUBAR, just setting it to ~/tmp for now ..."
mkdir -p ~/tmp
export TMP=~/tmp
export TMPDIR=~/tmp
env |grep TMP
