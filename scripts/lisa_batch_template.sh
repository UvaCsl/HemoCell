#PBS -lwalltime=5:00 
#PBS -lnodes=1:cores16:ppn=15 
NODES=1 
PPNODES=15 
CASE=pipeflow
CASE_FOLDER=/home/vazizi/hemocell/cases/"$CASE"
EXECUTABLE=pipeflow
CONFIG=config.xml

#Load Modules
module load mpicopy 
module load openmpi/gnu
cp -r "$CASE_FOLDER" "$TMPDIR"
cd "$TMPDIR"/"$CASE"/

#Create Rsync script
echo \#\!/bin/bash > rsync.sh
echo rsync -axr "$TMPDIR"/"$CASE"/"tmp" "$CASE_FOLDER" >> rsync.sh #Specify what you want synced to where
chmod +x rsync.sh

#initialize
mpicopy "$TMPDIR"/"$CASE"

#periodically sync
while true; do sleep 1; mpirun -n "$NODES" -npernode 1 "$TMPDIR"/"$CASE"/rsync.sh; done &

#run case
mpirun -n `expr "$NODES" \* "$PPNODES"` -npernode "$PPNODES" "$EXECUTABLE" "$CONFIG"

#kill the periodic sync, wait for it to finish
kill %1 && sleep 60

#finalize by syncing one last time 
mpirun -n "$NODES" -npernode 1 "$TMPDIR"/"$CASE"/rsync.sh
~                                                         
