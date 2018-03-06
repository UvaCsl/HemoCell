#PBS -lwalltime=3:00:00:00 
#PBS -lnodes=1:cores16:ppn=15 
NODES=1 
PPNODES=15 
CASE=/home/vazizi/hemocell/cases/pipeflow/
EXECUTABLE=pipeflow
CONFIG=config.xml
module load mpicopy 
cp -r "$CASE" "$TMPDIR"
cd "$TMPDIR"/"$CASE"/

#Create Rsync script
echo \#\!/bin/bash > rsync.sh
echo rsync -axr "$TMPDIR"/"$CASE"/"tmp" "$CASE" >> rsync.sh #Specify what you want synced to where
chmod +x rsync.sh

#initialize
mpicopy "$TMPDIR"/"$CASE"

#periodically sync
while true; do sleep 3600; mpirun -n "$NODES" -npernode 1 "$TMPDIR"/"$CASE"/rsync.sh; done &

#run case
mpirun -n `expr "$NODES" \* "$PPNODES"` -npernode "$PPNODES" "$EXECUTABLE"

#kill the periodic sync, wait for it to finish
kill %1 && sleep 60

#finalize by syncing one last time 
mpirun -n "$NODES" -npernode 1 "$TMPDIR"/"$CASE"/rsync.sh
~                                                         
