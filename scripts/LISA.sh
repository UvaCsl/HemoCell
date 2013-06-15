#PBS -S /bin/bash
#PBS -lwalltime=00:05:00 -lnodes=1:cores8:ppn=8

#PBS -lwalltime=72:00:00 -lnodes=1:cores8:ppn=8
#PBS -lwalltime=120:00:00 -lnodes=1:cores8:ppn=8

# ========================= LISA Batch Submit Script ========================= #
#    The purpose of this script is to submit multiple jobs to the cluster of
# LISA and have up-to-date results to the User Directory (UD).
# 
#    On LISA, jobs are running on Worker Nodes (WN) and the User Directory 
# "/home/lmountra/" is accessible from the WN. The Input and the Output of the 
# Job is performed in the Scratch system of the WN and is synchronized with the 
# UD every 1hour ($SLEEPTIME). Upon synchronization it takes care that the files
# to be synced are not opened by another process.(functions gzip_rsync_sleep, 
# safe_gzip, IsClosed).


# //// For multiple cases with different initial conditions submit with:
#####  for i in {01..15}; do qsub -v JOBID=$i LISA.sh ; done ######
#                 ////


# === Variables that might need change ====
export PATH=$PATH:/home/lmountra/bin/
export RANDID=${JOBID-${RANDOM}} # Identification of the job, randomized
export SUBJECT="RBC benchmark $RANDID" # Subject of the email sent to USER
export JobName=RBCBenchmark$RANDID # Used for directory making purposed and identifying the job
export ParentDir=~/RBCBenchmark/ # PATH where initial conditions (cases), results and snapshots are to be found and copied

export EXE=/home/lmountra/ficsion/ficsion # Path to the ficsion executable
export CREATECONFIG=/home/lmountra/ficsion/scripts/conficsion.py # Usually unnecessary, but check the rest of the script
export INITIALCONFIG=/home/lmountra/ficsion/config.xml
export LNS=/home/lmountra/ficsion/lib/

# ===== 
ulimit -v unlimited
ulimit -s unlimited
export OMP_NUM_THREADS=1
export PROCS=8
SLEEPTIME=60 # Wallclock time interval between each RSYNC (1h)


export ResultDir=${ParentDir}/results/${JobName}
export ScratchDir=${TMPDIR}/${JobName}

mkdir -p ${ResultDir} ${ScratchDir}
cd ${ScratchDir}

echo ${ScratchDir}

# ======== Sent initial mail that the Job has started ==========
(
echo "Job ${JobName} with JobID ${PBS_JOBID} started at `date`";
echo 
echo "ssh ${HOSTNAME} ";
echo "cd ${ScratchDir} ";
) | mail $USER -s "$SUBJECT"

# ======== Useful functions ==========
IsClosed() {
    lsof -c $(basename ${EXE}) -c tar | grep ${1} &> /dev/null
    echo $?
}

safe_gzip () {
        if [ $(IsClosed ${1}) -eq 1 ]; then
            gzip -f ${1}
        fi
}

wait_until_nproc_drops_below () {
    NPROC=$1
    echo "JOBS" `jobs | wc -l`
    while [ `jobs | wc -l` -ge $NPROC ]; do
        sleep 120
    done
}

gzip_rsync_sleep () {
    NPROC=$1
    cd ${ScratchDir}
    echo "STARTED WITH " `jobs | wc -l`
    while [ `jobs | wc -l` -ge $NPROC ]; do
        for fl in ${ScratchDir}/*/*/*.*; do
            safe_gzip ${fl}
        done
        rsync -zvr ${ScratchDir}/ ${ResultDir}/
        sleep ${SLEEPTIME}
    done
}

stretchForces="--stretchForce 0 15 19 30 38 47 67 88 108 130 150 172 192"
# ======== Create Config Files ==========
(
mkdir ${ScratchDir}/configurations;
cd ${ScratchDir}/configurations;
cp ${INITIALCONFIG} config_test.xml
${CREATECONFIG} ${stretchForces} --flowType 3 --shearRate 0
rm config_test.xml
);

# ========= Main Loop ================ # 
# Change According to your needs,      #
#              but try to fit 8 cases! #
# ==================================== #
for configFile in ${ScratchDir}/configurations/*xml; do
    (
        config=$(basename ${configFile} .xml);
        jobid=${config}
        export jobdir=${ScratchDir}/${jobid}/
        mkdir -p ${jobdir}/tmp/
        cd ${jobdir};
        ln -s ${LNS}
        cp ${configFile} ${config}.xml
        ${EXE} ${config}.xml &> ${jobid}.log
    ) &
    gzip_rsync_sleep $PROCS  &> ${ScratchDir}/gzipped.log 
done

# Leave unchanged if you want to keep the sync
wait
rsync -zvr ${ScratchDir}/ ${ResultDir}/ &>>  ${ScratchDir}/gzipped.log 
( echo "Job $PBS_JOBID ended at `date`.";  ) | mail $USER -s "$SUBJECT"




