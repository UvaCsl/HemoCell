#PBS -S /bin/bash
#PBS -lwalltime=00:05:00 -lnodes=1:cores8:ppn=8

#PBS -lwalltime=120:00:00 -lnodes=1:cores8:ppn=8



#PBS -lwalltime=72:00:00 -lnodes=1:cores8:ppn=8

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
#####  for i in 75 150 200 1000; do qsub -v KBEND=$i LISA.sh ; done ######
#                 ////


# === Variables that might need change ====
export PATH=$PATH:/home/lmountra/bin/
export RANDID=${SHEARRATE-${RANDOM}} # Identification of the job, randomized
export SUBJECT="LISA RUN: DSFD -- RBC Shear Flow $RANDID" # Subject of the email sent to USER
export JobName=DSFD_SHEARFLOW # Used for directory making purposed and identifying the job
export ParentDir=/home/lmountra/DSFD/ShearFlow/ # PATH where initial conditions (cases), results and snapshots are to be found and copied

export EXE=/home/lmountra/ficsion/ficsion # Path to the ficsion executable
export CREATECONFIG=/home/lmountra/ficsion/scripts/conficsion.py # Usually unnecessary, but check the rest of the script
export LNS=/home/lmountra/ficsion/lib/

echo $ParentDir
# ===== 
ulimit -v unlimited
ulimit -s unlimited
export OMP_NUM_THREADS=1
export PROCS=8
SLEEPTIME=3600 # Wallclock time interval between each RSYNC (1h)


export ResultDir=${ParentDir}
export ScratchDir=${TMPDIR-/scratch/}/${JobName}

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

rsync_sleep () {
    NPROC=$1
    cd ${ScratchDir}
    echo "STARTED WITH " `jobs | wc -l`
    SYNCINTERVAL=0
    if [ `jobs | wc -l` -ge $NPROC ]; then
        rsync -zvr ${ScratchDir}/ ${ResultDir}/;
    fi
    while [ `jobs | wc -l` -ge $NPROC ]; do
        if [ ${SYNCINTERVAL} -ge ${SLEEPTIME} ]; then
            killall -9 rsync;
            rsync -zvr ${ScratchDir}/ ${ResultDir}/;
            SYNCINTERVAL=0;
        fi
# Sleep 1 min
        sleep 60; 
        SYNCINTERVAL=$[${SYNCINTERVAL} + 60];
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

# ======== Create Config Files ==========
(
mkdir -p ${ScratchDir}/configurations;
cd ${ScratchDir}/configurations;
cp /home/lmountra/DSFD/TTT.xml config_test.xml
${CREATECONFIG} --rbcModel 0 1 --shearRate ${SHEARRATE} --caseId TTT --etaM 0.0 0.001e-6 
rm config_test.xml
cp /home/lmountra/DSFD/Yao.xml config_test.xml
${CREATECONFIG} --rbcModel 0 1 --shearRate ${SHEARRATE} --caseId Yao --etaM 0.0 0.001e-6 
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
        mkdir -p ${jobdir}/{tmp,lib}/
        cd ${jobdir};
        cp -r ${LNS}/*  ${jobdir}/lib/
        cp ${configFile} ${config}.xml
        echo ${EXE} ${config}.xml
        ${EXE} ${config}.xml &> ${jobid}.log
    ) &
    rsync_sleep $PROCS  &>> ${ScratchDir}/gzipped.log
done

# Leave unchanged if you want to keep the sync
rsync_sleep 1  &>> ${ScratchDir}/gzipped.log
wait
rsync -zvr ${ScratchDir}/ ${ResultDir}/ &>>  ${ScratchDir}/gzipped.log 
( echo "Job $PBS_JOBID ended at `date`.";  ) | mail $USER -s "$SUBJECT"


# Post Processing can be done with:
##  for d in rbcM*; do (cd $d; echo We are on $d; for i in *; do (cd $i; echo $i; ~/ficsion/scripts/stretchDeformation.py >> ../stretchLog.dat;); done); done

