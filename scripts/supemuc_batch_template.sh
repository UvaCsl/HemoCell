#@ job_type = parallel
#@ node = 25
#@ total_tasks = 400
#@ island_count = 1

#@ output = job_400.out
#@ error  = job_400.err
#@ wall_clock_limit = 00:20:00

#@ notification = always
#@ notify_user = s.a.alowayyed@uva.nl

#@ class = test

#@ energy_policy_tag = my_energy_tag 
#@ minimize_time_to_solution = yes

#@ queue

mkdir $SCRATCH/hemocell
cp pipeflow *.pos *.stl *.xml $SCRATCH/hemocell 
echo "copy done .. "

cd $SCRATCH/hemocell

mpiexec -n 400 pipeflow $SCRATCH/hemocell/tmp/checkpoint.xml
