#@ job_type = parallel
#@ total_tasks = 256
#@ island_count = 1

#@ output = 256.out
#@ error  = 256.err
#@ wall_clock_limit = 00:30:00

#@ notification = always
#@ notify_user = s.a.alowayyed@uva.nl

#@ class = test

#@ energy_policy_tag = my_energy_tag 
#@ minimize_time_to_solution = yes

#@ queue

mpiexec -n 256 ../performance_testing ../configs/config_256.xml
