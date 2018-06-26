#@ job_type = parallel
#@ total_tasks = 32
#@ island_count = 1

#@ output = 32.out
#@ error  = 32.err
#@ wall_clock_limit = 00:30:00

#@ notification = always
#@ notify_user = s.a.alowayyed@uva.nl

#@ class = test

#@ energy_policy_tag = my_energy_tag 
#@ minimize_time_to_solution = yes

#@ queue

mpiexec -n 32 ../performance_testing ../configs/config_32.xml
