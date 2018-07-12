#@ job_type = parallel
#@ total_tasks = 1024
#@ island_count = 1

#@ output = 1024.out
#@ error  = 1024.err
#@ wall_clock_limit = 00:30:00

#@ notification = always
#@ notify_user = s.a.alowayyed@uva.nl

#@ class = test

#@ energy_policy_tag = my_energy_tag 
#@ minimize_time_to_solution = yes

#@ queue

mpiexec -n 1024 ../performance_testing ../configs_timestep_5/config_1024.xml
