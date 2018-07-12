#@ job_type = parallel
#@ total_tasks = 128
#@ island_count = 1

#@ output = 128.out
#@ error  = 128.err
#@ wall_clock_limit = 00:30:00

#@ notification = always
#@ notify_user = s.a.alowayyed@uva.nl

#@ class = test

#@ energy_policy_tag = my_energy_tag 
#@ minimize_time_to_solution = yes

#@ queue

mpiexec -n 128 ../performance_testing ../configs_timestep_5/config_128.xml
