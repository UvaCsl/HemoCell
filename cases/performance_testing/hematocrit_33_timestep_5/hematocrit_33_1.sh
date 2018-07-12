#@ job_type = parallel
#@ total_tasks = 1
#@ island_count = 1

#@ output = 1.out
#@ error  = 1.err
#@ wall_clock_limit = 00:30:00

#@ notification = always
#@ notify_user = s.a.alowayyed@uva.nl

#@ class = test

#@ energy_policy_tag = my_energy_tag 
#@ minimize_time_to_solution = yes

#@ queue

mpiexec -n 1 ../performance_testing ../configs_timestep_5/config_1.xml
