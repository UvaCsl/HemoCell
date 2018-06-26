#@ job_type = parallel
#@ total_tasks = 2
#@ island_count = 1

#@ output = 2.out
#@ error  = 2.err
#@ wall_clock_limit = 00:30:00

#@ notification = always
#@ notify_user = s.a.alowayyed@uva.nl

#@ class = test

#@ energy_policy_tag = my_energy_tag 
#@ minimize_time_to_solution = yes

#@ queue

mpiexec -n 2 ../performance_testing ../configs/config_2.xml
