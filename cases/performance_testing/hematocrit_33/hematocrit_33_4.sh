#@ job_type = parallel
#@ total_tasks = 4
#@ island_count = 1

#@ output = 4.out
#@ error  = 4.err
#@ wall_clock_limit = 00:30:00

#@ notification = always
#@ notify_user = s.a.alowayyed@uva.nl

#@ class = test

#@ energy_policy_tag = my_energy_tag 
#@ minimize_time_to_solution = yes

#@ queue

mpiexec -n 4 ../performance_testing ../configs/config_4.xml
