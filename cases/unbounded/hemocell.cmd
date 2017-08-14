#@ job_type = parallel
#@ node = 1
#@ total_tasks = 1
#@ island_count = 1

#@ output = job.out
#@ error = job.err
#@ wall_clock_limit = 00:30:00

#@ notification = error
#@ notify_user = s.a.alowayyed@uva.nl

#@ class = test

#@ energy_policy_tag = my_energy_tag 
#@ minimize_time_to_solution = yes

#@ queue

poe ./unbounded config.xml