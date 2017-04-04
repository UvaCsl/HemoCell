#!/usr/bin/python
from os import environ
from sys import stderr

slots = 23

try:
  nodelist = environ['SLURM_NODELIST']
except KeyError:
  print >> stderr, "Environment variable SLURM_NODELIST not present, execute this within a job"
  exit(0)

final_nodelist = []
if "[" not in nodelist:
  final_nodelist.append(nodelist)
else:
  base = nodelist.split("[")[0]
  entries = nodelist.split('[')[1].rstrip("]").split(',')
  for entry in entries:
    if "-" in entry:
      first,second = entry.split('-')
      for i in range(int(first),int(second)+1):
          final_nodelist.append(base + str(i))
    else:
      final_nodelist.append(base + str(entry))

for node in final_nodelist:
  print node + " slots=" + str(slots)
