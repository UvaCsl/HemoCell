#!/usr/bin/python

import sys

npernode = 23
total = 384

print >> sys.stderr, "Using a total of " + str(((total-1)/npernode) + 2) + " nodes, 1 for the master process, " + str(((total-1)/npernode) + 1) + " for the other processes" 
print >> sys.stderr, "Total number of processes: " + str(total) + ". Processes per node : " + str(npernode)

print "rank 0=+n0 slot=1"

node = 1
processor = 1
for rank in range(1,total):
  print "rank " + str(rank) + "=+n" + str(node) + " slot=" + str(processor)
  processor += 1
  if processor > npernode:
    processor = 1
    node += 1

