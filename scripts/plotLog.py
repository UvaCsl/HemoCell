#!/usr/bin/env python 
# encoding: utf-8
#====================================================================
# Name        :  
# Version    : 
# Author      :  Lampros Mountrakis (L.Mountrakis@uva.nl)
# Date        : 
# Description : Very simple script to plot the behaviour of several
#               variables while ficsion is running.
#
#====================================================================

import pylab as pl
import numpy as np
from StringIO import StringIO
from time import sleep
import sys


args = sys.argv[1:]
nArgs = len(args)-1
fname = args[0]


# fname = '/home/lmount/devel/ficsion/tmp/plbCells.log'

try:
    float(args[1])
    axes = map(int,args[1:])
    tAxes = None
except:
    tAxes = int(args[2])
    axes = map(int,args[3:])



i=0
pl.ion()
fig = pl.figure()
pl.show()
delim = ';'
while True:
    i+= 1
    try:
        f = open(fname)
        fields = f.readline().split(delim)
        c = StringIO(f.read())
        f.close()
        z = np.loadtxt(c, delimiter=delim)
        fig = pl.figure(1)
        pl.clf()
        for x in range(len(axes)):
            pl.subplot(len(axes),1,x+1)
            V = z[:,axes[x]]
            t = z[:,axes[tAxes]] if tAxes != None else np.arange(z.shape[0])
            pl.plot(t, V)
            pl.title(fields[axes[x]])
        pl.show()
        #pl.savefig('test.png');sys.exit()
        sleep(0.1)
        sys.stdout.write("\r%d iterations, %d lines" % (i,z.shape[0]) )
        sys.stdout.flush()
    except:
        print 'Failed!'




