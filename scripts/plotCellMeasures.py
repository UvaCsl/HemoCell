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

fname = '/home/lmount/devel/ficsion/tmp/plbCells.log'

i=0
pl.ion()
fig = pl.figure()
pl.show()
while True:
    i+= 1
    try:
        f = open(fname)
        c = StringIO(f.read())
        f.close()
        z = np.loadtxt(c, delimiter=';')
        t = z[:,0]
        V = z[:,1]
        S = z[:,2]
        TS = z[:,3]
        ED = z[:,4]
        A = z[:,5]
        EDlu = z[:,6]
        maxEDlu = z[:,7]
        xx = z[:,8]
        xy = z[:,9]
        xz = z[:,10]

        fig = pl.figure(1)
        pl.clf()
        pl.subplot(511)
        pl.plot(t, V)
        pl.ylabel('Volume %')
        pl.subplot(512)
        pl.plot(t, S)
        pl.ylabel('Surface %')
        pl.subplot(513)
        pl.plot(t, A)
        pl.ylabel('Angle degrees')
        pl.subplot(514)
        pl.plot(t, maxEDlu,'g')
        pl.ylabel('Max Edge Distance (LU)')
        pl.subplot(5,3,13)
        pl.plot(t, xx,'g')
        pl.title('x Position (LU)')
        pl.subplot(5,3,14)
        pl.plot(t, xy,'g')
        pl.title('y Position (LU)')
        pl.subplot(5,3,15)
        pl.plot(t, xz,'g')
        pl.title('z Position (LU)')
        pl.show()
        sleep(0.1)
        print i
    except:
        print 'Failed!'




