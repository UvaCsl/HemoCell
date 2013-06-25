#!/usr/bin/env python 
# encoding: utf-8
#====================================================================
# Name        :  
# Version    : 
# Author      :  Lampros Mountrakis (L.Mountrakis@uva.nl)
# Date        : 
# Description : 
#
#
#====================================================================

import pylab as pl
import numpy as np
from StringIO import StringIO

sfile = open('tmp/stretchDeformation.log')
force = float(sfile.readline().split()[-1])
head = sfile.readline().split(';')
text = sfile.read().replace(';',' ')
c = StringIO(text)
data = np.loadtxt(c)
t = data[:,0]
D_A = data[:,2]*1e6
D_T = data[:,3]*1e6
sfile.close()


pl.plot(t, D_A,label='$D_A$')
pl.plot(t, D_T,label='$D_T$')
pl.xlabel('time [$s$]')
pl.ylabel('Diameter [$\mu{}m$]')
pl.title("Force = %03d pN"%(force*1e12))
pl.legend(loc='best')
#pl.show()
pl.savefig("./force.%03d.png"%(force*1e12))


print "%f %f %f"%(force*1e12, D_A[-1], D_T[-1],)


