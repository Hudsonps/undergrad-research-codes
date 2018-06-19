# -*- coding: utf-8 -*-
from pylab import *
from matplotlib import *
import numpy as np
from math import *

Wname = "4.00"

imp = file("/home/hudson/Desktop/dotinterage/dotimpuro/mag_susceptimp%s.plt"%Wname)

puro = file("/home/hudson/Desktop/dotinterage/dotpuro/mag_suscept%s.plt"%Wname)

temp_ = []
ximp_ = []
xpuro_ = []
xdiff_= []

for line in imp.readlines():
  sline = line.split()
  temp_.append(float(sline[0]))
  ximp_.append(float(sline[1]))
  
  
for iline,line in enumerate(puro.readlines()):
  sline = line.split()
  xdiff_.append(-float(sline[1]) + ximp_[iline])
  xpuro_.append(float(sline[1]))
  

Fig = pyplot.figure()
fig = Fig.add_subplot(111)


fig.semilogx(temp_, ximp_, color = 'red', label='Impureza', linewidth = 2.)
fig.semilogx(temp_, xpuro_, label='Sem impureza', color = 'blue', linewidth = 2.0)

fig.legend()

xlabel("$T$")
ylabel("$k_B T \chi$")

fig.text(2e-3, 0.35,'$\epsilon_d = -0.1$\n$U=-2\epsilon_d$\n$V=0.05$\n$W=4$\n$N=3$')



pylab.savefig("mag_susceptWinf.eps", format = "eps")
