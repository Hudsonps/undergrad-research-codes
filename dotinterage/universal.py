# -*- coding: utf-8 -*-
from pylab import *
from matplotlib import *
import numpy as np
from math import *

Wname = "-0.00"

TK = 1.65e-5

imp = file("/home/hudson/Desktop/dotinterage/dotimpuro/mag_susceptimp%s.plt"%Wname)

puro = file("/home/hudson/Desktop/dotinterage/dotpuro/mag_suscept%s.plt"%Wname)

univ = file("/home/hudson/Desktop/dotinterage/univ.s")

temp_ = []
ximp_ = []
xpuro_ = []
xdiff_= []

tempuniv_ = []  ##Como as temperaturas no arquivo da curva universal são diferentes dos outros, uso um vetor especial para este caso.
xuniv_ = []

for line in imp.readlines():
  sline = line.split()
  temp_.append(float(sline[0]))
  ximp_.append(float(sline[1]))
  
  
for iline,line in enumerate(puro.readlines()):
  sline = line.split()
  xdiff_.append(-float(sline[1]) + ximp_[iline])
  xpuro_.append(float(sline[1]))
  
for iline,line in enumerate(univ.readlines()):
  sline = line.split()
  tempuniv_.append(float(sline[0])/float(TK))
  xuniv_.append(float(sline[1]))  
  

Fig = pyplot.figure()
fig = Fig.add_subplot(111)


fig.semilogx(temp_, xdiff_, color = 'red', label='$k_B T (\chi_{{i}} - \chi_{{p}})$', linewidth = 2.)
fig.semilogx(tempuniv_, xuniv_, label='Curva universal', color = 'blue', linewidth = 2.0)

fig.legend()

xlabel("$T$")
ylabel("$k_B T \chi$")

fig.text(2e-3, 0.15,'$\epsilon_d = -0.1$\n$U=-2\epsilon_d$\n$V=0.05$\n$W=0$\n$N=3$')

pylab.savefig("universal.eps", format = "eps")
