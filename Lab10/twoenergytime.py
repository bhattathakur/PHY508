# This program prints the energy of the plot vs Monte Carlo Time

import numpy as np
import matplotlib.pyplot as plt

#read the datfiles
datahigh=np.loadtxt('datahigh.dat')   #data for small T nmcs=1
datalow=np.loadtxt('datalow.dat')     #data for big T nmcs=1
highTnmcs=np.loadtxt('highnmcs.dat',unpack=True)[0] #data with high T and 10 nmcs
lowTnmcs=np.loadtxt('lownmcs.dat',unpack=True)[0]   #data with low T and 10 nmc
hi=datahigh[:,0]
li=datalow[:,0]
#labels
label=[r'T=$\infty$,nmcs=1',r'T=$0$,nmcs=1']
labels=[r'T=$\infty$, nmcs=100',r'T=$0$, nmcs=100']
#Monte Carlotime
x=np.arange(len(hi))

plt.figure()
plt.subplot(211)
plt.plot(x,hi,label=label[0])
plt.plot(x,li,label=label[1])
plt.title("Energy as the function of MonteCarlo Time")
plt.legend()
plt.xlabel("MonteCarlo Time")
plt.ylabel("Energy")

plt.subplot(212)
plt.plot(x,highTnmcs,label=labels[0])
plt.plot(x,lowTnmcs,label=labels[1])
plt.xlabel("MonteCarlo Time")
plt.ylabel("Energy")
plt.legend()
plt.tight_layout()
plt.savefig('evsmctime.pdf')
plt.savefig('evsmctime.png')
plt.show()

