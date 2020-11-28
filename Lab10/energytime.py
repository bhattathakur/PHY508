# This program prints the energy of the plot vs Monte Carlo Time

import numpy as np
import matplotlib.pyplot as plt

#datahigh=np.loadtxt('datahighT.dat')   #data for small T
datahig=np.loadtxt('datahigh.dat',unpack=True)
datalow=np.loadtxt('datalow.dat',unpack=True)
#plt.plot(data[0])
#datalow=np.loadtxt('datalowT.dat') #data for big T
#hi=datahigh[:,0]
#li=datalow[:,0]
#label=[r'T=$\infty$',r'T=$0$']
##print("len e",len(e))
##x0=np.arange(len(e))
#x=np.arange(len(hi))
#plt.plot(x,hi,label=label[0])
#plt.plot(x,li,label=label[1])
#plt.xlabel("MonteCarlo Time")
#plt.ylabel("Energy")
#plt.legend()
plt.show()

