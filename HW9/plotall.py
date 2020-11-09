#This program plots the specific heat capacity (C_v), and Susceptibility (X) as the function of temperature (T) with ising model and 
#metropolis alogrithm for three different values of lattice size for L=4 (ising model), L=4,8 (Metropolis algorithm)

import numpy as np
import matplotlib.pyplot as plt
import time

ylabel=['',r'$C_{v}$',r'$\chi$']
savefile=['','CvvsT3','chivsT3']
titles=['',r'Specific heat $C_{v}$ vs Temperature Plot',r'Susceptibility $\chi$ vs Temperature Plot'] 
plt.rcParams['xtick.top']=True
plt.rcParams['ytick.right']=True
#datafiles l4enumdata4.dat,mcl4data.dat,mcl8data.dat
inputfiles=["l4enumdata4.dat","mcl4data.dat","mcl8data.dat"]
labels=["L=4 (enum)","L=4 (mc)","L=8 (mc)"]


for j in range(1,3):
    for i in range(len(inputfiles)):
        #inputfile=inputfiles[i]+".dat"                    #import data file T,C_v, X
        #inputfile=inputfiles[i]
        #print("iputfile name:\t",inputfile)
        data=np.loadtxt(inputfiles[i],skiprows=1,unpack=True) #skipped first row,unpacked the columns T,C_v,X
        xdata=data[0]    #T column
        ydata=data[j]    #y column
        plt.plot(xdata,ydata,'-',label=labels[i])

    plt.xlabel('T')
    plt.ylabel(ylabel[j])
    plt.title(titles[j],pad=10)   #pad for offset from axis
    plt.legend()
    plt.grid(True)
    plt.minorticks_on()
    plt.grid(which='major',linestyle=':',linewidth='0.5',color='k')
    plt.grid(which='minor',linestyle=':',linewidth='0.25',color='k')
    plt.savefig(savefile[j]+".pdf")
    plt.savefig(savefile[j]+".png")
    plt.show()
    time.sleep(1)
