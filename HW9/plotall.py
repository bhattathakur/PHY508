#This program plots the energy (E), specific heat capacity (C_v), and Susceptibility (X) as the function of temperature (T)
#for three different values of lattice size for L=2,3,4

import numpy as np
import matplotlib.pyplot as plt
import time

ylabel=['','E',r'$C_{v}$',r'$\chi$']
savefile=['','EvsT','CvvsT','chivsT']
titles=['','Energy vs Temperature Plot',r'Specific heat $C_{v}$ vs Temperature Plot',r'Susceptibility $\chi$ vs Temperature Plot'] 
plt.rcParams['xtick.top']=True
plt.rcParams['ytick.right']=True

for j in range(1,4):
    for i in range(2,5):
        inputfile="data"+str(i)+".dat"                    #import data file T,E, C_v, X
        print("iputfile name:\t",inputfile)
        data=np.loadtxt(inputfile,skiprows=1,unpack=True) #skipped first row,unpacked the columns T,E,C_v,X
        xdata=data[0]    #T column
        ydata=data[j]    #y column
        plt.plot(xdata,ydata,'-',label='L= '+str(i))

    plt.xlabel('T')
    plt.ylabel(ylabel[j])
    plt.title(titles[j],pad=10)   #pad for offset from axis
    plt.legend()
    plt.grid(True)
    plt.minorticks_on()
    plt.grid(which='major',linestyle=':',linewidth='0.5',color='k')
    plt.grid(which='minor',linestyle=':',linewidth='0.25',color='k')
    plt.savefig(savefile[j]+".pdf")
    plt.show()
    time.sleep(1)
