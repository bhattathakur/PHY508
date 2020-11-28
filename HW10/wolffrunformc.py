#This program run the Monte Carlo simulation based on Wolff single cluster algorithm for various input parameters and save the average values of <e> <e^2> <m> ,<m^2>, <average size> in a file

import os
import numpy as np
from scipy import stats

#L=[4,4,4,3,3]
L=[8,16]
#beta=[0.5,1,2,0.5,1.0]
beta=[0.5,1]
#Neql=[10**3]*5
Neql=[10**3]*2
#Nmcs=[10**4]*5
Nmcs=[10**2]*2
#Nbin=[10**3]*5
Nbin=[10**2]*2
seed=[100]*2
latt='sqlatt_PBC'
datafile="outputdata.out"  #output from cpp run

for i in range(len(L)):
    print("run\t",i+1)
    f=open("param.dat","w")  #input to the cpp progrm
    outfile='wolfoutput'+str(i)+'.dat' #stores the data from 1 run
    if i==0:print("{}\t{}\t{}\t{}\t{}\t{}".format("L","beta","Neql","Nmcs","Nbin","SEED","latt_"))
    print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(L[i],beta[i],Neql[i],Nmcs[i],Nbin[i],seed[i],latt))
    #write in param.dat
    f.write("{}\t{}\t{}\t{}\t{}\t{}\t{:}".format(L[i],beta[i],Neql[i],Nmcs[i],Nbin[i],seed[i],latt))
    f.close()
    #use parama.dat as input
    os.system('./a.out')
    #copy the datainto different name
    command='cp '+datafile+' '+outfile
    os.system(command)
#print("Using the data from Wolff single cluster algorithm\n")
#print("Run calculatemc.py function\n")
#os.system('python calculatemc.py')

