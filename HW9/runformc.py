#This program run the Monte Carlo simulation for various input parameters and save the average values of <e> <e^2> <m> and <m^2> in a file

import os
import numpy as np
from scipy import stats

L=[4,4,4,3,3]
beta=[0.5,1,2,0.5,1.0]
Neql=[10**3]*5
Nmcs=[10**4]*5
Nbin=[10**3]*5
seed=[100]*5
latt='sqlatt_PBC'
datafile="outputdata.out"  #output from cpp run

for i in range(5):
    print("run\t",i+1)
    f=open("param.dat","w")  #input to the cpp progrm
    outfile='output'+str(i)+'.dat' #stores the data from 1 run
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
