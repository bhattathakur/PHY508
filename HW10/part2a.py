#This program finds the Auto-correlation function as the function of Monte-Carlo time and also gets the fit values and stores i#in the file.

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os


#Parameters
Nlin=4              #linear size of lattice
beta=0.440687       #1/T at critical point
Neql=100            #eql sweeps
Nmcs=1              #sweeps in a bin
Nbin=10000          #number of bin
SEED=10             #for MT
latt="sqlatt_PBC"   #lattice kind

#mc="wolff"
mc="met"
param="param.dat"             #parameter file
wolfdata="wolfoutput.out"     #output after running ./wolff.out
#metdata="metoutput.out"      #output after metropolis alogrithm
run=0
total_run=200
for L in [5,10,20,40,80]:
    print("run\t",run+1);run+=1;

    #write the parameters value in the file
    f=open(param,'w')
    f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(L,beta,Neql,Nmcs,Nbin,SEED,latt))
    f.close()

    #run the output of ising_wolf.cpp
    os.system('./met.out > met.log')

    #load the resultant output of the wolff algorithm
    datawolff=np.loadtxt(wolfdata,unpack=True)
    
    #magnetization and magnetization square
    m=datawolff[2]
    msq=datawolff[3]
    #print m and msq
    #print("m=\t{:.5f}\tmsq=\t{:.5f}\n".format(m,msq))

    #Atau function
    def Atau(tau):
        msqavg=np.sum(msq)/Nbin
        if(tau!=0):
            magshift=np.roll(m,-tau)
        else:
            magshift=m
        mtauavg=np.dot(m,magshift)/Nbin
        return abs(mtauavg/msqavg)

    #value of x
    x=np.arange(total_run)
    Alist=[Atau(n) for n in range(total_run)]

    #define the fit function
    def func(x,a):
        return np.exp(-x/a)

    #get fit values
    popt,pconv=curve_fit(func,x,Alist)    #fit function, x-value, y-value
    fit=func(x,popt[0])                   #use fit values to the function
    th=popt[0]
    #theta value
    theta=[th for i in range(total_run)]
    y=np.column_stack((x,Alist,fit,theta)) #np.column_stack(tuple)
    np.savetxt(mc+"data"+str(L)+".dat",y)  #save in a file
