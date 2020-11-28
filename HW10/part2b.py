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
param="param.dat"              #parameter file
#wolfdata="wolfoutput.out"     #output after running ./wolff.out
metdata="metoutput.out"        #output after metropolis alogrithm
run=0
#different variable
L=[]
zlist=[]
Atau_list=[]
Atau_sum=[]
total_run=300
Nlinrange=np.arange(4,65,4)
for Nlin in Nlinrange:#np.arange(4,65,4):
    L.append(Nlin)
    print("run\t",run+1);run+=1;

    #write the parameters value in the file
    f=open(param,'w')
    f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(Nlin,beta,Neql,Nmcs,Nbin,SEED,latt))
    f.close()

    #run the output of ising_wolf.cpp
    os.system('./met.out > met.log')

    #load the resultant output of the wolff algorithm
    datawolff=np.loadtxt(metdata,unpack=True)
    
    #magnetization and magnetization square
    m=datawolff[2]
    msq=datawolff[3]
    #print m and msq
    #print("m=\t{:.5f}\tmsq=\t{:.5f}\n".format(m,msq))

    msqavg=np.sum(msq)/Nbin
    for tau in range(total_run):
        if(tau!=0):
            magshift=np.roll(m,-tau)
        else:
            magshift=m
        mtauavg=np.dot(m,magshift)/Nbin
        atau=abs(mtauavg/msqavg)
        Atau_list.append(atau)
    Atau_sum.append(np.sum(Atau_list))

#fit function
def func(x,a,b):
    return b*np.power(x,a)

#fit values
popt,pconv=curve_fit(func,L,Atau_sum)              #fit values
fit=func(L,popt[0],popt[1])                        #Use fit values in the funciton
z=popt[0]

for Nlin in Nlinrange:
    zlist.append(z)
y=np.column_stack((L,Atau_sum,fit,zlist))
np.savetxt(mc+"2bdata.dat",y)
