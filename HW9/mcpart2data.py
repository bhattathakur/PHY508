#This program calculates the C_v, X for various values of temperatures T using the Monte Carlo based on Metropolis algorithm.

import os
import time
import numpy as np
from scipy import stats

L=4    #lattice size
Nsite=L**2 #number of sites=L^2
Neql=1000
Nmcs=100
Nbin=100
seed=100
latt='sqlatt_PBC'

datafile="mcl4data.dat"                 #stores the averages after each run
outputdatafile="outputdata.out"         #output from c++ program
out=open(datafile,'a')
for i in np.arange(0.1,25,0.1):#range of Temp
    print("run for T=\t",i+1)
    time.sleep(0.5)
    beta=1.0/(1.0*i)                 #changing Temp to beta
    f=open("param.dat","w")
    #write in param.dat
    f.write("{}\t{}\t{}\t{}\t{}\t{}\t{:}".format(L,beta,Neql,Nmcs,Nbin,seed,latt))
    f.close()

    #use parama.dat as input
    os.system('./a.out') #store the data in the file outputdata.out
    #open the resultant data file

    data=np.loadtxt(outputdatafile)     #load data file
    averages=np.mean(data,axis=0)
    e_expect=averages[0]
    esq_expect=averages[1]
    m_expect=averages[2]
    msq_expect=averages[3]
    #print("averages\n",averages)
    #stderror=stats.sem(data)
    #e_expect_error=stderror[0]
    #esq_expect_error=stderror[1]
    #m_expect_error=stderror[2]
    #msq_expect_error=stderror[3]

    ##calculatation of chi
    if i<2: m_expect=0  #showing the anomoalous behaviour at the lower temperatuers
    X=1.0*(Nsite*beta)*(msq_expect-m_expect**2)
    #print("X(mc)\t",X)
    #print("e(mc)\t",e_expect)
    #calculation of C_v
    C_v=1.0*(Nsite*beta*beta)*(esq_expect-e_expect**2)
    #print("C_v\t",C_v)
    #print("standard errors\t",stderror) 
    #writing in the file
    if i==0:out.write("{}\t{}\t{}\n".format("T","C_v","X"))
    out.write("{:0.3f}\t{:0.3f}\t{:0.3f}\n".format(i,C_v,X))
out.close()
