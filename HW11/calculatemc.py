#This program calcualtes C_v, U_L and errors  based on the Wolff single cluster algorithm Monete Carlo simulation and saves in a file.

import os
import numpy as np
from scipy import stats

#input parameters
L=[4,4,4,3,3]
beta=[0.5,1,2,0.5,1.0]
Neql=[10**3]*5
Nmcs=[10**3]*3+[10**2]*2
Nbin=[10**4]*5
seed=[100]*5
latt='sqlatt_PBC'
alldata="resultwolffmc.dat"  #stores all the values from the simulation
#check and create the file if doesnot exist
if os.path.exists(alldata):
    print(alldata+ " exists.")
    os.remove(alldata)
    print(alldata+ " removed.")

else:print(alldata+ " doesnot exist.")
output=open(alldata,'a')

for i in range(5):
    print("run\t",i+1)
    if i==0:output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:}\n".format("L","beta","Neql","Nmcs","Nbin","C_v (mc)","C_v error","U_L (mc)","U_L error"))
    outfile='wolff'+str(i)+'.dat' #stores the data from 1 run
    #open the resultant data file
    data=np.loadtxt(outfile)
    #averages
    averages=np.mean(data,axis=0)
    #print("averages\n",averages)
    e_expect=averages[0]
    esq_expect=averages[1]
    m_expect=averages[2]
    msq_expect=averages[3]
    #standard error
    #stderror=stats.sem(data)
    #e_expect_error=stderror[0]
    #esq_expect_error=stderror[1]
    #m_expect_error=stderror[2]
    #msq_expect_error=stderror[3]

    if(beta[i]>1.5):m_expect=0  #showing strage behaviour at the lower temperatures
    ##calculatation of chi
    Nsite=L[i]**2 #number of sites=L^2
    #X=1.0*((1.0/(Nsite))*beta[i])*(msq_expect-(m_expect)**2)
    X=1.0*(Nsite*beta[i])*(1.0*msq_expect-(1.0*m_expect)**2)
    print("beta\t",beta[i])
    print("X(mc)\t",X)
    print("e(mc)\t",e_expect)
    output.write("{}\t{:0.2f}\t{}\t{}\t{}\t{:<8.5f}\t{:6.5f}\n".format(L[i],beta[i],Neql[i],Nmcs[i],Nbin[i],X,e_expect))
    #print("stndard errors\t",stderror) 
output.close()
