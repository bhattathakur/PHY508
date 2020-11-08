import os
import numpy as np
from scipy import stats
L=[4,4,4,3,3]
beta=[0.5,1,2,0.5,1.0]
Neql=[10**3]*5
Nmcs=[10**4,10**4,10**4,10**3,10**3]
Nbin=[10**3]*5
seed=[100]*5
latt='sqlatt_PBC'
datafile="outputdata.out"
alldata="resultmc.dat"
for i in range(1):
    f=open("param.dat","w")
    if i==0:print("{}\t{}\t{}\t{}\t{}\t{}".format("L","beta","Neql","Nmcs","Nbin","SEED","latt_"))
    print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(L[i],beta[i],Neql[i],Nmcs[i],Nbin[i],seed[i],latt))
    #write in param.dat
    f.write("{}\t{}\t{}\t{}\t{}\t{}\t{:}".format(L[i],beta[i],Neql[i],Nmcs[i],Nbin[i],seed[i],latt))
    f.close()
    #use parama.dat as input
    os.system('./a.out')
    #open the resultant data file
    data=np.loadtxt(datafile)
    averages=np.mean(data,axis=0)
    print("averages\n",averages)
    stderror=stats.sem(data)
    print("standard errors\t",stderror) 
