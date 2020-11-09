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
datafile="outputdata.out"
#alldata="resultmc.dat"
for i in range(5):
    if i!=2:continue
    print("run\t",i+1)
    f=open("param.dat","w")
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
    #open the resultant data file
    data=np.loadtxt("output2.dat")
    averages=np.mean(data,axis=0)
    #e_expect=averages[0]
    #esq_expect=averages[1]
    #m_expect=averages[2]
    #msq_expect=averages[3]
    print("averages\n",averages)
    #stderror=stats.sem(data)
    #e_expect_error=stderror[0]
    #esq_expect_error=stderror[1]
    #m_expect_error=stderror[2]
    #msq_expect_error=stderror[3]
    ##calculatation of chi
    #Nsite=L[i]**2 #number of sites=L^2
    #X=1.0*(Nsite*beta)*(msq_expect-msq_expect**2)
    #print("X(mc)\t",X)
    #print("e(mc)\t",e_expect)
    #print("standard errors\t",stderror) 
