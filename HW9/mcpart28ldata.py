import os
import numpy as np
from scipy import stats
L=8
Nsite=L**2 #number of sites=L^2
Neql=1000
Nmcs=100
Nbin=100
seed=100
latt='sqlatt_PBC'
datafile="mcl8data.dat"
outputdatafile="outputdata.out"
out=open(datafile,'a')
#alldata="resultmc.dat"
for i in np.arange(0.1,25,0.1):#range of Temp
    print("run for T=\t",i+1)
    beta=1.0/i                 #changing Temp to beta
    print("run for T=\t",i+1)
    f=open("param.dat","w")
    #outfile='output'+str(i)+'.dat' #stores the data from 1 run
    #if i==0:print("{}\t{}\t{}\t{}\t{}\t{}".format("L","beta","Neql","Nmcs","Nbin","SEED","latt_"))
    #print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(L[i],beta[i],Neql[i],Nmcs[i],Nbin[i],seed[i],latt))
    #write in param.dat
    f.write("{}\t{}\t{}\t{}\t{}\t{}\t{:}".format(L,beta,Neql,Nmcs,Nbin,seed,latt))
    f.close()
    #use parama.dat as input
    os.system('./a.out') #store the data in the file outputdata.out
    #copy the datainto different name
    #command='cp '+datafile+' '+outfile
    #os.system(command)
    #open the resultant data file
    data=np.loadtxt(outputdatafile)
    averages=np.mean(data,axis=0)
    e_expect=averages[0]
    esq_expect=averages[1]
    m_expect=averages[2]
    msq_expect=averages[3]
    print("averages\n",averages)
    stderror=stats.sem(data)
    e_expect_error=stderror[0]
    esq_expect_error=stderror[1]
    m_expect_error=stderror[2]
    msq_expect_error=stderror[3]
    ##calculatation of chi
    X=1.0*(Nsite*beta)*(msq_expect-m_expect**2)
    print("X(mc)\t",X)
    #print("e(mc)\t",e_expect)
    #calculation of C_v
    C_v=1.0*(Nsite*beta*beta)*(esq_expect-e_expect**2)
    print("C_v\t",C_v)
    print("standard errors\t",stderror) 
    #writing in the file
    if i==0:out.write("{}\t{}\t{}\n".format("T","C_v","X"))
    out.write("{:0.3f}\t{:0.3f}\t{:0.3f}\n".format(i,C_v,X))
out.close()
