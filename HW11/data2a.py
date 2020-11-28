#This program generates the data for plotting <m^2> as 1/L Wolff single cluster algorithm Monete Carlo simulation and saves in a file.

#load libraries
import os
import numpy as np
from scipy import stats

#input parameters

#critical temperature=2.269
#temperatue above and below the critical temperaturs
temp=[1.25,1.5,1.75,2.0,2.269,2.75,3.0,3.25,3.5]
Neql=1000
Nmcs=100
Nbin=100
seed=100
latt='sqlatt_PBC'
alldata="m2l.dat"         #stores all the values from the simulation
wolfdata="wolfoutput.out" #output from c++ program

#check and create the file if doesnot exist
if os.path.exists(alldata):
    print(alldata+ " exists.")
    os.remove(alldata)
    print(alldata+ " removed.")

else:print(alldata+ " doesnot exist.")

#Create the original datafile for given temperature and lattice size using the other parameters as given.
def crete_data(temperature,size):
    beta=1.0/temperature
    f=open("param.dat","w")     #input to the cpp progrm
    print("{}\t{}\t{}\t{}\t{}\t{}".format("L","beta","Neql","Nmcs","Nbin","SEED","latt_"))
    print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(size,beta,Neql,Nmcs,Nbin,seed,latt))

    #write in param.dat
    f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(size,beta,Neql,Nmcs,Nbin,seed,latt))
    f.close()

    #run c++ program using parama.dat as input
    os.system('./a.out>m2l.log')

#iterates for various values of L
j=0
for temperature in temp:
    print("run\t",j+1)
    l=[];m2=[]
    for lattsize in range(4,21):
        print("length\t",lattsize)
        l.append(lattsize)               #Adding lattisize in the list
        #run c++ program and produce a file wolfoutput.out
        crete_data(temperature,lattsize)      #temperature,lattice-size
        
        #read the output file from the c++ program
        data=np.loadtxt(wolfdata,unpack=True)         #e,e^2,m,m^2,m^4
        msq=np.mean(data[3])                          #<m^2>
        m2.append(msq)                                #adding m^2 in m2 list
    lm2data=np.column_stack((l,m2)) #adding lists in the 2d array
    print("lmt2data\n",lm2data)
    filename="lm2data"+str(j)+".dat"
    head="temperature->"+str(temp[j])
    np.savetxt(filename,lm2data,fmt='%d %0.5f',header=head)
    j+=1
