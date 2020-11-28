#This program run the Monte Carlo simulation based on Wolff single cluster algorithm for various input parameters and save the average values of <e> <e^2> <m> ,<m^2>,<m^4>, <average size> in a file

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
datafile="wolfoutput.out"  #output from cpp run
synno=10                   #number of synthetic files
saveoutput="dataerror.txt" #final data with errors

#checks if file exists, delete if so
def checkfile(f):
    if os.path.exists(f):
        print(f,"\t exist")
        os.remove(f)
        print(f, "\t deleted")
    else:print(f,"\t doesnot exists...created\n")

#creates a synthetic files number of times for a file 
def synthetic_data(datafile,i):                      #index, temp,filename
    print("Inside synthetic data\n")
    orgdata=np.loadtxt(datafile)                     #read the data file
    rows,columns=np.shape(orgdata)
    for j in range(i):
        idx=datafile.find('.')
        newname=datafile[:idx]+str(j)+".dat"
        print("newname\t",newname)
        #if j==0:newname=newname+str(j)+".dat"
        checkfile(newname)
        #open the file for read
        f=open(newname,'a')
        #loops number of rows times, get the random row and store in the synethetic file
        for j in np.arange(rows):
            #print("j\t",j)
            ran_row=np.random.randint(rows)   #randint(n)->[0,n)
            ran_row_value=np.asarray(orgdata[ran_row,:])
            np.savetxt(f,ran_row_value,fmt="%10.5f",delimiter='\t',newline=' ')
            f.write("\n")
        f.close()
#####################################################################################################
#read the output file and calculate specific heat and binder parameter
def specific_heat_binder_parameter(infile,T):
   data=np.loadtxt(infile,unpack=True) 
   #averages of e,e^2, m, m^2, m^4
   e=np.mean(data[0])
   e2=np.mean(data[1])
   m=np.mean(data[2])
   m2=np.mean(data[3])
   m4=np.mean(data[4])

   c_v=(1.0*Nsite/T**2)*(e2-e**2)                    #specific heat
   U_L=(3./2.)*(1.-(m4/(3.*m2**2)))                  #Binder Parameter
   #print("c_v\t{:0.3f}\tU_L\t{:0.3f}".format(c_v,U_L))
   return c_v,U_L

#########################################################################################
#Create the original datafiles
def crete_data():
    for i in range(len(L)):
        print("run\t",i+1)
        f=open("param.dat","w")  #input to the cpp progrm
        outfile='wolff'+str(i)+'.dat' #stores the data from 1 run
        if i==0:print("{}\t{}\t{}\t{}\t{}\t{}".format("L","beta","Neql","Nmcs","Nbin","SEED","latt_"))
        print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(L[i],beta[i],Neql[i],Nmcs[i],Nbin[i],seed[i],latt))
        #write in param.dat
        f.write("{}\t{}\t{}\t{}\t{}\t{}\t{:}".format(L[i],beta[i],Neql[i],Nmcs[i],Nbin[i],seed[i],latt))
        f.close()
        #use parama.dat as input
        os.system('./a.out>demo.log')
        #copy the datainto different name
        command='cp '+datafile+' '+outfile
        os.system(command)

#files from the c++ program with above parameters
input_files=["wolff0.dat", "wolff1.dat", "wolff2.dat", "wolff3.dat", "wolff4.dat"]

#get specific heat and binder parameter
k=0
C_v=[];U_L=[]
C_v_error=[];U_L_error=[]

for file in input_files:
    #get synthetic files
    synthetic_data(file,synno)
    #f=np.loadtxt(file)
    Nsite=L[k]*L[k]
    t=1./beta[k]      #temperature
    c,u=specific_heat_binder_parameter(file,t)
    C_v.append(c);U_L.append(u)
    print("k={:2.2f}\tc_v={:2.5f}\tu_l={:2.5f}".format(k,c,u))
    idx=file.find('.')
    print("idx\t",idx)
    prefix=file[:idx]
    print("prefix\t",prefix)
    cvtemp=[];ultemp=[]
    for i in range(synno):
        datafile=prefix+str(i)+".dat"
        print("datafile\t",datafile)
        c,u=specific_heat_binder_parameter(datafile,t)
        print("k={:2.2f}\tc_v={:2.5f}\tu_l={:2.5f}".format(k,c,u))
        cvtemp.append(c);ultemp.append(u)
        #go to each synthetic file get cv,ul and append to the given
    print("c-v\t",cvtemp)
    print("u-l\t",ultemp)
    #get sem
    k+=1
    cvrms=stats.sem(cvtemp)
    C_v_error.append(cvrms)
    ulrms=stats.sem(ultemp)
    U_L_error.append(ulrms)
#print("."*40) 
#print("C_v\t",C_v)
#print("C_v_error\t",C_v_error)
#print("U_L\t",U_L)
#print("U_L_error\t",U_L_error)

#creating the column stack
dataerror=np.column_stack((L,beta,Neql,Nmcs,Nbin,C_v,C_v_error,U_L,U_L_error))
print("all data\n",dataerror)
#saving to the file
head="L beta Neql Nmcs Nbin C_v  C_v_error  U_L  U_L_error"
np.savetxt(saveoutput,dataerror,fmt='%d %0.2f %d %d %d %.5f %.5f %.5f %.5f',header=head)
