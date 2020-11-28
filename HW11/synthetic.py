"""This program reads the number of rows in the original file, randomly selects the rows for the times equal to the number of rows in the original data file, saves the new rows in the new synthetic file """
#####################################################################################################

#libraries
import numpy as np
import os
from scipy import stats

#parameters
#L=20                        #lattice size
#Nsite=L**2                  #number of sites=L^2
Neql=1000
Nmcs=100
Nbin=100
seed=100
latt='sqlatt_PBC'
datafile="wolfoutput.out"   #original output from c++ program
synfileno=10                #number of synthetic files for each file
suffix="synthetic"          #synthetic data suffix

#This function creates the param.dat file with given temperature and various other 
#fixed parameters values as above 

def createparam(temp,L):
    print("."*40)
    print("creating param.dat file")
    f=open("param.dat","w")
    #write in param.dat
    f.write("{}\t{}\t{}\t{}\t{}\t{}\t{:}".format(L,1.0/temp,Neql,Nmcs,Nbin,seed,latt))
    f.close()
    os.system('./a.out>test.log')       #create a datafile with above parameters

#checks if file exists, delete if so
def checkfile(f):
    if os.path.exists(f):
        print(f,"\t exist")
        os.remove(f)
        print(f, "\t deleted")
    else:print(f,"\t doesnot exists...created")

#returns root mean square of the given data
def root_mean_square(data):
    #data=np.array(data)
    #sum=0
    print("list inside sqrt\t",data)
    print("length:\t",len(data))
    square=[i**2 for i in data]
    sum1=sum(square)
    n=len(data)
    return (sum1/n)**0.5

#creates a synthetic file for a file of given suffix and index number

def synthetic_data(i,T):                  #index, temp
    orgdata=np.loadtxt(datafile)          #read the data file
    #print("run\t",i+1)
    outfile=suffix+str(T)+str(i)+".dat"
    print("outputfile\t",outfile)
    checkfile(outfile)

    #open the file for read
    f=open(outfile,'a')
    rows,columns=np.shape(orgdata)

    #loops number of rows times, get the random row and store in the synethetic file
    for j in np.arange(rows):
        #print("j\t",j)
        ran_row=np.random.randint(rows)   #randint(n)->[0,n)
        ran_row_value=np.asarray(orgdata[ran_row,:])
        #print("random row\t",ran_row_value)
        np.savetxt(f,ran_row_value,fmt="%10.5f",delimiter='\t',newline=' ')
        f.write("\n")
    f.close()

def specific_heat_binder_parameter(infile,T):
   data=np.loadtxt(infile,unpack=True) 
   #averages of e,e^2, m, m^2, m^4
   e=np.mean(data[0])
   e2=np.mean(data[1])
   m=np.mean(data[2])
   m2=np.mean(data[3])
   m4=np.mean(data[4])

   c_v=(1.0*Nsite/T**2)*(e2-e**2)      #specific heat
   U_L=(3./2.)*(1.-(m4/(3.*m2**2)))                  #Binder Parameter
   #print("c_v\t{:0.3f}\tU_L\t{:0.3f}".format(c_v,U_L))
   return c_v,U_L

#returns the root means square of the list
def rms(li):
    sq=[x**2 for x in li]
    re=sum(sq)/len(sq)
    return re**0.5


#synthetic file numbers
synrun=20
os.system("rm -f *.dat")
for size in [5,10,15,20,25,30]:
    if size!=30:continue
    Nsite=size*size
    t_v=[]
    c_v=[];cvrms=[]
    u_l=[];ulrms=[]
    print("run for size\t",size)
    for T in np.arange(0.1,25,0.1):
        createparam(T,size)
        #specific heat binder parameter
        c,u=specific_heat_binder_parameter(datafile,T)
        #print("{:0.2f}\t{:0.2f}".format(c,u))
        c_v.append(c);u_l.append(u)
        cv=[];ul=[]
        for i in range(synrun):
            synthetic_data(i,T)
            #read the data file
            outfile=suffix+str(T)+str(i)+".dat"
            #command='head '+outfile
            #os.system('head '+outfile)
            #data=np.loadtxt(outfile,unpack=True)
            data=np.loadtxt(outfile)
            #print(data[:5])
            print(np.mean(data[0]))

            #specific heat binder parameter
            c,u=specific_heat_binder_parameter(outfile,T)
            print("{:0.2f}\t{:0.2f}".format(c,u))
            cv.append(c);ul.append(u)
       # print("cv list\n",cv)
       # print("ul list\n",ul)

        #calculate the root mean square of c_v and u_l
        #cvr=rms(cv);ulr=rms(ul)
        cvr=stats.sem(cv);ulr=stats.sem(ul)
        cvrms.append(cvr);ulrms.append(ulr)
       # print("cvrms\n",cvr)
       # print("ulrms\n",ulr)
        t_v.append(T)
   # print("Final Here")    
   # print("*"*40)    
   # print("T_v list\n",t_v)
   # print("c_v list\n",c_v)
   # print("c_vrms list\n",cvrms)
   # print("u_l list\n",u_l)
   # print("u_lrms list\n",ulrms)
    #stack
    table=np.column_stack((t_v,c_v,cvrms,u_l,ulrms))
    #print("tabular\n",table)
    os.system("rm -f *.dat")
    np.savetxt("myfile"+str(size)+".txt",table,fmt="%1.4e")




