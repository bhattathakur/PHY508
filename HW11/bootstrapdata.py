"""This program reads the number of rows in the original file, randomly selects the rows for the times equal to the number of rows in the original data file, saves the new rows in the new synthetic file """
#####################################################################################################

#libraries
import numpy as np
import os

#parameters
L=4                         #lattice size
Nsite=L**2                  #number of sites=L^2
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

def createparam(temp):
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

def synthetic_data(T,i):
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

#get  the specific heat and binder parameter for the given input file with 5 columns (<e>,<e^2>,<m>,<m^2>,<m^4>)
# C_v=(N_s/T^2)(<e^2>-<e>^2)
# U_L=(3/2)(1-<m^4>/3<m^2>^2)

def specific_heat_binder_parameter(infile,T):
   data=np.loadtxt(infile,unpack=True) 
   #averages of e,e^2, m, m^2, m^4
   e=np.mean(data[0])
   e2=np.mean(data[1])
   m=np.mean(data[2])
   m2=np.mean(data[3])
   m4=np.mean(data[4])

   c_v=(1.0*Nsite/T**2)*(1.*e2-e**2)      #specific heat
   U_L=(3./2.)*(1.-1.0*m4/(3.*m2**2))                  #Binder Parameter

   return c_v,U_L

#returns the error values in C_v and U_L
def get_error_cv():
    c_vlist=[]
    u_llist=[]
    for i in range(len(synfileno)):                    #number of synthetic files
        c_v,u_l=specific_heat_binder_parameter(synfile,T)
        c_vlist.append(c_v)
        u_llist.append(u_l)
        
    return root_mean_square(c_vlist),root_mean_square(u_l)

#cv,ul
cv=[]
ul=[]
c_verror=[]
u_lerror=[]
Te=[]
    
#main program
for T in range(1,16,5):
    Te.append(T)
    createparam(T)                       #create a param.dat file, run c++ and 
                                         #create the wolfoutput.dat file 
                                         #for given temperature and predefined parameters
    c_v,u_l=specific_heat_binder_parameter(datafile,T) #cv,ul from the original file
    print("c_v\n",c_v)
    print("u_l\n",u_l)
    cv.append(c_v)
    ul.append(u_l)
    c_vlist=[]
    u_llist=[]
    for i in range(synfileno):
        print("sythetic-file\t",i+1)
        synthetic_data(T,i)                               #synthetic file with given temperature
        synfile=suffix+str(T)+str(i)+".dat"               #synthetic file created
        cv,ul=specific_heat_binder_parameter(synfile,T)   #gives specific heat and binder parameter
        print("cv\n",cv)
        print("ul\n",ul)
        c_vlist.append(cv)
        u_llist.append(ul)
    #get root mean square
    print("c_v\n",c_vlist)
    print("u_l\n",u_llist)
    cv_rms=root_mean_square(c_vlist)
    print("cv-rms\n",cv_rms)
    ul_rms=root_mean_square(u_llist)
    print("cv-rms\n",ul_rms)
    print("c_v\n",c_vlist)
    c_verror.append(cv_rms)
    u_lerror.append(ul_rms)
    print("deleting all the .dat files....")
    #os.system('rm *.dat')
#create a column_stack    
#all=np.column_stack((Te,c_v,c_verror,u_l,u_lerror))
print("{}\t{}\t{}\t{}\t{}".format("T","C_v","C_v-error","U_L","U_L-error"))
print("Te\t",Te)
print("C_v\t",cv)

#save the file
#np.save("allthedata.dat",all)
