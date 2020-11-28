import numpy as np
import os
import sys    

Len=[4,4,4,3,3]
beta=[0.5,1,2,0.5,1.0]
#returns the spin of the configuration
def get_mag(spin):    
    return np.sum(spin)

#returns the energy of the configuration
def get_ener(spin):    
    #convert spin into L*L vector
    sp=np.reshape(spin,(L,L))

    #shift sp along x by 1 unit
    spx=np.roll(sp,1,axis=1)

    #shift sp along y by 1 unit
    spy=np.roll(sp,1,axis=0)

    #flatten shifted x
    spxflat=spx.flatten()

    #flatten shifted y
    spyflat=spy.flatten()

    #dot product of spin and flatten x
    dotx=-np.dot(spin,spxflat)

    #dot product of spin and flatten y
    doty=-np.dot(spin,spyflat)

    #total sum
    total=dotx+doty
    return total


#file to store L, beta, <e>
lbefile='lbetae.dat'
os.system('rm '+lbefile)
lbe=open(lbefile,'a')
for run in range(5):
    T=1.0/beta[run]
    L=Len[run]
    print("run\n",run+1)
    #file to save the data
    datafile='enumdata'+str(run)+'.dat'
    #open the file to write the T,E,C_v, X after each T
    f=open(datafile,'w')
    #check and create the file if doesnot exist
    #if os.path.exists(datafile):
    #    print(datafile+ " exists.")
    #    os.remove(datafile)
    #    print(datafile+ " removed.")

    #else:print(datafile+ "doesnot exist.")
    if run==0:f.write("{:}\t{:}\t{:}\t{:}\t{:}\n".format("T","E","C_v","U_L","X")) #header to the data file skipped in the plot
    if run==0:lbe.write("{}\t{}\t\t{}\t\t\t{}\n".format("L","beta","C_v","U_L"))

    Nsite =L*L #sites number
    Nstate=int(pow(2,Nsite))  #possible configuration size
    #print(40*'.')
    #declear z, e, e^2, m, m^2,energy to 0 to iterate the sum
    z=e=esq=m=msq=mqd=energy=0
    print("run for T\t",T)
    for state in range(Nstate):
        svec=np.array(list(np.binary_repr(state).zfill(Nsite))).astype(np.int8)
        #binary_repr returns the binary representation of the input number as a string
        #str.zfill(width)return a copy of the string left filled with ASCII '0'
        #digits to make a string of length width
        spin=2*svec-1    #gives the spin +-1
        #calculation of energy and magnetization of a state
        mag=get_mag(spin)
        ener=get_ener(spin)
        factor=np.exp(-ener/T)   #used later

        #Calculation of partition function
        z+=factor
        #calculation of energy
        energy+=ener*factor
        #e and m,e=E/Nsite, m=M/Nsite
        e+=1.0*(ener/Nsite)*factor
        m+=1.0*(mag/Nsite)*factor
        #e-square and m-square
        esq+=1.0*((ener/Nsite)**2)*factor
        msq+=1.0*((mag/Nsite)**2)*factor
        #m-qd
        mqd+=1.0*((mag/Nsite)**4)*factor

    #Calculation of expectation values
    e_expect=1.0*e/z                #<e>
    #print("<e>\t",e_expect)
    esq_expect=1.0*esq/z            #<e^2>
    m_expect=1.0*m/z                #<m>
    msq_expect=1.0*msq/z            #<m^2>
    mqd_expect=1.0*mqd/z            #<m^4>

    #calculation of C_v, X,U_L and E
    C_v=1.0*(Nsite/T**2)*(esq_expect-e_expect**2)       #specific heat
    U_L=(3.0/2.0)*(1-(mqd_expect/(3.*msq_expect**2)))   #Binder parameter
    X=1.0*(Nsite/T)*(msq_expect-m_expect**2)            #susceptibility
    E=1.0*energy/z                                      #Energy

    #write in the file
    f.write("{:0.3f}\t{:0.3f}\t{:0.3f}\t{:0.3f}\t{:0.3f}\n".format(T,E,C_v,U_L,X))
    lbe.write("{:}\t{:0.2f}\t{:8.5f}\t{:>0.5f}\n".format(L,beta[run],C_v,U_L))
    f.close() #close the file
lbe.close()
