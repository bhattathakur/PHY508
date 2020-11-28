import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os

fig = plt.figure()

#parameters
Nlin=32
beta=0.04
#beta=-1
Neql=1
Nmcs=1
Nbin=1
SEED=100
latt="sqlatt_PBC"
#removing data.out
os.system('rm data.out')

#writing the parameters in the file 
parameterfile='param.dat'
paramfile=open(parameterfile,'w');
paramfile.write("{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n".format(Nlin,beta,Neql,Nmcs,Nbin,SEED,latt))
#print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(Nlin, beta, Neql, Nmcs, Nbin, SEED,latt))
paramfile.close()

#run the a.out once to get the value of data.out also produces conf.out
os.system('./a.out')

#load the given data and show in the image
data=np.loadtxt('conf.out')
im=plt.imshow(data,interpolation='nearest',animated=True)

#update function
def updatefig(i):
    os.system('cp conf.out conf.in')
    Neql=0      #for this condition conf.in is read
    paramfile=open(parameterfile,'w') #file to store parameters
    paramfile.writelines("{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n".format(Nlin,beta,Neql,Nmcs,Nbin,SEED+i,latt))
    paramfile.close()

    os.system('./a.out') #runs the c++ program
    data=np.loadtxt('conf.out')
    im.set_array(data)
    plt.title(str(i))
    return im

ani=animation.FuncAnimation(fig, updatefig,np.arange(600), interval=10)
plt.show()
