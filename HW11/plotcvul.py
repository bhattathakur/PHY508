"""
This program plots the specific heat C_v and Binder paramter U_L along with associted errors. Errors are calulated by the 
Bootstrap method.
"""
#load libraries
import numpy as np
import matplotlib.pyplot as plt
import os

#datafile suffix
suffix="myfile"
xlabels="T"
ylabels=[r"$C_{v}$",r"$U_{L}$"]
save=["cvasT.pdf","ulasT.pdf"]
colors=['b','g','y','c','m','k']
loc="(2.2379,0.9295)"
for j in range(2):
    for i in range(5,31,5):
        #if(i!=5):continue
        datafile=suffix+str(i)+".txt"               #datafile

        #print(suffix+str(i)+".txt")
        #check if the file exists
        if(os.path.exists(datafile)):print(datafile+" exists")
        else: print(datafile+" doesnot exist")

        #load the datafile
        T,c_v,c_verror,u_l,u_lerror=np.loadtxt(datafile,unpack=True)

        #print the data
        #print("T\n",T)
        #print("c_v\n",c_v)
        #print("c_v\n",c_v)
        col=int(i/5)-1
        #print("col\t",col)
        form='o-'+colors[col]
        #print("form\t",form)
        if j==0:plt.errorbar(T,c_v,fmt=form,ms=0.95,yerr=c_verror,ecolor='red',elinewidth=1,capsize=1,barsabove=True,label='L='+str(i))

        else:
            plt.errorbar(T,u_l,fmt=form,ms=0.95,yerr=u_lerror,ecolor='red',elinewidth=1,capsize=1,barsabove=True,label='L='+str(i))
            #plt.xlim(0.0,3.0)
            #plt.ylim(0.75,1.05)
            plt.text(2.3,0.93,loc)
    #grids
    plt.legend()
    plt.grid(True)
    plt.minorticks_on()
    plt.grid(which='major',linestyle=":",linewidth='0.5',color='k')
    plt.grid(which='minor',linestyle=":",linewidth='0.25',color='k')
    plt.xlabel(xlabels)
    plt.ylabel(ylabels[j])
    #plt.set_markersize(0.01)
    plt.xlabel(xlabels)
    plt.savefig(save[j])
    plt.show()

        
        
