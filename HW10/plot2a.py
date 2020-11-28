import numpy as np
import matplotlib.pyplot as plt


#Load the data
#mc='wolff'
mc="met"
for i in [5,10,20,40,80]:#np.arange(4,41,4):
    datafile=mc+'data'+str(i)+'.dat'
    data=np.loadtxt(datafile,unpack=True)
    theta=round(np.mean(data[3]),2)
    plt.plot(data[0],data[1],".",label="L="+str(i)+"\t"+r"$\Theta$="+str(theta))
    #fit function
    #plt.plot(data[0],data[2],"--",label="fit")
    plt.plot(data[0],data[2],"--")

plt.legend()
plt.grid(True)
plt.xlabel(r"$\tau$")
plt.ylabel(r"$A_{0}(\tau$)")
#plt.title("Auto-correlation function as a function of "+r"$\tau$"+" using Wolff's algorithm")
plt.title("Auto-correlation function as a function of "+r"$\tau$"+" using Metropolis algorithm")
plt.savefig(mc+".pdf")
#plt.savefig(mc+".png")
plt.show()
        



