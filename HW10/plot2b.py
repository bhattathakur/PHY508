import numpy as np
import matplotlib.pyplot as plt


#Load the data
wolfdata='wolff2bdata.dat'
metdata='met2bdata.dat'

#reading two datafiles
data0=np.loadtxt(wolfdata,unpack=True) #L, theta fittheta, z
data1=np.loadtxt(metdata,unpack=True)

#getting the value of z
wolf_z=round(np.mean(data0[3]),4)
met_z=round(np.mean(data1[3]),4)

#plotting and fitting wolff data
#plt.plot(data0[0],data0[1],".",label=r"$z_{Wolff}$ = "+str(wolf_z))
#plt.plot(data0[0],data0[2],"--")

#plotting and fitting metropolis data
plt.plot(data1[0],data1[1],".",label=r"$z_{ss-met}$ = "+str(met_z))
plt.plot(data1[0],data1[2],"--")

#plt.yscale("log")
#plt.xscale("log")
plt.legend()
plt.grid(True)
plt.xlabel(r"$L$")
plt.ylabel(r"$\Theta$")
#plt.title("Plot of autocorrection time as a function of L using Wolff's algorithm")
plt.title("Plot of autocorrection time as a function of L using Metropolis algorithm")
#plt.title("Autocorrection time as a function of L with Wolff & Metropolis algorithm")
plt.savefig("met2b.pdf")
plt.savefig("met2b.png")
plt.show()
        



