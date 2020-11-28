#create multiple subplots with <m^2> vs 1/L vs  for different value of temperature

#libraries
import numpy as np
import matplotlib.pyplot as plt

#datafiles=[ "lm2data0.dat", "lm2data1.dat", "lm2data2.dat", "lm2data3.dat", "lm2data4.dat", "lm2data5.dat", "lm2data6.dat",\
#        "lm2data7.dat", "lm2data8.dat"]

temperature=[1.25,1.5,1.75,2.0,2.269,2.75,3.0,3.25,3.5]     #temperatureo
colors=['b','g','r','c','m','y','k','b','m']                #subplots color
lattsize=np.arange(4,21)                                    #lattice sizes
xname=r'$\frac{1}{L}$'                                      #xvalue 1/L
yname=r'$\langle m^{2}\rangle$'                             #yvalue <m^2>

#rows=3,columns=3 for subplots

fig,axs=plt.subplots(3,3,constrained_layout=True,figsize=(12,8))
for row in range(3):
    for column in range(3):
        num=row*3+column                                    #file number
        #print("test\t",num)
        filename="lm2data"+str(num)+".dat"                  #data file
        #print("filename\t",filename)
        data=np.loadtxt(filename,unpack=True)               #input datafile
        x=1./data[0]        #1/L
        y=data[1]           #<m^2>
        axs[row,column].plot(x,y,'.'+str(colors[num]),label='T='+str(temperature[num])) #plot in each subplot
        if num==4:axs[row,column].set_facecolor("lightcyan") #different color to the critical temperature

#set same labels,legend and grid  for all subplots
for ax in axs.flat:
    ax.set(xlabel=xname,ylabel=yname)
    ax.legend()
    ax.grid()

#plt.tight_layout()
fig.suptitle(yname+" as the function of "+xname+" for various temperature"+r" T")
plt.savefig("m2vsl.pdf",dpi=800)
plt.savefig("m2vsl.svg",dpi=800)
plt.show()


