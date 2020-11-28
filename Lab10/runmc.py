#This program runs the a.out - the result of the c++ program (ising_ssf.cpp) and calcuates <e> <e^2> <m> and <m^2> in a file data.out. 
#Steps used: for a single run of a.out very high T (very low beta/)very low T(high beta) was choosen and random configuration conf.out is obtained
# Nmcs=1, bin=1, Neql=1, SEED=1, square PBC lattice is choosen
# conf.out is copied as conf.in
# For subsequent runs beta is change to reasonalbe value with other remaining same

import os
import numpy as np
import matplotlib.pyplot as plt
run=100
#run the mc for multiple times

alldata="data.out"  #stores all the values from the simulation when run ./a.out
#check and create the file if doesnot exist
if os.path.exists(alldata):
    print(alldata+ " exists.")
    os.remove(alldata)
    print(alldata+ " removed.")
else:print(alldata+ " doesnot exist.")
for i in range(run):
    print("run\t",i+1)
    os.system('./a.out')
    os.system('cp conf.out conf.in')
