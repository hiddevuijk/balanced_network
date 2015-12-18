import numpy as np
import matplotlib.pyplot as plt
import sys


N = int(np.genfromtxt("N.csv",delimiter=';'))
tmax = int(np.genfromtxt("tmax.csv",delimiter=';'))
spike = np.genfromtxt("spike.csv",delimiter=';')
Ein = np.genfromtxt("Ein.csv",delimiter=';')
Iin = np.genfromtxt("Iin.csv",delimiter=';')
Tin = Ein+Iin



grid = (4,1)

plt.subplot2grid(grid,(0,0),colspan=1,rowspan=3)
plt.plot(Ein,color='black')
plt.plot(Iin,color='black')
plt.plot(Tin,color='black')
plt.plot([1]*len(Tin),linestyle='--',color='black')
plt.ylabel("current")

plt.tick_params(axis='x',which='both',labelbottom='off')

plt.subplot2grid(grid,(3,0),colspan=1,rowspan=1)
plt.plot(spike,color='black')
plt.xlabel("time")
plt.ylabel("spikes")
plt.tick_params(axis='y',which='both',left='off',labelleft='off',right='off')
tk = range(0,tmax+1,10*N)
lb = range(0,tmax/N + 1,10)
plt.xticks(tk,lb)

plt.show()

