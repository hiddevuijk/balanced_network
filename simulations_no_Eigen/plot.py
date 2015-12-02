import numpy as np
import matplotlib.pyplot as plt
import sys


nwe = np.genfromtxt("Eactivity.csv", delimiter=';')
nwi = np.genfromtxt("Iactivity.csv", delimiter=';')
Ein = np.genfromtxt("Ein.csv", delimiter=';')
Iin = np.genfromtxt("Iin.csv", delimiter=';')
spike = np.genfromtxt("spike.csv",delimiter=';')
the = np.genfromtxt("the.csv",delimiter=';')
#thi = np.genfromtxt("thi.csv",delimiter=';')
inode = int(np.genfromtxt("inode.txt"))

plt.subplot(211)
plt.plot(Ein/the[inode],color='blue')
plt.plot(Iin/the[inode],color='red')
plt.plot((Ein+Iin)/the[inode],color='black')
plt.plot([1]*len(Ein),color='black',linestyle=':')

plt.subplot(212)
plt.plot(spike)

plt.show()
exit(0)





