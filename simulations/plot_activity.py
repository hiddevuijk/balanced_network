import numpy as np
import matplotlib.pyplot as plt
import sys

nwe = np.genfromtxt("Eactivity.csv",delimiter=';')
nwi = np.genfromtxt("Iactivity.csv",delimiter=';')

plt.plot(nwe,color='blue')
plt.plot(nwi,color='red')
plt.show()


