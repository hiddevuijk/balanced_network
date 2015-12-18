import numpy as np
import matplotlib.pyplot as plt
import sys

nwe = np.genfromtxt("Eactivity.csv",delimiter=';')
nweh = np.genfromtxt("Eactivityh.csv",delimiter=';')



plt.plot(nwe,color='blue',label='Jee=Jie=1.0')
plt.plot(nweh,color='red',label='Jee=Jei=1.5')
plt.ylim([0,1])
plt.legend()

try:
	plt.savefig(sys.argv[1])
except:
	plt.show()



