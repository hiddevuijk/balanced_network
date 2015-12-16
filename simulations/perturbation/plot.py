import numpy as np
import matplotlib.pyplot as plt
import sys

nwe = np.genfromtxt("Eactivity.csv",delimiter=';')

plt.plot(nwe)
plt.ylim([0,1])
try:
	plt.savefig(argv[1] + "png")
except:
	plt.show()



