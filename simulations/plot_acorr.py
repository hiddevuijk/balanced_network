import numpy as np
import matplotlib.pyplot as plt
import sys

acorr = np.genfromtxt("acorr.csv",delimiter=';')

plt.plot(acorr)
plt.show()

