import numpy as np
import matplotlib.pyplot as plt

acorr = np.genfromtxt("acorr.csv",delimiter=';')

plt.plot(acorr)
plt.show()



