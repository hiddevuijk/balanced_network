import numpy as np
import matplotlib.pyplot as plt
import sys

nwe_acorr = np.genfromtxt("nwe_act_acorr.csv",delimiter=';')
nwe = np.genfromtxt("nwe_act.csv",delimiter=';')

nwe_acorr0 = nwe_acorr[0]
nwe_acorr /= nwe_acorr0

plt.plot(nwe_acorr[:30000],color='blue')

plt.show()


