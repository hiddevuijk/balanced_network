import numpy as np
import matplotlib.pyplot as plt

me = np.genfromtxt("me.csv",delimiter=';')
mi = np.genfromtxt("mi.csv",delimiter=';')
mo = np.genfromtxt("mo.csv",delimiter=';')


plt.plot(mo,me,color='blue')
plt.plot(mo,mi,color='red')
plt.show()


