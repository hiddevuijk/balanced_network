import numpy as np
import matplotlib.pyplot as plt

me = np.genfromtxt("me.csv",delimiter=';')
mi = np.genfromtxt("mi.csv",delimiter=';')
mo = np.genfromtxt("mo.csv",delimiter=';')
mme = np.genfromtxt("mme.csv",delimiter=';')
mmi = np.genfromtxt("mmi.csv",delimiter=';')
mmo = np.genfromtxt("mmo.csv",delimiter=';')


plt.plot(mo,me,color='blue')
plt.plot(mo,mi,color='red')
plt.plot(mmo,mme,color='blue',linestyle='--')
plt.plot(mmo,mmi,color='red',linestyle='--')
plt.show()


