import numpy as np
import matplotlib.pyplot as plt
import sys

inname = sys.argv[1];

nwe_acorr = np.genfromtxt(inname, delimiter=';')
nwe_acorr = nwe_acorr[:100000]
c0 = nwe_acorr[0];
plt.plot(nwe_acorr/c0)

plt.savefig(sys.argv[2]+".png")




