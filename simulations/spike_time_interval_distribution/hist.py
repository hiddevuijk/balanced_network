import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.optimize import curve_fit
from scipy.misc import factorial

save = False
name = ''
if len(sys.argv) == 2:
	name = sys.argv[1]
	save = True


bins = 20

nm = open("names.txt",'r')
names = []
while 1:
	line = nm.readline()
	if not line: break
	names.append(line[:-1])
stis=[]
for n in names:
	stis.append(np.genfromtxt(n,delimiter=';'))

n = len(stis)

def poisson(x,lamb):
	return np.exp(-1.*lamb*x)*lamb


for i in range(n):
	#l=1/mean = lambda in exp distr
	l = 1.*len(stis[i][:-1])/sum(stis[i][:-1])
	plt.subplot(1,n,i)
	entries, bin_edges, patches=plt.hist(stis[i][:-1],bins, normed=True)
	x = np.linspace(0,bin_edges[-1],1000)
	plt.plot(x,poisson(x,l))
	plt.title(names[i])

if not save:
	plt.show()
else:
	figure = plt.gcf()
	figure.set_size_inches(5.5*n,8)
	plt.savefig(name)


