import numpy as np
import matplotlib
matplotlib.use('pdf')
from matplotlib import pyplot as plt
from matplotlib import rc

rc('text', usetex=True)

min_x = -10
max_x = 10
size_x = 21
x = np.linspace(min_x,max_x, num=size_x)
sig = np.zeros(size_x)
np.savetxt('xs.txt',x, fmt='%5.2f')
m = 3
c = 5
mu, sigma = m*x + c, 1.0
sig[:] = sigma
np.savetxt('sigmas.txt',sig, fmt='%5.2f')
print mu.shape;
l = np.prod(x.shape)
y = np.zeros(l)
np.random.seed(1234)

for i in range(l):
	y[i] = np.random.normal(mu[i],sigma,1)

np.savetxt('obs.txt',y, fmt='%5.2f')

fig = plt.figure()
ax = fig.add_subplot(111)
line1, = ax.plot(x,mu,linewidth=2)
line2, = ax.plot(x,y,'r', linestyle='None', marker="*")
ax.set_xlim([x[0]-1, x[-1]+1])
ax.legend([line1, line2], [r'$\mathrm{Model}$', r'$\mathrm{Data}$'], loc=2)
fig.savefig('truth_data')



#y = np.random.normal(mu,sigma,1)
#print y.shape;

