import sys
import re
from StringIO import StringIO
import numpy as np
import matplotlib

matplotlib.use('pdf')
import matplotlib.pyplot as plt

chain = str(sys.argv[2])

with open(str(sys.argv[1])+'/ip_raw_chain_sub'+chain+'.m', 'r') as fin:
    header = fin.readline()

    # Sort out the dimension of the 2D array of samples
    dims = re.findall(r'\d+', header)
    dims = [int(dim) for dim in dims]
    dims[0] -= 1

    # Ignore first line because laziness
    # for i in range(20):
    samples = fin.readline()

    # Now read all and ignore last line
    samples = fin.readlines()
    samples = samples[0:-1]

    # Now join all up and use magic numpy read help
    samples = ' '.join(samples)
    samples = np.loadtxt(StringIO(samples))

num_samples = samples.shape[0]
N = 60000

mins = np.array([0.95, 0.9, 0.9, 0.9, 0.9])
maxs = np.array([1.05, 1.1, 1.1, 1.1, 1.1])

samples[:,3] = mins[3] + (samples[:,3] * (maxs[3] - mins[3]))
samples[:,4] = mins[4] + (samples[:,4] * (maxs[4] - mins[4]))

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.scatter(samples[N:,3], samples[N:,4])
ax.set_xlabel('k_xkwlx')
ax.set_ylabel('k_cd')
fig.savefig('joint.pdf')

ax.cla()

eta = 0.906 * samples[:,3] + 0.478 * samples[:,4]

ax.hist(eta[N:], bins=20)
fig.savefig('eta.pdf')

print "eta mean:", np.mean(eta[N:])
print "eta std:", np.std(eta[N:], ddof=1)  # should be 0.00791
