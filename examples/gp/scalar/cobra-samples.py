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

mins = np.array([0.95, 0.9, 0.9, 0.9, 0.9])
maxs = np.array([1.05, 1.1, 1.1, 1.1, 1.1])

filenames = ('k_tmasl', 'k_tmoml', 'k_tnrgl', 'k_xkwlx', 'k_cd',
             'emulatorMean',
             'emulatorPrecision',
             'emulatorCorrelation1',
             'emulatorCorrelation2',
             'emulatorCorrelation3',
             'emulatorCorrelation4',
             'emulatorCorrelation5',
             'emulatorCorrelation6',
             'discrepancyPrecision',
             'discrepancyCorrelation1')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
for i, name in enumerate(filenames):
    savename = name + '_' + chain + '.pdf'
    if i < 5:
        theta_samples = mins[i] + (samples[:,i] * (maxs[i] - mins[i]))
    else:
        theta_samples = samples[:,i]
    ax.plot(theta_samples)
    fig.savefig(savename)
    ax.cla()
