import os
import shutil
import numpy as np
import matplotlib as mpl
mpl.use('pdf')
from matplotlib import pyplot as plt
from matplotlib import rc
from scipy.stats import norm
import seaborn as sns

#current_palette = sns.color_palette()
#sns.palplot(current_palette)

rc('text', usetex=True)
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14

# Extract priors
os.system("awk -F 'm_minValues = ' '/m_minValues/{print $2}' outputData/display_env_sub0.txt > step1.txt")
os.system("awk -F ' ,' '/volume/{print $1}' step1.txt > mins.txt")
os.system("awk -F 'm_maxValues = ' '/m_minValues/{print $2}' outputData/display_env_sub0.txt > step2.txt")
os.system("awk -F ' ,' '/volume/{print $1}' step2.txt > maxs.txt")

# Remove the 'Plots' dir if it exists
if (os.path.exists('./Plots')):
	if (os.path.exists('./Plots/truth_data.pdf')):
		os.system("mv Plots/truth_data.pdf .")
	shutil.rmtree('Plots')

f1 = open('outputData/sip_lineSlope_raw_chain.txt','r')
f2 = open('outputData/sip_lineSlope_filtered_chain.txt','r')
f3 = open('mins.txt','r')
f4 = open('maxs.txt','r')


# Read number of chain samples and dimensions

line1 = f1.readline()
p = line1.split()
num_samplesr = int(p[0])
num_dims = int(p[1])
#print num_dims
mcr = []
ccr = []

line1 = f2.readline()
p = line1.split()
num_samplesf = int(p[0])
mcf = []
ccf = []

# Store prior limits and posterior samples as numpy arrays

line1 = f3.readline()
line2 = f4.readline()
p1 = line1.split()
p2 = line2.split()
mins = []
maxs = []

for i in range(1,num_dims+1):
	mins.append(float(p1[i-1]))	
	maxs.append(float(p2[i-1]))
	minsv = np.array(mins)	
	maxsv = np.array(maxs)

for i in range(1,num_samplesr+1):
	datar = f1.readline()
	p = datar.split()
	mcr.append(float(p[0]))
	ccr.append(float(p[1]))	
	mcrv = np.array(mcr)
	ccrv = np.array(ccr)	

for i in range(1,num_samplesf+1):
	dataf = f2.readline()
	p = dataf.split()
	mcf.append(float(p[0]))
	ccf.append(float(p[1]))	
	mcfv = np.array(mcf)
	ccfv = np.array(ccf)	

# MCMC CHAINS

fig = plt.figure()
ax = fig.add_subplot(221)
ax.plot(mcrv,'b')
ax.set_xlabel(r'$\mathrm{Samples}$',fontsize=12)
ax.set_ylabel(r'$\mathrm{Slope,~m}$',fontsize=12)
ax.set_title(r'$\mathrm{Raw~Chain}$', fontsize=12)
ax.set_xlim([0, num_samplesr])
ax.set_ylim([minsv[0], maxsv[0]])
x0,x1 = ax.get_xlim()
y0,y1 = ax.get_ylim()
ax.set_aspect((x1-x0)/(y1-y0))

ax = fig.add_subplot(222)
ax.plot(mcfv,'r')
ax.set_xlabel(r'$\mathrm{Samples}$',fontsize=12)
ax.set_ylabel(r'$\mathrm{Slope,~m}$',fontsize=12)
ax.set_xlim([0, num_samplesf])
ax.set_ylim([minsv[0], maxsv[0]])
ax.set_title(r'$\mathrm{Filtered~Chain}$', fontsize=12)
x0,x1 = ax.get_xlim()
y0,y1 = ax.get_ylim()
ax.set_aspect((x1-x0)/(y1-y0))

ax = fig.add_subplot(223)
ax.plot(ccrv,'b')
ax.set_xlabel(r'$\mathrm{Samples}$',fontsize=12)
ax.set_ylabel(r'$\mathrm{y-intercept,~c}$',fontsize=12)
ax.set_xlim([0, num_samplesr])
ax.set_ylim([minsv[1], maxsv[1]])
ax.set_title(r'$\mathrm{Raw~Chain}$', fontsize=12)
x0,x1 = ax.get_xlim()
y0,y1 = ax.get_ylim()
ax.set_aspect((x1-x0)/(y1-y0))

ax = fig.add_subplot(224)
ax.plot(ccfv,'r')
ax.set_xlabel(r'$\mathrm{Samples}$',fontsize=12)
ax.set_ylabel(r'$\mathrm{y-intercept,~c}$',fontsize=12)
ax.set_xlim([0, num_samplesf])
ax.set_ylim([minsv[1], maxsv[1]])
ax.set_title(r'$\mathrm{Filtered~Chain}$', fontsize=12)
x0,x1 = ax.get_xlim()
y0,y1 = ax.get_ylim()
ax.set_aspect((x1-x0)/(y1-y0))

plt.tight_layout()
fig.savefig('mcmc_chains.pdf')

# HISTOGRAM PLOTS

fig = plt.figure()
ax = fig.add_subplot(221)
plt.hist(mcrv, edgecolor='none')
ax.set_ylabel(r'$\mathrm{Frequency}$',fontsize=12)
ax.set_xlabel(r'$\mathrm{Slope,~m}$',fontsize=12)
ax.set_title(r'$\mathrm{Raw~Chain}$', fontsize=12)
ax.set_xlim([minsv[0], maxsv[0]])
x0,x1 = ax.get_xlim()
y0,y1 = ax.get_ylim()
ax.set_aspect((x1-x0)/(y1-y0))

ax = fig.add_subplot(222)
plt.hist(mcfv,facecolor='r', edgecolor='none')
ax.set_ylabel(r'$\mathrm{Frequency}$',fontsize=12)
ax.set_xlabel(r'$\mathrm{Slope,~m}$',fontsize=12)
ax.set_title(r'$\mathrm{Filtered~Chain}$', fontsize=12)
ax.set_xlim([minsv[0], maxsv[0]])
x0,x1 = ax.get_xlim()
y0,y1 = ax.get_ylim()
ax.set_aspect((x1-x0)/(y1-y0))

ax = fig.add_subplot(223)
plt.hist(ccrv, edgecolor='none')
ax.set_ylabel(r'$\mathrm{Frequency}$',fontsize=12)
ax.set_xlabel(r'$\mathrm{y-intercept,~c}$',fontsize=12)
ax.set_title(r'$\mathrm{Raw~Chain}$', fontsize=12)
ax.set_xlim([minsv[1], maxsv[1]])
x0,x1 = ax.get_xlim()
y0,y1 = ax.get_ylim()
ax.set_aspect((x1-x0)/(y1-y0))

ax = fig.add_subplot(224)
plt.hist(ccfv,facecolor='r', edgecolor='none')
ax.set_ylabel(r'$\mathrm{Frequency}$',fontsize=12)
ax.set_xlabel(r'$\mathrm{y-intercept,~c}$',fontsize=12)
ax.set_title(r'$\mathrm{Filtered~Chain}$', fontsize=12)
ax.set_xlim([minsv[1], maxsv[1]])
x0,x1 = ax.get_xlim()
y0,y1 = ax.get_ylim()
ax.set_aspect((x1-x0)/(y1-y0))

plt.tight_layout()
fig.savefig('histogram_plots.pdf')

# KDE plots

fig = plt.figure()
ax = fig.add_subplot(221)
ax = sns.distplot(mcrv, fit=norm, kde=False, hist=False, color="b")
ax.set_xlabel(r'$\mathrm{Slope,~m}$',fontsize=12)
ax.set_title(r'$\mathrm{Raw~Chain}$', fontsize=12)
#ax.set_xlim([minsv[0], maxsv[0]])
x0,x1 = ax.get_xlim()
y0,y1 = ax.get_ylim()
ax.set_aspect((x1-x0)/(y1-y0))

ax = fig.add_subplot(222)
ax = sns.distplot(mcfv, fit=norm, kde=False, hist=False, color="r")
ax.set_xlabel(r'$\mathrm{Slope,~m}$',fontsize=12)
ax.set_title(r'$\mathrm{Filtered~Chain}$', fontsize=12)
#ax.set_xlim([minsv[0], maxsv[0]])
x0,x1 = ax.get_xlim()
y0,y1 = ax.get_ylim()
ax.set_aspect((x1-x0)/(y1-y0))

ax = fig.add_subplot(223)
ax = sns.distplot(ccrv, fit=norm, kde=False, hist=False, color="b")
ax.set_xlabel(r'$\mathrm{y-intercept,~c}$',fontsize=12)
ax.set_title(r'$\mathrm{Raw~Chain}$', fontsize=12)
#ax.set_xlim([minsv[0], maxsv[0]])
x0,x1 = ax.get_xlim()
y0,y1 = ax.get_ylim()
ax.set_aspect((x1-x0)/(y1-y0))

ax = fig.add_subplot(224)
ax = sns.distplot(ccfv, fit=norm, kde=False, hist=False, color="r")
ax.set_xlabel(r'$\mathrm{y-intercept,~c}$',fontsize=12)
ax.set_title(r'$\mathrm{Filtered~Chain}$', fontsize=12)
#ax.set_xlim([minsv[0], maxsv[0]])
x0,x1 = ax.get_xlim()
y0,y1 = ax.get_ylim()
ax.set_aspect((x1-x0)/(y1-y0))

plt.tight_layout()
fig.savefig('kde_plots.pdf')


# JOINT PDF

fig = plt.figure()
ax = fig.add_subplot(111)
ax = sns.kdeplot(mcrv, ccrv, shade=True)
ax.set_xlabel(r'$\mathrm{Slope,~m}$',fontsize=12)
ax.set_ylabel(r'$\mathrm{y-intercept,~c}$', fontsize=12)
#ax.set_xlim([minsv[0], maxsv[0]])
x0,x1 = ax.get_xlim()
y0,y1 = ax.get_ylim()
ax.set_aspect((x1-x0)/(y1-y0))
fig.savefig('joint_pdf')

## SFP PLOTS

# Read number of chain samples and dimensions

f1 = open('outputData/sfp_lineSlope_qoi_seq_post.txt')
line1 = f1.readline()
p = line1.split()
num_samples = int(p[0])
num_QoIs = int(p[1])
#print num_dims
qoi = []

# Construct the QoI vector

for i in range(1,num_samples+1):
	data = f1.readline()
	p = data.split()
	qoi.append(float(p[0]))
	qoiv = np.array(qoi)

# QoI MC CHAIN

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(qoiv,'b')
ax.set_xlabel(r'$\mathrm{Samples}$',fontsize=12)
ax.set_ylabel(r'$\mathrm{Qoi~(y~at~x = 3.0)}$',fontsize=12)
ax.set_title(r'$\mathrm{MC~Chain}$', fontsize=12)
ax.set_xlim([0, num_samples])
#ax.set_ylim([minsv[0], maxsv[0]])
x0,x1 = ax.get_xlim()
y0,y1 = ax.get_ylim()
ax.set_aspect((x1-x0)/(y1-y0))
fig.savefig('qoi_mcChain')

# KDE plot

fig = plt.figure()
ax = fig.add_subplot(111)
ax = sns.distplot(qoiv, fit=norm, kde=False, hist=False, color="b")
ax.set_xlabel(r'$\mathrm{Qoi~(y~at~x = 3.0)}$',fontsize=12)
#ax.set_xlim([minsv[0], maxsv[0]])
x0,x1 = ax.get_xlim()
y0,y1 = ax.get_ylim()
ax.set_aspect((x1-x0)/(y1-y0))
fig.savefig('qoi_pdf')

# Sensitivity Plot

f1 = open('outputData/sense_mc.txt','r')
line1 = f1.readline()
p = line1.split()
num_samples = int(p[0])
dims = int(p[1])
qoi_mc = []

for i in range(1,num_samples+1):
	data1 = f1.readline()
	p = data1.split()
	qoi_mc.append(float(p[0]))
	qoi_mcv = np.array(qoi_mc)

f1.close()


f2 = open('outputData/sense_m.txt','r')
data2 = f2.readlines()[1:]
f2.close()
qoi_m = []

for line in data2:
	p = line.split()
	qoi_m.append(float(p[0]))
	qoi_mv = np.array(qoi_m)

#np.savetxt('qoi_mv.txt',qoi_mv)


f3 = open('outputData/sense_c.txt','r')
data3 = f3.readlines()[1:]
f3.close()
qoi_c = []

for line in data3:
	p = line.split()
	qoi_c.append(float(p[0]))
	qoi_cv = np.array(qoi_c)

m1 = (np.var(qoi_mv))/(np.var(qoi_mcv))
c1 = float(np.var(qoi_cv))/float(np.var(qoi_mcv))

print m1, c1, m1+c1

### Bar-graph for sensitivity

index = np.arange(2)
s = np.array([m1, c1])
bar_width = 0.07

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim([0,1])
ax.set_ylim([0,1])
plt.bar(0.25, m1, bar_width, color = 'b', edgecolor = 'none')
plt.bar(0.75, c1, bar_width, color = 'r', edgecolor = 'none')
ax.set_ylabel(r'$\mathrm{First~Order~Sensitivity}$', fontsize=14)
ax.set_title(r'$\mathrm{Sensitivity~Analysis}$', fontsize=14)
plt.xticks([0.25+bar_width/2, 0.75+bar_width/2], (r'$\mathrm{m}$', r'$\mathrm{c}$'))
fig.savefig('sensitivity_plot') 

source = './'
os.makedirs('Plots')
dest = './Plots'
files = os.listdir(source)

for f in files:
	if (f.endswith('pdf')):
		shutil.move(f,dest)

os.system("rm step1.txt step2.txt")


