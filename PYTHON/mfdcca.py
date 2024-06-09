import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
import time
import sys
def printf(format, *args):
    sys.stdout.write(format % args)

import numpy as np
import mfdcca

MAX_BOX=200		# maximum 200 points on logarithmic scale...
MAXQ=201		# max q resolution -10,...,10 dq=0.1

# Load data.
# os.chdir("D:/laptop_07_02_2010/D/papers/MDFA/MFDCCALIB/PYTHON")

series = 'BIOMIAL_MULTIFRACTAL_65536_0.3_0.4'	# binomial multifractal model
#series = 'sug_eta_returns'			#sugar/etanol returns

data1 = np.loadtxt('../data/'+series+'.dat',skiprows=1, usecols=0)
data2 = np.loadtxt('../data/'+series+'.dat',skiprows=1, usecols=1)
total = data1.shape[0]
printf("Loaded two columnns of %d data points from %s\n", total, '../data/'+series+'.dat')

# Parameters.
dq = 0.1
qmin = -10.000001	# small offset to avoid 0
qmax = 10
n = int(np.abs(qmax - qmin)/dq)+1

# arrays for results
H = np.empty(n)
tau = np.empty(n)
f = np.empty(n)
alpha = np.empty(n)
dmse = np.ndarray(shape = [MAXQ, MAX_BOX])	# maximum 200 points on logarithmic scale, and 201 maximum q values
nfit=1		# polynomial degree 1 (linear), 2 or 3
sw=0		# sw=1 for sliding segments, takes longer, use with care...

rs=np.zeros(MAX_BOX, dtype=int)			# segment size list
# prepare multiplicative scale (equidistant on log-log plot) using library
nrs=0
minseg=4		# minimum segment size
maxseg=int(total/4)	# maximum segment size
boxratio=np.power(2,1.0/8)	# segment ratio (multiplicative scale)
#...OR prepare your own segment scale
#rs=np.arange(10, 101, 1) 	
#nrs=len(rs)

# version of MFDCCA algorithm:
# 1 - Original MF-DXA W.-X. Zhou, Phys. Rev. E 77, 066211 (2008).
# 2 - ABS
# 3 - MFCCA Phys. Rev. E 89, 023305 (2014)
# 4 - Plus sum
# 5 - Minus sum
# 6 - Plus Box
# 7 - Minus box
# 8 - PP
# 9 - PM
# 10 - MP
# 11 - MM
ver=1

#run the MFDCCA C dynamic link library code
t0 = time.time()
perc=mfdcca.calc(data1,data2,total,ver,H,tau,f,alpha,dmse,rs,nrs,qmin,qmax,dq,minseg,maxseg,boxratio,nfit,sw)
t1 = time.time()
runtime = t1-t0
printf("MFDCCA for q from %.2f to %.2f with step %.2f runtime=%.4f sec\n", qmin, qmax, dq, runtime)
printf("used %.2f%% data pairs",perc)

nr=np.count_nonzero(rs)
# Do the plots
fig = plt.figure(tight_layout=True)
gs = gridspec.GridSpec(2, 2)

qlist = []
for q in np.arange(qmin,qmax,dq):
    qlist.append(q)
q = np.asarray(qlist)

ax = fig.add_subplot(gs[0, 0])		# plot tau
ax.plot(q,tau, marker='o')
ax.set_ylabel('tau(q)')
ax.set_xlabel('q')

ax = fig.add_subplot(gs[0, 1])		# plot fluctuations
for i in range(0, n):
    ax.loglog(rs[0:nr], dmse[i,0:nr], 'o', label = q[i])
ax.set_ylabel('F(s)')
ax.set_xlabel('s')

ax = fig.add_subplot(gs[1, 0])		# plot f(alpha)
ax.scatter(alpha, f, s=40, facecolors='none', edgecolors='b')
ax.set_ylabel('f(alpha)')
ax.set_xlabel('alpha')

ax = fig.add_subplot(gs[1, 1])		#plot H(q)
ax.plot(q,H, marker='o')
ax.set_ylabel('H(q)')
ax.set_xlabel('q')

plt.show()
fig.savefig(series + str(ver)+'.png')

