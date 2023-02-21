import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
import re

dirname = sys.argv[1]
step = int(sys.argv[2])

x = np.linspace(0, np.pi * 4, 100)

fig = plt.figure(figsize=(10,10)) 
gs = gridspec.GridSpec(2, 2, wspace=0.30,hspace=0.2) 
ax = []
for i in range(4):
    ax.append( plt.subplot(gs[i]) )

ylabel = [r"$\rho$",r"$v_x$",r"$v_y$",r"$B_y$"]
for i,ax0 in enumerate(ax):
    ax0.minorticks_on()
    ax0.set_xlim(-0.5,0.5)
    ax0.set_ylabel(ylabel[i])

ax[2].set_xlabel(r"$x$")
ax[3].set_xlabel(r"$x$")


anasol = np.loadtxt("briowu_nonregsol.dat")
tout_ana = 0.1
x_ana = anasol[:,0] - 0.5
den_ana = anasol[:,1]
vx_ana = anasol[:,3]
vy_ana = anasol[:,4]
by_ana = anasol[:,7]

# フレームごとの Artist を作成する。

icount = 0
foutname = dirname + "/snap%05d.dat"%(step)
print("reading " + foutname)
with open(foutname, 'r') as data_file: 
    line = data_file.readline() 
    attributes = re.search(r'# (\S+)', line) 
    time = float(attributes.group(1))

data = np.loadtxt(foutname)

x = data[:,0]
den = data[:,1]
vx = data[:,2]
vy = data[:,3]
By = data[:,7]

# グラフを作成する。
ax[0].plot(x, den, 'o-',c="r",label="numerical")
ax[0].plot(x_ana*time/tout_ana, den_ana, '-',c="b",label="exact")
ax[1].plot(x, vx, 'o-',c="r",label="numerical")
ax[1].plot(x_ana*time/tout_ana, vx_ana, '-',c="b",label="exact")
ax[2].plot(x, vy, 'o-',c="r",label="numerical")
ax[2].plot(x_ana*time/tout_ana, vy_ana, '-',c="b",label="exact")
ax[3].plot(x, By, 'o-',c="r",label="numerical")
ax[3].plot(x_ana*time/tout_ana, by_ana, '-',c="b",label="exact")

ax[0].text(0,1.10,r"$\mathrm{time} = %.2f$"%(time),horizontalalignment="center")

if icount == 0: 
    ax[0].legend()


icount +=1

#plt.savefig(dirname + '/snap%05d.png'%(step),bbox_inches="tight", pat_inches=0.0,dpi=1000)
plt.savefig(dirname + '/snap%05d.pdf'%(step),bbox_inches="tight")

#plt.show()
#plt.close()
