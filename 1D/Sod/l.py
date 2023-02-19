import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import ArtistAnimation
from matplotlib import gridspec
import re

plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 20

dirname = "lax"
nmax = 40


x = np.linspace(0, np.pi * 4, 100)


frames = []  

anasol = np.loadtxt("sod_ana.dat")
tout_ana = 0.2
x_ana = anasol[:,0]
den_ana = anasol[:,1]
vel_ana = anasol[:,2]
pre_ana = anasol[:,3]

icount = 0
for istep in range(0,nmax+1):
    fig = plt.figure(figsize=(8,10)) 
    gs = gridspec.GridSpec(3, 1, wspace=0.30,hspace=0.2) 
    ax = []
    for i in range(3):
        ax.append( plt.subplot(gs[i]) )
    
    ylabel = [r"$\rho$",r"$v$",r"$P$"]
    for i,ax0 in enumerate(ax):
        ax0.minorticks_on()
        ax0.set_xlim(-0.5,0.5)
        ax0.set_ylabel(ylabel[i])
    
    ax[2].set_xlabel(r"$x$")

    foutname = dirname + "/snap" + "%05d.dat"%(istep)
    print("reading " + foutname)
    with open(foutname, 'r') as data_file: 
        line = data_file.readline() 
        attributes = re.search(r'# (\S+)', line)
    time = float(attributes.group(1))
    print(time)

    data = np.loadtxt(foutname)

    x = data[:,0]
    den = data[:,1]
    vel = data[:,2]
    pre = data[:,3]

    
    ax[0].plot(x, den, 'o-',c="r",label="numerical")
    ax[0].plot(x_ana*time/tout_ana, den_ana, '-',c="b",label="exact")
    ax[1].plot(x, vel, 'o-',c="r",label="numerical")
    ax[1].plot(x_ana*time/tout_ana, vel_ana, '-',c="b",label="exact")
    ax[2].plot(x, pre, 'o-',c="r",label="numerical")
    ax[2].plot(x_ana*time/tout_ana, pre_ana, '-',c="b",label="exact")

    ax[0].text(0,1.10,r"$\mathrm{time} = %.2f$"%(time),horizontalalignment="center")

    ax[0].legend()

    plt.savefig(dirname + "/snap%05d.pdf"%(istep),bbox_inches="tight")
    icount +=1
    plt.clf()
    plt.close()

#ani = ArtistAnimation(fig, frames, interval=50)
#
#ani.save(dirname + "/pyanime.mp4", writer="imagemagick")
#plt.show()
#plt.close()
