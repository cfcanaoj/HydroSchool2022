import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import ArtistAnimation
from matplotlib import gridspec
import re

plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 30

x = np.linspace(0, np.pi * 4, 100)

fig = plt.figure(figsize=(10,8)) 
gs = gridspec.GridSpec(1, 1, wspace=0.35,hspace=0.2) 

ax0 = plt.subplot(gs[0])
ax0.minorticks_on()
ax0.set_xlim(1,5000)
ax0.set_xlabel(r"$N$")
ax0.set_ylabel(r"$\epsilon$")
ax0.set_xscale("log")
ax0.set_yscale("log")

data= np.loadtxt("error1st.dat")

    # グラフを作成する。
ax0.plot(data[:,0], data[:,1], 'o-',c="r")
ax0.plot(data[:,0], 1e-3*data[:,0]**(-1), '-',c="k",label=r"$\epsilon\propto N^{-1}$")
ax0.legend()

#ax0 = plt.subplot(gs[1])
#ax0.minorticks_on()
#ax0.set_xlim(0,1.0)
#ax0.set_xlabel(r"$t/t_0$")
#ax0.set_ylabel(r"$\rho(x=0)$")
#
#filename = ["sound8.dat","sound32.dat","sound128.da"]
#
#for num in [8,32,128]: 
#    data= np.loadtxt("sound%d.dat"%(num))
#
#    # グラフを作成する。
#    ax0.plot(data[:,0]*np.sqrt(5.0/3.0), data[:,1], '-',label=r"$N=%d$"%(num),linewidth=3)
#
#x = np.linspace(0,1.0,100)
#ax0.plot(x, 1.0 + 1e-4*np.sin(2.0*np.pi*x), '-',c="k")
#
#ax0.legend()



# mp4 画像として保存する。
plt.savefig("sound_err.pdf",bbox_inches="tight")
plt.show()
#plt.close()
