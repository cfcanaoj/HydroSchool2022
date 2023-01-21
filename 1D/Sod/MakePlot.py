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

fig = plt.figure(figsize=(20,8)) 
gs = gridspec.GridSpec(1, 2, wspace=0.20,hspace=0.2) 

ax0 = plt.subplot(gs[0])
ax0.minorticks_on()
ax0.set_xlim(-0.5,0.5)
ax0.set_xlabel(r"$x$")
ax0.set_title("1st order")


anasol = np.loadtxt("sod_ana.dat")
tout_ana = 0.2
x_ana = anasol[:,0]
den_ana = anasol[:,1]
vel_ana = anasol[:,2]
pre_ana = anasol[:,3]

istep = 40
time = 0.2

data_lax = np.loadtxt("snap/lax%05d.dat"%(istep))
data_hll = np.loadtxt("snap/hll%05d.dat"%(istep))
data_hllc = np.loadtxt("snap/hllc%05d.dat"%(istep))

    # グラフを作成する。
ax0.plot(data_lax[:,0], data_lax[:,1], 'o-',c="r",label="Lax")
ax0.plot(data_hll[:,0], data_hll[:,1], 'o-',c="b",label="HLL")
ax0.plot(data_hllc[:,0], data_hllc[:,1], 'o-',c="g",label="HLLC")
ax0.plot(x_ana*time/tout_ana, den_ana, '-',c="k",label="exact")
ax0.legend()

ax0 = plt.subplot(gs[1])
ax0.minorticks_on()
ax0.set_xlim(-0.5,0.5)
ax0.set_xlabel(r"$x$")
ax0.set_title("2nd order")


data_lax = np.loadtxt("snap/lax2nd%05d.dat"%(istep))
data_hll = np.loadtxt("snap/hll2nd%05d.dat"%(istep))
data_hllc = np.loadtxt("snap/hllc2nd%05d.dat"%(istep))

    # グラフを作成する。
ax0.plot(data_lax[:,0], data_lax[:,1], 'o-',c="r",label="Lax")
ax0.plot(data_hll[:,0], data_hll[:,1], 'o-',c="b",label="HLL")
ax0.plot(data_hllc[:,0], data_hllc[:,1], 'o-',c="g",label="HLLC")
ax0.plot(x_ana*time/tout_ana, den_ana, '-',c="k",label="exact")
ax0.legend()

# mp4 画像として保存する。
plt.savefig("sod_lax_hll_hllc.pdf",bbox_inches="tight")
plt.show()
#plt.close()
