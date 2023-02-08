import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import re

plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 20

fig = plt.figure()  
plt.xlim(0, 10)     
plt.ylim(1e-4, 1)
plt.yscale("log")
plt.xlabel("time") 
plt.ylabel(r"$\delta v_x$") 
plt.ylim(1e-4,1e-1)
plt.minorticks_on()

rho1 = 1.0
rho2 = 2.0
grav = 0.1
kwave = 2*np.pi/0.5

Gam = np.sqrt( grav*kwave*(rho2 - rho1)/(rho2 + rho1) )


data = np.loadtxt("vx_evo.dat")
plt.plot(data[:,0], data[:,1], linewidth=4,label=r"$B_0=0$")
plt.plot(data[:,0], 2.5e-4*np.exp( Gam*data[:,0] ), '--', color="black")

B0 = 4.4603102903819275E-002
Gam = np.sqrt( grav*kwave*( (rho2 - rho1)/(rho2 + rho1) - 2*B0**2*kwave/(grav*(rho1 + rho2)) ) )

data = np.loadtxt("vx_evo_B0.5_hdc.dat")
plt.plot(data[:,0], data[:,1], linewidth=4,label=r"$B_0=0.5 B_\mathrm{crit}$")
plt.plot(data[:,0], 3e-4*np.exp( Gam*data[:,0] ), '--', color="black")

data = np.loadtxt("vx_evo_B1.0_hdc.dat")
plt.plot(data[:,0], data[:,1], linewidth=4,label=r"$B_0=B_\mathrm{crit}$")

plt.legend()

plt.savefig("vx_evo_rt.pdf",bbox_inches="tight", pat_inches=1.0,dpi=1000)
#ani.save(fname_anime, writer="imagemagick")
plt.show()
