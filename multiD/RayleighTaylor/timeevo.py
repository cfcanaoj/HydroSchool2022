import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import re

fig = plt.figure()  
plt.xlim(0, 10)     
plt.ylim(1e-4, 1)
plt.yscale("log")
plt.xlabel("time") 
plt.ylabel(r"$\delta v_x$") 
plt.minorticks_on()

rho1 = 1.0
rho2 = 2.0
grav = 0.1
kwave = 2*np.pi/0.5

Gam = np.sqrt( grav*kwave*(rho2 - rho1)/(rho2 + rho1) )


data = np.loadtxt("vx_evo.dat")
plt.plot(data[:,0], data[:,1], label="simulation")
plt.plot(data[:,0], 3e-4*np.exp( Gam*data[:,0] ), label="linear analysis")


plt.legend()

plt.savefig("vx_evo.pdf",bbox_inches="tight", pat_inches=1.0,dpi=1000)
#ani.save(fname_anime, writer="imagemagick")
plt.show()
