import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import re
from matplotlib import gridspec

plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 20


fig = plt.figure(figsize=(14, 8)) 
gs = gridspec.GridSpec(2, 2, wspace=0.2,hspace=0.2) 

plt.xlim(0, 1)     
data = np.loadtxt("daiwoodward_exact2.dat")

label = [r"$\rho$", r"$P$", r"$v_y$", r"$B_y$"]
data_sol = [ data[:,1], data[:,2], data[:,4], data[:,7]]
ymin = [0.9, 0.9, -0.45, 0.9]
ymax = [4.2, 170, 0.45, 6]

xpos = [0.11576053157111393,0.4910218753698951,0.5579008,0.6146880167865735,0.86557591513473]

for ipl in range(4): 
    ax0 = plt.subplot(gs[ipl]) 
    ax0.set_xlim(0, 1)
    ax0.set_ylim(ymin[ipl], ymax[ipl])
    ax0.set_xlabel(r"$x$") 
    ax0.minorticks_on()
    ax0.plot(data[:,0],data_sol[ipl],linewidth=3)

    ax0.text(0.05,ymin[ipl] + (ymax[ipl] - ymin[ipl])*0.8, label[ipl])

    for xx in xpos:
        ax0.plot( [xx,xx],[ymin[ipl],ymax[ipl]], '--', color='black', zorder=0)



#plt.plot([0.32957362497368775,0.32957362497368775],[-1.1,1.1],'--',color="black")
#plt.plot([0.473241830023867,0.473241830023867],[-1.1,1.1],'--',color="black")
#plt.plot([0.49253975058004595,0.49253975058004595],[-1.1,1.1],'--',color="black")
#plt.plot([0.564841,0.564841],[-1.1,1.1],'--',color="black")
#plt.plot([0.6325009038473688,0.6325009038473688], [-1.1,1.1],'--',color="black")
#plt.plot([0.86557591513473,0.86557591513473], [-1.1,1.1],'--',color="black")

#size_list = [ 10,20,30,40,50]
#for i, size in enumerate(size_list):
##plt.text(i*0.1, i*0.1, f"{size} pt", size =size  ) #  fontsize or size = size
#plt.text(0.32957362497368775, 1.15, f"FR", fontsize=20,horizontalalignment="center") #  fontsize or size = size
#plt.text(0.473241830023867*0.98, 1.15, f"RD", fontsize=20,horizontalalignment="center") #  fontsize or size = size
#plt.text(0.49253975058004595*1.02, 1.15, f"SS", fontsize=20,horizontalalignment="center") #  fontsize or size = size
#plt.text(0.564841, 1.15, f"CD", fontsize=20,horizontalalignment="center") #  fontsize or size = size
#plt.text(0.6325009038473688, 1.15, f"SS", fontsize=20,horizontalalignment="center") #  fontsize or size = size
#plt.text(0.86557591513473, 1.15, f"FR", fontsize=20,horizontalalignment="center") #  fontsize or size = size
 

#plt.plot(data[:,0],data[:,1],linewidth=3,label=r"$\rho$")
#plt.plot(data[:,0],data[:,2],linewidth=3,label=r"$P$")
#plt.plot(data[:,0],data[:,3],linewidth=3,label=r"$v_x$")
#plt.plot(data[:,0],data[:,7],linewidth=3,label=r"$B_y$")


#plt.legend(fontsize=20)

plt.savefig("daiwoodward2.pdf",bbox_inches="tight")
#ani.save(fname_anime, writer="imagemagick")
plt.show()
