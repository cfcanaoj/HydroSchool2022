import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import re

plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 30


fig = plt.figure(figsize=(14, 8)) 
plt.xlim(0, 1)     
plt.ylim(-1.1, 1.1)
plt.xlabel(r"$x$") 
plt.minorticks_on()


#argsor = np.argsort(data[:,0])
#print(data[argsor,0])

plt.plot([0.32957362497368775,0.32957362497368775],[-1.1,1.1],'--',color="black")
plt.plot([0.473241830023867,0.473241830023867],[-1.1,1.1],'--',color="black")
plt.plot([0.49253975058004595,0.49253975058004595],[-1.1,1.1],'--',color="black")
plt.plot([0.564841,0.564841],[-1.1,1.1],'--',color="black")
plt.plot([0.6325009038473688,0.6325009038473688], [-1.1,1.1],'--',color="black")
plt.plot([0.86557591513473,0.86557591513473], [-1.1,1.1],'--',color="black")

#for i, size in enumerate(size_list):
#plt.text(i*0.1, i*0.1, f"{size} pt", size =size  ) #  fontsize or size = size
plt.text(0.32957362497368775, 1.15, f"FR", fontsize=20,horizontalalignment="center") #  fontsize or size = size
plt.text(0.473241830023867*0.98, 1.15, f"IS", fontsize=20,horizontalalignment="center") #  fontsize or size = size
plt.text(0.49253975058004595*1.02, 1.15, f"SR", fontsize=20,horizontalalignment="center") #  fontsize or size = size
plt.text(0.564841, 1.15, f"CD", fontsize=20,horizontalalignment="center") #  fontsize or size = size
plt.text(0.6325009038473688, 1.15, f"SS", fontsize=20,horizontalalignment="center") #  fontsize or size = size
plt.text(0.86557591513473, 1.15, f"FR", fontsize=20,horizontalalignment="center") #  fontsize or size = size
 

data = np.loadtxt("briowu_nonregsol.dat")
plt.plot(data[:,0],data[:,1],linewidth=3,label=r"$\rho$")
plt.plot(data[:,0],data[:,2],linewidth=3,label=r"$P$")
plt.plot(data[:,0],data[:,3],linewidth=3,label=r"$v_x$")
plt.plot(data[:,0],data[:,7],linewidth=3,label=r"$B_y$")

#fout = open("l.dat","w")
#for i in range(data[:,0].shape[0]):
#    fout.write("{0} ".format(data[argsor[i],0]+0.5))
#    for k in range(1,10): 
#        fout.write("{0} ".format(data[argsor[i],k]))
#    fout.write("\n")
#fout.close()

    


plt.legend(fontsize=20)

plt.savefig("briowu_nonregsol.pdf",bbox_inches="tight")
#ani.save(fname_anime, writer="imagemagick")
plt.show()
