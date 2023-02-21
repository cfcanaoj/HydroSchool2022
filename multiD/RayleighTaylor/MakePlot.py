import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import re

dirname = sys.argv[1]
step_s = int(sys.argv[2])
step_e = int(sys.argv[3])

def makedirs(path):
    if not os.path.isdir(path):
        os.makedirs(path)

fig = plt.figure()  
plt.xlim(0, 1)     
plt.ylim(0, 1)
plt.xlabel("x axis") 
plt.ylabel("y axis") 

xmin = -0.25
xmax =  0.25
ymin = -0.75
ymax =  0.75


graph_list = [] 
for istep in range(step_s,step_e+1):
    foutname = dirname + "/snap%05d.dat"%(istep) 
    print("making plot ",foutname)
    with open(foutname, 'r') as data_file:
        line = data_file.readline();
        attributes1 = re.findall("\d+\.\d+", line)

        line = data_file.readline();
        attributes2 = re.findall("\d+", line)

    time = float(attributes1[0]) 
    nx = int(attributes2[0])
    ny = int(attributes2[1])

    data = np.loadtxt(foutname)

    x = data[:,0].reshape(ny,nx)
    y = data[:,1].reshape(ny,nx)
    den = data[:,2].reshape(ny,nx)
    vx = data[:,3].reshape(ny,nx)
    vy = data[:,4].reshape(ny,nx)
    pre = data[:,5].reshape(ny,nx)

    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)

    pg00 = plt.text(0.5*(xmin+xmax),ymax*1.1,r"$\mathrm{time}=%.2f$"%(time),horizontalalignment="center")

    im=plt.imshow(den[:,:],extent=(xmin,xmax,ymin,ymax),origin="lower",vmin=1,vmax=2)

    if istep == step_s: 
        plt.colorbar(im,orientation="vertical")
    graph_list.append([pg00,im])               

#    makedirs(dirname + "/pdffile")
#    plt.savefig(dirname + "/pdffile/rt%05d.pdf"%(istep),bbox_inches="tight", pat_inches=1.0,dpi=1000)
    makedirs(dirname + "/pngfile")
    plt.savefig(dirname + "/pngfile/rt%05d.png"%(istep),bbox_inches="tight", pat_inches=1.0,dpi=1000)

#ani = animation.ArtistAnimation(fig, graph_list, interval=200) 
#print("making animation file", fname_anime)
#ani.save(fname_anime, writer="imagemagick")
plt.show()
