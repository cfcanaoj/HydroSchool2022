import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import re

dirname = sys.argv[1]
step_s = int(sys.argv[2])

def makedirs(path):
    if not os.path.isdir(path):
        os.makedirs(path)

fig = plt.figure()  
plt.xlim(0, 1)     
plt.ylim(0, 1)
plt.xlabel("x axis") 
plt.ylabel("y axis") 

xmin = -0.5
xmax =  0.5
ymin = -1.0
ymax =  1.0

for istep in range(step_s,step_s+1):
    foutname = dirname + "/snap%05d.bin"%(istep) 
    fp = open(foutname,"rb")

    time = np.fromfile(fp,np.float64,1)[0]
    nx   = np.fromfile(fp,np.int32,1)[0]
    ny   = np.fromfile(fp,np.int32,1)[0]
    nvar = np.fromfile(fp,np.int32,1)[0]
    xv   = np.fromfile(fp,np.float64,nx)
    yv   = np.fromfile(fp,np.float64,ny)
    Q    = np.fromfile(fp,np.float32,nx*ny*nvar).reshape(ny,nx,nvar)

    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)

    plt.text(0.5*(xmin+xmax),ymax*1.1,r"$\mathrm{time}=%.2f$"%(time),horizontalalignment="center")

    im=plt.imshow(Q[:,:,9],extent=(xmin,xmax,ymin,ymax),origin="lower",vmin=0,vmax=1)

    plt.colorbar(im,orientation="vertical")

    makedirs(dirname + "/pdffile")
    plt.savefig(dirname + "/pdffile/rt%05d.pdf"%(istep),bbox_inches="tight", pat_inches=1.0,dpi=1000)

#ani = animation.ArtistAnimation(fig, graph_list, interval=200) 
#print("making animation file", fname_anime)
#ani.save(fname_anime, writer="imagemagick")
plt.show()
