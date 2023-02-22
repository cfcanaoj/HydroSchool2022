import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import re

dirname = sys.argv[1]
step_s = int(sys.argv[2])
step_e = int(sys.argv[3])

fig = plt.figure()  
plt.xlim(0, 1)     
plt.ylim(0, 1)
plt.xlabel("x axis") 
plt.ylabel("y axis") 

xmin = -0.5
xmax =  0.5
ymin = -0.5
ymax =  0.5

fname_anime = "animation.mp4"


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

    print(data[:,0].shape)

    x = data[:,0].reshape(ny,nx)
    y = data[:,1].reshape(ny,nx)
    den = data[:,2].reshape(ny,nx)
    vx = data[:,3].reshape(ny,nx)
    vy = data[:,4].reshape(ny,nx)
    vz = data[:,5].reshape(ny,nx)
    pre = data[:,6].reshape(ny,nx)
    bx = data[:,7].reshape(ny,nx)
    by = data[:,8].reshape(ny,nx)
    bz = data[:,9].reshape(ny,nx)
    vor = data[:,10].reshape(ny,nx)
    Bpot = data[:,11].reshape(ny,nx)

    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)

    pg00 = plt.text(0.5*(xmin+xmax),ymax*1.1,r"$\mathrm{time}=%.2f$"%(time),horizontalalignment="center")

    vor_disp = np.sqrt(np.average(vor[1:nx-1,1:nx-1]**2))


    im=plt.imshow(vor/vor_disp,extent=(xmin,xmax,ymin,ymax),origin="lower",cmap=plt.cm.bwr,vmin=-5,vmax=5)

    if istep == step_s: 
        plt.colorbar(im,orientation="vertical")
    graph_list.append([pg00,im])               


ani = animation.ArtistAnimation(fig, graph_list, interval=200) 
print("making animation file", fname_anime)
ani.save(dirname + "/" + fname_anime, writer="imagemagick")
#plt.show()
