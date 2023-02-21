import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import re

fig = plt.figure()  
plt.xlim(0, 1)     
plt.ylim(0, 1)
plt.xlabel("x axis") 
plt.ylabel("y axis") 

xmin = -0.5
xmax =  0.5
ymin = -0.5
ymax =  0.5

dirname = "hlld"
base = "snap"
suffix = ".dat"
fname_anime = "animation.mp4"

step_s = 0
step_e = 50

graph_list = [] 
for istep in range(step_s,step_e+1):
    foutname = dirname + "/" + base + "%05d"%(istep) + suffix
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
    pre = data[:,5].reshape(ny,nx)

    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)

    pg00 = plt.text(0.5*(xmin+xmax),ymax*1.1,r"$\mathrm{time}=%.2f$"%(time),horizontalalignment="center")

    dx = x[0,1] - x[0,0]
    dy = y[1,0] - y[0,0]

    dvxdx = np.gradient( vx, axis=1 )/dx
    dvydy = np.gradient( vy, axis=0 )/dy
    dvxdy = np.gradient( vx, axis=0 )/dy
    dvydx = np.gradient( vy, axis=1 )/dx

#    for j in range(1,ny-1):
#        for i in range(1,nx-1):
#            dvxdx[j,i] = 0.5*( vx[j,i+1] - vx[j,i-1] )/dx
#            dvydy[j,i] = 0.5*( vy[j+1,i] - vy[j-1,i] )/dy

    omega = dvydx - dvxdy
#    omega = dvxdx + dvydy

    ome_disp = np.sqrt(np.average(omega[1:nx-1,1:nx-1]**2))
    print(ome_disp, np.sqrt(np.average((dvxdx+dvydy)[1:nx-1,1:nx-1]**2)))


    im=plt.imshow(omega/ome_disp,extent=(xmin,xmax,ymin,ymax),origin="lower",cmap=plt.cm.bwr,vmin=-5,vmax=5)

    if istep == step_s: 
        plt.colorbar(im,orientation="vertical")
    graph_list.append([pg00,im])               


ani = animation.ArtistAnimation(fig, graph_list, interval=200) 
print("making animation file", fname_anime)
ani.save(dirname + "/" + fname_anime, writer="imagemagick")
#plt.show()
