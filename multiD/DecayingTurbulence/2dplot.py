import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import re

fig = plt.figure()  

xmin = -0.5; xmax = 0.5
ymin = -0.5; ymax = 0.5


plt.xlim(xmin, xmax)     
plt.ylim(ymin, ymax)
plt.xlabel("x axis") 
plt.ylabel("y axis") 


graph_list = [] 
for istep in range(0,1):
#    foutname = "snap/snap%05d_uw.dat"%(istep)
    foutname = "openmp_test_4cpu/kh%05d.dat"%(istep)
#    foutname = "l"
#    with open(foutname, 'r') as data_file:
#        line = data_file.readline();
#        attributes1 = re.search(r'time=   (\S+)', line)
#        attributes2 = re.search(r'nx=         (\S+)', line)
#        attributes3 = re.search(r'ny=         (\S+)', line)
#    time = float(attributes1.group(1)) 
#    nx = int(attributes2.group(1)) 
#    ny = int(attributes3.group(1)) 

    data = np.loadtxt(foutname)
    nx = 128*1
    ny = 128*1

    x = data[:,0].reshape(ny,nx)
    y = data[:,1].reshape(ny,nx)
    den = data[:,2].reshape(ny,nx)
    vx = data[:,3].reshape(ny,nx)
    vy = data[:,4].reshape(ny,nx)
    pre = data[:,5].reshape(ny,nx)

    print(x.min(),x.max())

    plt.xlim(-0.5,0.5)
    plt.ylim(-0.5,0.5)

    dx = x[0,1] - x[0,0]
    dy = y[1,0] - y[0,0]

    dvxdy = np.gradient( vx, axis=1 )/dy
    dvydx = np.gradient( vy, axis=0 )/dx

    omega = dvydx - dvxdy

#    im=plt.imshow(omega[:,:],extent=(xmin,xmax,ymin,ymax),origin="lower",cmap=plt.cm.bwr,vmin=-1.0,vmax=1.0)
    im=plt.imshow(omega[:,:],extent=(xmin,xmax,ymin,ymax),origin="lower",cmap=plt.cm.bwr)

    cb = plt.colorbar(im)
    cb.set_label('vorticity')

plt.savefig("vorticity.pdf",bbox_inches="tight")
plt.show()
#    data = np.loadtxt(foutname)
#    print data
#    xv = data[:,0]
#    U  = data[:,1]


#    graph = plt.plot(xv, U, 'o-', color="black") 
#    graph_list.append(graph)               

#ani = animation.ArtistAnimation(fig, graph_list, interval=200) 
#plt.show()
