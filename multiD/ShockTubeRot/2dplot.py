import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import re

fig = plt.figure()  
plt.xlim(0, 1)     
plt.ylim(0, 1)
plt.xlabel("x axis") 
plt.ylabel("y axis") 

graph_list = [] 
for istep in range(2,3):
#    foutname = "snap/snap%05d_uw.dat"%(istep)
    foutname = "snap/snap%05d_256.dat"%(istep)
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
    nx = 128*2
    ny = 128*2

    x = data[:,0].reshape(ny,nx)
    y = data[:,1].reshape(ny,nx)
    den = data[:,2].reshape(ny,nx)
    vx = data[:,3].reshape(ny,nx)
    vy = data[:,4].reshape(ny,nx)
    pre = data[:,5].reshape(ny,nx)

    print(x.min(),x.max())

    plt.xlim(-0.5,0.5)
    plt.ylim(-0.5,0.5)

    im=plt.imshow(den[:,:],extent=(-0.5,0.5,-0.5,0.5),origin="lower",vmin=1,vmax=3)
#    im=plt.imshow(den[:,:],extent=(-1,1,-1,1),origin="lower")

    plt.colorbar(im)

plt.show()
#    data = np.loadtxt(foutname)
#    print data
#    xv = data[:,0]
#    U  = data[:,1]


#    graph = plt.plot(xv, U, 'o-', color="black") 
#    graph_list.append(graph)               

#ani = animation.ArtistAnimation(fig, graph_list, interval=200) 
#plt.show()
