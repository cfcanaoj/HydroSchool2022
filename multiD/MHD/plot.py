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
icount = 0
#for istep in range(1,101):
for istep in range(100,101):
    foutname = "output/bin%05d_uw.dat"%(istep)
    print(foutname)
    with open(foutname, 'r') as data_file:
        line = data_file.readline();
        attributes1 = re.search(r'time = (\S+)', line)
    print(attributes1)
#        attributes2 = re.search(r'nx=         (\S+)', line)
#        attributes3 = re.search(r'ny=         (\S+)', line)
#    time = float(attributes1.group(1)) 
#    nx = int(attributes2.group(1)) 
#    ny = int(attributes3.group(1)) 

    data = np.loadtxt(foutname)
    nx = int(np.sqrt(data[:,0].shape[0]))
    ny = nx

    den = data[:,2].reshape(ny,nx)
    vx = data[:,3].reshape(ny,nx)
    vy = data[:,4].reshape(ny,nx)
    pre = data[:,5].reshape(ny,nx)


    im=plt.imshow(den[:,:],extent=(0,1,0,1),origin="lower",vmin=0.5,vmax=2.4,animated=True)

    if icount==0:
        plt.colorbar(im)
    graph_list.append([im])
    icount+=1

#ani = animation.ArtistAnimation(fig, graph_list, interval=200) 
plt.show()
#plt.show()
#    data = np.loadtxt(foutname)
#    print data
#    xv = data[:,0]
#    U  = data[:,1]


#    graph = plt.plot(xv, U, 'o-', color="black") 
#    graph_list.append(graph)               

#ani = animation.ArtistAnimation(fig, graph_list, interval=200) 
#plt.show()
