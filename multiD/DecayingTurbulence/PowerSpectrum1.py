import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import re

def PowerSpectrum1( A, kwave, kbins ):


    Ak = np.fft.fft2( A )/(nx)

    Ekink = np.abs(Ak)

    print("normalization ",np.sum(np.abs(Ak))/np.sum(A))

    index = np.digitize(kwave.flat, kbins)

    N = len(kbins)

    Epower = np.zeros(N-1)
    for n in range(1,N):
        Epower[n-1] = np.sum( Ekink.flat[index == n] )

    return Epower


def PowerSpectrum( Ax, Ay, kwave, kbins ):


    Axk = np.fft.fft2( Ax )/(nx)
    Ayk = np.fft.fft2( Ay )/(ny)

    Ekink = 0.5*( np.abs(Axk)**2 + np.abs(Ayk)**2 )

    print("normalization ",np.sum(np.abs(Axk)**2 + np.abs(Ayk)**2)/np.sum(Ax**2 + Ay**2))
    print(np.sum(np.abs(Axk)**2 + np.abs(Ayk)**2))

    index = np.digitize(kwave.flat, kbins)

    N = len(kbins)

    Epower = np.zeros(N-1)
    for n in range(1,N):
        Epower[n-1] = np.sum( Ekink.flat[index == n] )

    return Epower


fig = plt.figure()  
#plt.xlim(0, 1)     
#plt.ylim(0, 1)
plt.xlabel(r"$k/(2\pi)$") 
plt.ylabel("power spectrum") 

xmin = -0.5
xmax =  0.5
ymin = -0.5
ymax =  0.5

dirname = "hlld"

step_s = 50
step_e = 50

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
    print("nx = ",nx)
    print("ny = ",ny)

    x = data[:,0].reshape(ny,nx)
    y = data[:,1].reshape(ny,nx)
    vx = data[:,3].reshape(ny,nx)
    vy = data[:,4].reshape(ny,nx)

    kx = np.fft.fftfreq(nx)*nx/(xmax - xmin)
    ky = np.fft.fftfreq(ny)*ny/(ymax - ymin)
    kx2d, ky2d = np.meshgrid( kx, ky, indexing='ij')
    kwave = np.sqrt(kx2d**2 + ky2d**2)

    kmin = np.min(1.0/(xmax - xmin))
    kmax = np.min(0.5*nx/(xmax - xmin))

    kbins = np.arange(kmin,kmax,kmin)
    N = len(kbins)
    kbinc = 0.5*( kbins[1:] + kbins[:-1] )


    Epower = PowerSpectrum( vx, vy, kwave, kbins )


    plt.xscale('log')
    plt.yscale('log')
    plt.plot(kbinc,Epower,'o-')
    i0 = np.where(kbinc > 10)[0][0]

    plt.plot(kbinc,Epower[i0]*2.0*(kbinc/10)**(-5.0/3.0),'-')
    plt.plot(kbinc,Epower[i0]*2.0*(kbinc/10)**(-3.0),'-')



    output = dirname + "/spect%05d.dat"%(istep) 
    fout = open(output,"w")
    for i in range(N-1): 
        fout.write("{0} {1}\n".format(kbinc[i], Epower[i]))
    fout.close()


plt.show()
