import numpy as np
import os

import matplotlib.pyplot as plt
import matplotlib.tri as tri

#preparation for the animation
from matplotlib.animation import ArtistAnimation
#from IPython.display import HTML
#from IPython.display import Image

def main():
    gam = 1.4
#    gam = 5.0/3.0
#    time = 0.15
    time = 0.2
    x0 = 0.0

    rhoL = 1.0
    velL = 0.0
    preL = 1.0

    rhoR = 0.125
    velR = 0.0
    preR = 0.1

#    rhoL = 1.0
#    velL = 0.0
#    preL = 1.0
#
#    rhoR = 0.125
#    velR = 0.0
#    preR = 0.1

    x, rho, vel, pre = ExactRiemanSolution(time, x0, gam, rhoL, velL, preL, rhoR, velR, preR)

    figA = plt.figure(figsize=(8,6))
    plt.rcParams['font.size']=20
    artistA = []
    ax1 = figA.add_subplot(1,1,1)

    ax1.set_xlim(-0.5,0.5)
    ax1.set_xlabel(r'$x$')
    im1 = ax1.plot(x,rho,'-',color='red',linewidth=4,label=r"$\rho$")
    im2 = ax1.plot(x,vel,'-',color='green',linewidth=4,label=r"$v$")
    im3 = ax1.plot(x,pre,'-',color='blue',linewidth=4,label=r"$P$")

    for i in range(len(x)):
        print x[i], rho[i], vel[i], pre[i]

#    figA.savefig('riemann_sol.png')
#    Image(filename='riemann_sol.png')

def ExactRiemanSolver(gam, rhoL, velL, preL, rhoR, velR, preR ):
    gamp1 = gam + 1.0
    gamm1 = gam - 1.0

    cL = np.sqrt(gam*preL*rhoL) 
    cR = np.sqrt(gam*preR*rhoR)

#    testP = 0.5*(preL + preR)
    testP = ( cR*preL + cL*preR - cR*cL*(velR - velL) )/(cR + cL)
    if testP<0: testP = 1e-8
    hantei = 10
    zR = 0.0
    zL = 0.0
    vR = 0.0
    vL = 0.0
    while abs(hantei) > 1e-4:
        if testP >= preL:   #shock
            wsqL = 0.5*cL**2*(gamp1*testP/preL + gamm1)/gam 
            wL = np.sqrt(wsqL) 
            zL = 2.0*wsqL*wL/(wsqL + cL**2)
        else: #rarefaction
            wL = gamm1/(2.0*gam)*(1.0 - testP/preL)/(1.0 - (testP/preL)**(gamm1/(2.0*gam)) )*cL
            zL = cL*(testP/preL)**(1.0 - gamm1/(2.0*gam))

        if testP >= preR:  #shock
            wsqR = 0.5*cR**2*(gamp1*testP/preR + gamm1)/gam 
            wR = np.sqrt(wsqR) 
            zR = 2.0*wsqR*wR/(wsqR + cR**2)
        else: #rarefunction
            wR = gamm1/(2.0*gam)*(1.0 - testP/preR)/(1.0 - (testP/preR)**(gamm1/(2.0*gam)) )*cR
            zR = cR*(testP/preR)**(1.0 - gamm1/(2.0*gam))


        vR = velR + (testP - preR)/(wR)
        vL = velL - (testP - preL)/(wL)

        hantei = zR*zL*( vR - vL )/( (zR + zL)*testP )
        testP = testP*(1.0 - hantei)



    return testP, (zL*vL + zR*vR)/(zL + zR)

def ExactRiemanSolution(time, x0, gam, rhoL, velL, preL, rhoR, velR, preR ):
    gamp1 = gam + 1.0
    gamm1 = gam - 1.0

    prest, velst = ExactRiemanSolver(gam, rhoL, velL, preL, rhoR, velR, preR )

    aL = np.sqrt(gam*preL/rhoL)
    aR = np.sqrt(gam*preR/rhoR)

    is_shock_L = prest > preL
    is_shock_R = prest > preR

    if is_shock_L: #shock
        rhostL = rhoL*(prest/preL + gamm1/gamp1)/( gamm1/gamp1*prest/preL + 1.0)
        sL = velL - aL*np.sqrt( ( gamp1*prest/preL + gamm1 )/(2.0*gam) )

    else: #rarefuction
        rhostL = rhoL*(prest/preL)**(1.0/gam)
        sL_head = velL - aL
        sL_tail = velst - np.sqrt(gam*prest/rhostL)

    if is_shock_R: #shock
        rhostR = rhoR*(prest/preR + gamm1/gamp1)/( gamm1/gamp1*prest/preR + 1.0)
        sR = velR + aR*np.sqrt( ( gamp1*prest/preR + gamm1 )/(2.0*gam) )

    else: #rarefuction
        rhostR = rhoR*(prest/preR)**(1.0/gam)
        sR_head = velR + aR
        sR_tail = velst + np.sqrt(gam*prest/rhostR)

#    print "rhostR =",  rhostR
#    print "rhostL =",  rhostL

    x_riem = []
    rho_riem = []
    vel_riem = []
    pre_riem = []

    n_rarefac = 64

    if is_shock_L:
        x_riem += [sL*time-0.5, sL*time, sL*time]
        rho_riem += [rhoL, rhoL, rhostL]
        vel_riem += [velL, velL, velst]
        pre_riem += [preL, preL, prest]
    else:
        x_riem += [sL_head*time - 0.5, sL_head*time]
        rho_riem += [rhoL,rhoL]
        vel_riem += [velL,velL]
        pre_riem += [preL,preL]

        dx = -(sL_head - sL_tail)*time/(n_rarefac-1.0)
        for i in range(n_rarefac):
            x0 = sL_head*time + i*dx
            x_riem.append(x0)
            rho_riem.append( rhoL*(2.0/gamp1 + gamm1/(gamp1*aL)*(velL - x0/time) )**(2.0/gamm1) )
            vel_riem.append( 2.0/gamp1*(aL + gamm1/2.0*velL + x0/time ) )
            pre_riem.append( preL*(2.0/gamp1 + gamm1/(gamp1*aL)*(velL - x0/time) )**(2.0*gam/gamm1) )

    #contact discontinuity
    x_riem += [velst*time, velst*time] 
    rho_riem += [rhostL, rhostR]
    vel_riem += [velst, velst]
    pre_riem += [prest, prest]

    if is_shock_R:
        x_riem += [sR*time, sR*time, sR*time + 0.5]
        rho_riem += [rhostR, rhoR, rhoR]
        vel_riem += [velst, velR, velR]
        pre_riem += [prest, preR, preR]
    else:
        dx = (sR_head - sR_tail)*time/(n_rarefac-1.0)
        for i in range(n_rarefac):
            x0 = sR_tail*time + i*dx
            x_riem.append(x0)
            rho_riem.append( rhoR*(2.0/gamp1 - gamm1/(gamp1*aR)*(velR - x0/time) )**(2.0/gamm1) )
            vel_riem.append( 2.0/gamp1*(-aR + gamm1/2.0*velR + x0/time ) )
            pre_riem.append( preR*(2.0/gamp1 - gamm1/(gamp1*aR)*(velR - x0/time) )**(2.0*gam/gamm1) )
        x_riem.append(sR_head*time+0.5)
        rho_riem.append(rhoR)
        vel_riem.append(velR)
        pre_riem.append(preR)


    return x_riem, rho_riem, vel_riem, pre_riem


def PreparationAnimation():
    global figA, artistA, ax1
    #set up for animation
    figA = plt.figure(figsize=(8,6))
    plt.rcParams['font.size']=20
    artistA = []
    ax1 = figA.add_subplot(1,1,1)

    #set up for animation
    figA = plt.figure(figsize=(8,6))
    plt.rcParams['font.size']=20
    artistA = []
    ax1 = figA.add_subplot(1,1,1)

def AddPictureToAnimation(time,x,rho,vel,pre): 
    global ax1, artistA
    ax1.set_xlim(-0.5,0.5)
    ax1.set_xlabel(r'$x$')
    im1 = ax1.plot(x,rho,'-',color='red',linewidth=4,label=r"$\rho$")
    im2 = ax1.plot(x,vel,'-',color='green',linewidth=4,label=r"$v$")
    im3 = ax1.plot(x,pre,'-',color='blue',linewidth=4,label=r"$P$")

    ax1.legend(loc='upper right') 
    
    figA.tight_layout()
    artistA.append(im1+im2+im3)

if __name__ == "__main__":
    main()
