#R1 (Wf+) R2 (Wa+) R3 (Ws+) R4 (W0) R5 (Ws-) R6 (Wa-) R7 (Wf-) R8
import numpy as np

gam = 5.0/3.0
Bx = 0.75

x0 = 0.5
time = 0.10

IDN = 0
IPR = 1
IVX = 2
IVY = 3
IVZ = 4
IBY = 5
IBZ = 6

Q1 = np.array([1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0])
Q2 = np.array([6.5673e-1, 4.9618e-1, 6.5790e-1, -2.6700e-1, 0.0, 5.4966e-1, 0.0])
Q3 = np.array([6.5673e-1, 4.9618e-1, 6.5790e-1, -1.6235e0, 0.0, -5.4966e-1, 0.0])
Q4 = np.array([6.6535e-1, 5.0709e-1, 6.4841e-1, -1.6053e0, 0.0, -5.3799e-1, 0.0])
Q5 = np.array([2.74e-1  , 5.0709e-1, 6.4841e-1, -1.6053e0, 0.0, -5.3799e-1, 0.0])
Q6 = np.array([1.1571e-1, 8.7927e-2, -2.7717e-1, -1.9847e-1, 0.0, -8.8576e-1, 0.0])
Q7 = np.array([1.1571e-1, 8.7927e-2, -2.7717e-1, -1.9847e-1, 0.0, -8.8576e-1, 0.0])
Q8 = np.array([0.125, 0.1, 0, 0.0, 0.0, -1.0, 0.0])

ca1 = Bx/np.sqrt(Q1[IDN])
ca1_tot = np.sqrt( Bx**2 + Q1[IBY]**2 + Q1[IBZ]**2)/np.sqrt(Q1[IDN])
a1 = np.sqrt(gam*Q1[IPR]/Q1[IDN])
cf1 = np.sqrt( 0.5*( ca1_tot**2 + a1**2 ) + np.sqrt( 0.25*(ca1_tot**2 + a1**2)**2 - a1**2*ca1**2 ) )
cs1 = np.sqrt( 0.5*( ca1_tot**2 + a1**2 ) - np.sqrt( 0.25*(ca1_tot**2 + a1**2)**2 - a1**2*ca1**2 ) )

print( 0.0, Q1[0], Q1[1], Q1[2], Q1[3], Q1[4], Bx, Q1[5], Q1[6] )
print( 0.5 + (Q1[IVX] - cf1)*time, Q1[0], Q1[1], Q1[2], Q1[3], Q1[4], Bx, Q1[5], Q1[6] )

iend = 1000
Bperp = np.linspace(Q1[IBY],Q2[IBY],iend)
dBperp = Bperp[1] - Bperp[0]

den = Q1[IDN]
vx  = Q1[IVX]
vy  = Q1[IVY]
pre = Q1[IPR]
eps = -1.0
for i in range(iend):
    deno = den
    vxo  = vx
    vyo  = vy
    preo  = pre
    a0 = np.sqrt(gam*pre/den)
    ca0 = Bx/np.sqrt(den)
    ca0_tot = np.sqrt(Bx**2 + Bperp[i]**2)/np.sqrt(den)
    cf0 = np.sqrt( 0.5*( ca0_tot**2 + a0**2 ) + np.sqrt( 0.25*(ca0_tot**2 + a0**2)**2 - a0**2*ca0**2 ) )

    den = deno + den*(cf0**2 - ca0**2)/(cf0**2*Bperp[i])*dBperp
    vx = vxo + (cf0**2 - ca0**2)/(eps*cf0*Bperp[i])*dBperp
    vy = vyo - Bx/(eps*cf0*deno)*dBperp
    pre = preo + eps*cf0*deno*(vx - vxo) - Bperp[i]*dBperp

    print( 0.5 + (vx - cf0)*time, den, pre, vx, vy, 0.0, Bx, Bperp[i], 0.0)

#print("den",den,Q2[IDN])
#print("vx",vx/Q2[IVX])
#print("vy",vy/Q2[IVY])
#print("pre",pre/Q2[IPR])

ca2 = Bx/np.sqrt(Q2[IDN])
ca2_tot = np.sqrt( Bx**2 + Q2[IBY]**2 + Q2[IBZ]**2)/np.sqrt(Q2[IDN])
a2 = np.sqrt(gam*Q2[IPR]/Q2[IDN])
cf2 = np.sqrt( 0.5*( ca2_tot**2 + a2**2 ) + np.sqrt( 0.25*(ca2_tot**2 + a2**2)**2 - a2**2*ca2**2 ) )
cs2 = np.sqrt( 0.5*( ca2_tot**2 + a2**2 ) - np.sqrt( 0.25*(ca2_tot**2 + a2**2)**2 - a2**2*ca2**2 ) )

print( 0.5 + (Q2[IVX] - cf2)*time, Q2[0], Q2[1], Q2[2], Q2[3], Q2[4], Bx, Q2[5], Q2[6] )

#Q1[5:6] = Q1[5:6]/np.sqrt(4.0*np.pi)
#Q2[5:6] = Q2[5:6]/np.sqrt(4.0*np.pi)
#Q3[5:6] = Q3[5:6]/np.sqrt(4.0*np.pi)
#Q4[5:6] = Q4[5:6]/np.sqrt(4.0*np.pi)
#Q5[5:6] = Q5[5:6]/np.sqrt(4.0*np.pi)
#Q6[5:6] = Q6[5:6]/np.sqrt(4.0*np.pi)
#Q7[5:6] = Q7[5:6]/np.sqrt(4.0*np.pi)
#Q8[5:6] = Q8[5:6]/np.sqrt(4.0*np.pi)

Sfp = (Q2[IDN]*Q2[IVX] - Q1[IDN]*Q1[IVX])/(Q2[IDN] - Q1[IDN])
Sap = Q2[IVX] - Bx/np.sqrt(Q2[IDN])
Ssp = (Q4[IDN]*Q4[IVX] - Q3[IDN]*Q3[IVX])/(Q4[IDN] - Q3[IDN])

Sc = Q4[IVX]

Ssm = (Q5[IDN]*Q5[IVX] - Q6[IDN]*Q6[IVX])/(Q5[IDN] - Q6[IDN])
Sam = Q6[IVX] + Bx/np.sqrt(Q6[IDN])
Sfm = (Q7[IDN]*Q7[IVX] - Q8[IDN]*Q8[IVX])/(Q7[IDN] - Q8[IDN])

xfp = 0.5+Sfp*time
xap = 0.5+Sap*time
xsp = 0.5+Ssp*time
xc = 0.5+Sc*time
xsm = 0.5+Ssm*time
xam = 0.5+Sam*time
xfm = 0.5+Sfm*time

#print( xfp, Q1[0], Q1[1], Q1[2], Q1[3], Q1[4], Bx, Q1[5], Q1[6] )
#print( xfp, Q2[0], Q2[1], Q2[2], Q2[3], Q2[4], Bx, Q2[5], Q2[6] )
print( xap, Q2[0], Q2[1], Q2[2], Q2[3], Q2[4], Bx, Q2[5], Q2[6] )
print( xap, Q3[0], Q3[1], Q3[2], Q3[3], Q3[4], Bx, Q3[5], Q3[6] )
print( xsp, Q3[0], Q3[1], Q3[2], Q3[3], Q3[4], Bx, Q3[5], Q3[6] )
print( xsp, Q4[0], Q4[1], Q4[2], Q4[3], Q4[4], Bx, Q4[5], Q4[6] )
print( xc , Q4[0], Q4[1], Q4[2], Q4[3], Q4[4], Bx, Q4[5], Q4[6] )
print( xc , Q5[0], Q5[1], Q5[2], Q5[3], Q5[4], Bx, Q5[5], Q5[6] )
print( xsm , Q5[0], Q5[1], Q5[2], Q5[3], Q5[4], Bx, Q5[5], Q5[6] )
print( xsm , Q6[0], Q6[1], Q6[2], Q6[3], Q6[4], Bx, Q6[5], Q6[6] )
print( xam , Q6[0], Q6[1], Q6[2], Q6[3], Q6[4], Bx, Q6[5], Q6[6] )
print( xam , Q7[0], Q7[1], Q7[2], Q7[3], Q7[4], Bx, Q7[5], Q7[6] )
#print( xfm , Q7[0], Q7[1], Q7[2], Q7[3], Q7[4], Bx, Q7[5], Q7[6] )
#print( xfm , Q8[0], Q8[1], Q8[2], Q8[3], Q8[4], Bx, Q8[5], Q8[6] )

ca2 = Bx/np.sqrt(Q7[IDN])
ca2_tot = np.sqrt( Bx**2 + Q7[IBY]**2 + Q7[IBZ]**2)/np.sqrt(Q7[IDN])
a2 = np.sqrt(gam*Q7[IPR]/Q7[IDN])
cf2 = np.sqrt( 0.5*( ca2_tot**2 + a2**2 ) + np.sqrt( 0.25*(ca2_tot**2 + a2**2)**2 - a2**2*ca2**2 ) )
cs2 = np.sqrt( 0.5*( ca2_tot**2 + a2**2 ) - np.sqrt( 0.25*(ca2_tot**2 + a2**2)**2 - a2**2*ca2**2 ) )

print( 0.5 + (Q7[IVX] + cf2)*time, Q7[0], Q7[1], Q7[2], Q7[3], Q7[4], Bx, Q7[5], Q7[6] )

iend = 100
Bperp = np.linspace(Q8[IBY],Q7[IBY],iend)
dBperp = Bperp[1] - Bperp[0]

den = Q8[IDN]
vx  = Q8[IVX]
vy  = Q8[IVY]
pre = Q8[IPR]
eps = 1.0
for i in range(iend):
    deno = den
    vxo  = vx
    vyo  = vy
    preo  = pre
    a0 = np.sqrt(gam*pre/den)
    ca0 = Bx/np.sqrt(den)
    ca0_tot = np.sqrt(Bx**2 + Bperp[i]**2)/np.sqrt(den)
    cf0 = np.sqrt( 0.5*( ca0_tot**2 + a0**2 ) + np.sqrt( 0.25*(ca0_tot**2 + a0**2)**2 - a0**2*ca0**2 ) )

    den = deno + den*(cf0**2 - ca0**2)/(cf0**2*Bperp[i])*dBperp
    vx = vxo + (cf0**2 - ca0**2)/(eps*cf0*Bperp[i])*dBperp
    vy = vyo - Bx/(eps*cf0*deno)*dBperp
    pre = preo + eps*cf0*deno*(vx - vxo) - Bperp[i]*dBperp

    print( 0.5 + (vx + cf0)*time, den, pre, vx, vy, 0.0, Bx, Bperp[i], 0.0)

ca2 = Bx/np.sqrt(Q8[IDN])
ca2_tot = np.sqrt( Bx**2 + Q8[IBY]**2 + Q8[IBZ]**2)/np.sqrt(Q8[IDN])
a2 = np.sqrt(gam*Q8[IPR]/Q8[IDN])
cf2 = np.sqrt( 0.5*( ca2_tot**2 + a2**2 ) + np.sqrt( 0.25*(ca2_tot**2 + a2**2)**2 - a2**2*ca2**2 ) )
cs2 = np.sqrt( 0.5*( ca2_tot**2 + a2**2 ) - np.sqrt( 0.25*(ca2_tot**2 + a2**2)**2 - a2**2*ca2**2 ) )

print( 0.5 + (Q8[IVX] + cf2)*time, Q8[0], Q8[1], Q8[2], Q8[3], Q8[4], Bx, Q8[5], Q8[6] )
print( 1.0 , Q8[0], Q8[1], Q8[2], Q8[3], Q8[4], Bx, Q8[5], Q8[6] )

