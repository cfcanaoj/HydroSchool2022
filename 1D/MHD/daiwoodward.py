#R1 (Wf+) R2 (Wa+) R3 (Ws+) R4 (W0) R5 (Ws-) R6 (Wa-) R7 (Wf-) R8
import numpy as np

gam = 5.0/3.0
Bx = 2.0/np.sqrt(4.0*np.pi)

IDN = 0
IPR = 1
IVX = 2
IVY = 3
IVZ = 4
IBY = 5
IBZ = 6

Q1 = np.array([1.08, 0.95, 1.2, 0.01, 0.5, 3.6, 2.0])
Q2 = np.array([1.4903, 1.6558e0, 6.0588e-1, 1.1235e-1, 5.5686e-1, 5.0987e0, 2.8326e0])
Q3 = np.array([1.4903, 1.6558e0, 6.0588e-1, 2.2157e-1, 3.0125e-1, 5.5713e0, 1.7264e0])
Q4 = np.array([1.6343, 1.9317e0, 5.7538e-1, 4.7601e-2, 2.4734e-1, 5.0074e0, 1.5517e0])
Q5 = np.array([1.4735, 1.9317e0, 5.7538e-1, 4.7601e-2, 2.4734e-1, 5.0074e0, 1.5517e0])
Q6 = np.array([1.3090, 1.5844e0, 5.3432e-1, -1.8411e-1, 1.7554e-1, 5.7083e0, 1.7689e0])
Q7 = np.array([1.3090, 1.5844e0, 5.3432e-1, -9.4572e-2, -4.7286e-2, 5.3452e0, 2.6726e0])
Q8 = np.array([1.0, 1.0, 0, 0.0, 0.0, 4.0, 2.0])

Q1[5:7] = Q1[5:7]/np.sqrt(4.0*np.pi)
Q2[5:7] = Q2[5:7]/np.sqrt(4.0*np.pi)
Q3[5:7] = Q3[5:7]/np.sqrt(4.0*np.pi)
Q4[5:7] = Q4[5:7]/np.sqrt(4.0*np.pi)
Q5[5:7] = Q5[5:7]/np.sqrt(4.0*np.pi)
Q6[5:7] = Q6[5:7]/np.sqrt(4.0*np.pi)
Q7[5:7] = Q7[5:7]/np.sqrt(4.0*np.pi)
Q8[5:7] = Q8[5:7]/np.sqrt(4.0*np.pi)

Sfp = (Q2[IDN]*Q2[IVX] - Q1[IDN]*Q1[IVX])/(Q2[IDN] - Q1[IDN])
Sap = Q2[IVX] - Bx/np.sqrt(Q2[IDN])
Ssp = (Q4[IDN]*Q4[IVX] - Q3[IDN]*Q3[IVX])/(Q4[IDN] - Q3[IDN])

Sc = Q4[IVX]

Ssm = (Q5[IDN]*Q5[IVX] - Q6[IDN]*Q6[IVX])/(Q5[IDN] - Q6[IDN])
Sam = Q6[IVX] + Bx/np.sqrt(Q6[IDN])
Sfm = (Q7[IDN]*Q7[IVX] - Q8[IDN]*Q8[IVX])/(Q7[IDN] - Q8[IDN])

x0 = 0.5
time = 0.2
xfp = 0.5+Sfp*time
xap = 0.5+Sap*time
xsp = 0.5+Ssp*time
xc = 0.5+Sc*time
xsm = 0.5+Ssm*time
xam = 0.5+Sam*time
xfm = 0.5+Sfm*time

print( 0.0, Q1[0], Q1[1], Q1[2], Q1[3], Q1[4], Bx, Q1[5], Q1[6] )
print( xfp, Q1[0], Q1[1], Q1[2], Q1[3], Q1[4], Bx, Q1[5], Q1[6] )
print( xfp, Q2[0], Q2[1], Q2[2], Q2[3], Q2[4], Bx, Q2[5], Q2[6] )
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
print( xfm , Q7[0], Q7[1], Q7[2], Q7[3], Q7[4], Bx, Q7[5], Q7[6] )
print( xfm , Q8[0], Q8[1], Q8[2], Q8[3], Q8[4], Bx, Q8[5], Q8[6] )
print( 1.0 , Q8[0], Q8[1], Q8[2], Q8[3], Q8[4], Bx, Q8[5], Q8[6] )


