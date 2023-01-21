module globals
    integer, parameter ::  nx=128;
    integer, parameter ::  ngh=2;
    integer, parameter ::  ntot=nx + 2*ngh;
    integer, parameter :: is = ngh + 1 ! index of the leftmost gird
    integer, parameter :: ie = ngh + nx  ! index of the rightmost gird
    integer, parameter :: IDEN=1
    integer, parameter :: IMX=2
    integer, parameter :: IENE=3
    integer, parameter :: NVARS=3

    integer, parameter :: IVX=2
    integer, parameter :: IPRE=3
    real(8), parameter :: CFL = 0.05d0
    real(8), parameter :: gam = 5.0d0/3.0d0
    real(8), parameter :: gamp1 = gam + 1.0d0
    real(8), parameter :: gamm1 = gam - 1.0d0

    real(8), parameter :: xmin = 0.0d0
    real(8), parameter :: xmax = 1.0d0
end module

program Godunov_method
    use globals
    real(8) :: xf(ntot+1), xv(ntot)
    real(8) :: U(ntot,NVARS) ! conserved quantities at x_{i}
    real(8) :: W(ntot,NVARS) ! primitive quantities at x_{i}
    real(8) :: F(ntot+1,NVARS) ! numerical flux at x_{i+1/2}
    real(8) :: Uo(ntot,NVARS) ! conserved quantities at x_{i}
    real(8) :: time, dt
    integer :: itime

    call GenerateGrid( xf, xv )
    call SetInitialCondition( xv, xf, U )
    call Cons2Prim( U, W )
    call ApplyBoundaryCondition( W )

    time = 0.0d0
    itime = 0
    timeout = 2.5

    do while( time < timeout ) 
          Uo(:,:) = U(:,:)
          dt = NewDt( xf, W )

!          call CalculateNumericalFlux( xf, xv, W, F )
!          call UpdateConservedVariables( 0.5d0*dt, xf, Uo, F, U )
!          call Cons2Prim( U, W )
!          call ApplyBoundaryCondition( W )

!          call CalculateNumericalFlux( xf, xv, W, F )
!          call UpdateConservedVariables( 0.5d0*dt, xf, U, F, U )
!          call Cons2Prim( U, W )
!          call ApplyBoundaryCondition( W )

          call CalculateNumericalFlux( xf, xv, W, F )
          call UpdateConservedVariables( dt, xf, Uo, F, U )
          call Cons2Prim( U, W )
          call ApplyBoundaryCondition( W )

          time = time + dt
          itime = itime + 1
!          print*,time,W(2,IDEN)-1.0d0
    enddo

    do i=is,ie
       print*, xv(i), W(i,IDEN), W(i,IVX), W(i,IPRE), F(i,IDEN) - F(i-1,IDEN)
   enddo


contains
!---------------------------------------------------------
!---------------------------------------------------------
Subroutine GenerateGrid( xf, xv )
use globals, only: ntot, is, xmin, xmax
implicit none
real(8), intent(out) :: xf(:), xv(:)
integer :: i
real(8) :: dx


    dx = (xmax - xmin)/dble(nx)
    do i=1, ntot+1
         xf(i) = xmin + dble(i-is)*dx
    end do
    do i=1, ntot
         xv(i) = 0.5d0*(xf(i) + xf(i+1))
    enddo

end Subroutine
!---------------------------------------------------------
!---------------------------------------------------------
Subroutine SetInitialCondition( xv, xf, U )
use globals
implicit none
real(8), intent(in)  :: xv(:), xf(:) 
real(8), intent(out) :: U(:,:)
integer :: i
real(8) :: amp, Lx, kwave,pi
real(8) :: den, vel, pre

    pi = acos(-1.0d0)

    amp = 1d-5
    Lx = xf(ie+1) - xf(is)
    kwave = 2.0d0*pi/Lx

    do i=is,ie
         den = 1.0d0 + amp*sin(kwave*xv(i))
         vel = sqrt(gam)*amp*sin(kwave*xv(i))
         pre = 1.0 + gam*amp*sin(kwave*xv(i))


!      if( xv(i) <0.5d0 ) then
!          den = 1.0d0
!          vel = 0.0d0
!          pre = 1.0d0
!      else 
!          den = 0.125d0
!          vel = 0.0d0
!          pre = 0.1d0
!      endif
         U(i,IDEN) = den
         U(i,IMX ) = den*vel
         U(i,IENE) = pre/(gam - 1.0d0) + 0.5d0*den*vel**2
    enddo

end Subroutine

!---------------------------------------------------------
!---------------------------------------------------------

Subroutine Cons2Prim( U, W )
use globals
implicit none
real(8), intent(in ) :: U(:,:)
real(8), intent(out) :: W(:,:)
real(8) :: deni
integer :: i


    do i=is,ie
        W(i,IDEN) = U(i,IDEN)
        deni = 1.0/U(i,IDEN)

        W(i,IVX) = deni*U(i,IMX)
        W(i,IPRE) = ( U(i,IENE) - 0.5d0*U(i,IMX)**2*deni )*(gam - 1.0d0)
    enddo

end Subroutine

!---------------------------------------------------------
!---------------------------------------------------------

real(8) function NewDt( xf, W )
use globals
implicit none
real(8), intent(in) :: xf(:), W(:,:)
integer :: i

    dt = 1d100
    do i=is, ie
         dt = min( dt, (xf(i+1) - xf(i))/(dsqrt(gam*W(i,IPRE)/W(i,IDEN)) + dabs(W(i,IVX)) ) )
    enddo

    NewDt = CFL*dt

end function

!---------------------------------------------------------
!---------------------------------------------------------

Subroutine ApplyBoundaryCondition( W )
use globals
implicit none
real(8), intent(out) :: W(:,:)
real(8) :: deni
integer :: i

! left boundary 
    do i=1,ngh
!        W(is-i,IDEN) = W(is+i-1,IDEN)
!        W(is-i,IVX ) = W(is+i-1,IVX )
!        W(is-i,IPRE) = W(is+i-1,IPRE)
       W(is-i,IDEN) = W(ie-i+1,IDEN)
        W(is-i,IVX ) = W(ie-i+1,IVX )
        W(is-i,IPRE) = W(ie-i+1,IPRE)
    enddo

! right boundary
    do i=1,ngh
!        W(ie+i,IDEN) = W(ie-i+1,IDEN)
!        W(ie+i,IVX ) = W(ie-i+1,IVX )
!        W(ie+i,IPRE) = W(ie-i+1,IPRE)

        W(ie+i,IDEN) = W(is+i-1,IDEN)
        W(ie+i,IVX ) = W(is+i-1,IVX )
        W(ie+i,IPRE) = W(is+i-1,IPRE)
    enddo

end Subroutine

!---------------------------------------------------------
!---------------------------------------------------------
Subroutine CalculateNumericalFlux( xf, xv, W, F )
use globals
implicit none
real(8), intent(in ) :: xf(:), xv(:), W(:,:)
real(8), intent(out) :: F(:,:)

integer :: i, n
real(8), allocatable :: WL(:,:), WR(:,:)
real(8) :: dwL, dwR, dwC, dW

allocate( WL(ntot,NVARS), WR(ntot,NVARS) )

!---( 1st order )---
do n=1, NVARS
   do i=is-1, ie+1
!       WL(i,n) = W(i,n)
!       WR(i,n) = W(i+1,n)
        dwL = W(i,n) - W(i-1,n)
        dwR = W(i+1,n) - W(i,n)
        dwC = 0.5d0*( dwL + dwR )
        if (dwL*dwR > 0) then
!            dW = sign(dwL)*min( 2.0*dabs(dwL), 2.0*dabs(dwR), dabs(dwC))
            dW = 2.0d0*dwL*dwR/(dwL+dwR+1d-20)
        else
            dW = 0.0d0
        endif
       WL(i,n) = W(i,n) + 0.5d0*dW
       WR(i-1,n) = W(i,n) - 0.5d0*dW
   enddo
enddo

! call GodunovNumericalFlux( is-1, ie, WL, WR, F )
 call HLLNumericalFlux( is-1, ie, WL, WR, F )

end Subroutine

Subroutine GodunovNumericalFlux( iss, iee, WL, WR, F )
use globals
implicit none
integer, intent(in) :: iss, iee
real(8), intent(in ) :: WL(:,:), WR(:,:)
real(8), intent(out) :: F(:,:)
real(8) :: wwL(NVARS), wwR(NVARS)

integer :: i
real(8) :: prest(ntot), velst(ntot)
real(8) :: aL, aR, S(3), rhostL, rhostR

     call ExactRiemanSolver(iss, iee, WL, WR, velst, prest)

     do i=iss, iee
        wwL(:) = WL(i,:)
        wwR(:) = WR(i,:)

        aL = dsqrt(gam*wwL(IPRE)/wwL(IDEN))
        aR = dsqrt(gam*wwR(IPRE)/wwR(IDEN))

        if ( prest(i) >= wwL(IPRE) ) then !shock
            rhostL = wwL(IDEN)*(prest(i)/wwL(IPRE) + gamm1/gamp1)/( gamm1/gamp1*prest(i)/wwL(IPRE) + 1.0d0)
            S(1) = wwL(IVX) - aL*dsqrt( ( gamp1*prest(i)/wwL(IPRE) + gamm1 )/(2.0d0*gam) )
        else !rarefuction 
            rhostL = wwL(IDEN)*(prest(i)/wwL(IPRE))**(1.0d0/gam)
            S(1) = 0.5d0*( wwL(IVX) - aL + velst(i) - dsqrt(gam*prest(i)/rhostL) )
        endif

        S(2) = velst(i) !const dicontinuity

        if(prest(i) >= wwR(IPRE)) then !shock
            rhostR = wwR(IDEN)*(prest(i)/wwR(IPRE) + gamm1/gamp1)/( gamm1/gamp1*prest(i)/wwR(IPRE) + 1.0d0)
            S(3) = wwR(IVX) + aR*dsqrt( ( gamp1*prest(i)/wwR(IPRE) + gamm1 )/(2.0d0*gam) )
        else !rarefuction
            rhostR = wwR(IDEN)*(prest(i)/wwR(IPRE))**(1.0d0/gam)
            S(3) = 0.5d0*( wwR(IVX) + aR + velst(i) + dsqrt(gam*prest(i)/rhostR) )
        endif

        if( S(1) >= 0) then ! W_L
            F(i,IDEN) = wwL(IDEN)*wwL(IVX)
            F(i,IMX ) = wwL(IPRE) + wwL(IDEN)*wwL(IVX)**2
            F(i,IENE) = ( gam*wwL(IPRE)/gamm1 + 0.5d0*wwL(IDEN)*wwL(IVX)**2 )*wwL(IVX)
        else if (S(3) <= 0) then ! W_R
            F(i,IDEN) = wwR(IDEN)*wwR(IVX)
            F(i,IMX ) = wwR(IPRE) + wwR(IDEN)*wwR(IVX)**2
            F(i,IENE) = ( gam*wwR(IPRE)/gamm1 + 0.5d0*wwR(IDEN)*wwR(IVX)**2 )*wwR(IVX)
        else if (S(2) >= 0) then !W_L*
            F(i,IDEN) = rhostL*velst(i)
            F(i,IMX)  = prest(i) + rhostL*velst(i)**2
            F(i,IENE) = ( gam*prest(i)/gamm1 + 0.5d0*rhostL*velst(i)**2 )*velst(i)
        else  !W_R*
            F(i,IDEN) = rhostR*velst(i)
            F(i,IMX)  = prest(i) + rhostR*velst(i)**2
            F(i,IENE) = ( gam*prest(i)/gamm1 + 0.5d0*rhostR*velst(i)**2 )*velst(i)
        endif
     enddo


end Subroutine

Subroutine HLLNumericalFlux( iss, iee, WL, WR, F )
use globals
implicit none
integer, intent(in) :: iss, iee
real(8), intent(in ) :: WL(:,:), WR(:,:)
real(8), intent(out) :: F(:,:)
real(8) :: wwL(NVARS), wwR(NVARS)

integer :: i
real(8) :: prest(ntot), velst(ntot)
real(8) :: cL, cR, eneL, eneR
real(8) :: FrhoL, FmomL, FeneL
real(8) :: FrhoR, FmomR, FeneR, SL, SR

     do i=iss, iee
        wwL(:) = WL(i,:)
        wwR(:) = WR(i,:)

        cL = dsqrt(gam*wwL(IPRE)/wwL(IDEN))
        cR = dsqrt(gam*wwL(IPRE)/wwL(IDEN))

        SL = min(wwL(IVX) - cL, wwR(IVX) - cR)
        SR = max(wwL(IVX) + cL, wwR(IVX) + cR)

        eneL = 0.5*wwL(IDEN)*wwL(IVX)**2 + wwL(IPRE)/(gam - 1.0)
        eneR = 0.5*wwR(IDEN)*wwR(IVX)**2 + wwR(IPRE)/(gam - 1.0)

        FrhoL = wwL(IDEN)*wwL(IVX)
        FmomL = wwL(IPRE) + wwL(IDEN)*wwL(IVX)**2
        FeneL = ( wwL(IPRE) + eneL )*wwL(IVX)

        FrhoR = wwR(IDEN)*wwR(IVX)
        FmomR = wwR(IPRE) + wwR(IDEN)*wwR(IVX)**2
        FeneR = ( gam*wwR(IPRE)/(gam-1.0) + 0.5*wwR(IDEN)*wwR(IVX)**2 )*wwR(IVX)


        if( SL > 0) then ! W_L
            F(i,IDEN) = FrhoL
            F(i,IMX ) = FmomL
            F(i,IENE) = FeneL
        else if (SR <= 0) then ! W_R
            F(i,IDEN) = FrhoR
            F(i,IMX ) = FmomR
            F(i,IENE) = FeneR
        else 
            F(i,IDEN) = (SR*FrhoL - SL*FrhoR + SR*SL*(wwR(IDEN) - wwL(IDEN)))/(SR - SL)
            F(i,IVX) = (SR*FmomL - SL*FmomR + SR*SL*(wwR(IDEN)*wwR(IVX) - wwL(IDEN)*wwL(IVX)) )/(SR - SL)
            F(i,IENE) = (SR*FeneL - SL*FeneR + SR*SL*(0.5*wwR(IDEN)*wwR(IVX)**2 + wwR(IPRE)/(gam - 1.0) &
               - 0.5*wwL(IDEN)*wwL(IVX)**2 - wwL(IPRE)/(gam - 1.0)))/(SR - SL)
        endif
     enddo


end Subroutine


Subroutine UpdateConservedVariables( dt, xf, Uold, F, U )
use globals
real(8), intent(in) :: dt, F(:,:)
real(8), intent(in) :: xf(:)
real(8), intent(inout) :: Uold(:,:)
real(8), intent(out) :: U(:,:)

    do n=1,NVARS
        do i=is,ie
             U(i,n) = Uold(i,n) - dt*( F(i,n) - F(i-1,n) )/(xf(i) - xf(i-1));
        enddo
    enddo

end Subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine ExactRiemanSolver(iss, iee, wL, wR, velst, prest)
use globals
implicit none
integer, intent(in) :: iss, iee
real(8), intent(in) :: wL(:,:), wR(:,:)
real(8), intent(out) :: velst(:), prest(:)
integer :: i
real(8) :: wwL(NVARS), wwR(NVARS)
real(8) :: cL, cR, zR, zL, vR, vL, mL, mR, msqL, msqR
real(8) :: testP, hantei

do i=iss, iee
    wwL(:) = wL(i,:)
    wwR(:) = wR(i,:)

    wwL(IDEN) = 1.0d0
    wwL(IVX) = 0.0d0
    wwL(IPRE) = 1.0d0
    wwR(IDEN) = 0.125d0
    wwR(IVX) = 0.0d0
    wwR(IPRE) = 0.1d0

    cL = dsqrt(gam*wwL(IPRE)*wwL(IDEN)) 
    cR = dsqrt(gam*wwR(IPRE)*wwR(IDEN))

    testP = ( cR*wwL(IPRE) + cL*wwR(IPRE) - cR*cL*(wwL(IVX) - wwR(IVX)) )/(cR + cL)
    if(testP<0) then
        testP = 1d-8
    endif

    zR = 0.0d0
    zL = 0.0d0
    vR = 0.0d0
    vL = 0.0d0
    hantei = 100.0d0
    do while( abs(hantei) > 1d-4 )
        if( testP >= wwL(IPRE) ) then !shock
            msqL = 0.5d0*cL**2*(gamp1*testP/wwL(IPRE) + gamm1)/gam 
            mL = sqrt(msqL) 
            zL = 2.0d0*msqL*mL/(msqL + cL**2)
        else  !rarefaction
            mL = gamm1/(2.0*gam)*(1.0 - testP/wwL(IPRE))/(1.0 - (testP/wwL(IPRE))**(gamm1/(2.0*gam)) )*cL
            zL = cL*(testP/wwL(IPRE))**(1.0 - gamm1/(2.0*gam))
        endif
        if( testP >= wwR(IPRE) ) then !shock
            msqR = 0.5d0*cR**2*(gamp1*testP/wwR(IPRE) + gamm1)/gam 
            mR = sqrt(msqR) 
            zR = 2.0d0*msqR*mR/(msqR + cR**2)
        else ! rarefaction
            mR = gamm1/(2.0*gam)*(1.0 - testP/wwR(IPRE))/(1.0 - (testP/wwR(IPRE))**(gamm1/(2.0*gam)) )*cR
            zR = cR*(testP/wwR(IPRE))**(1.0 - gamm1/(2.0*gam))
        endif


        vR = wwR(IVX) + (testP - wwR(IPRE))/(mR)
        vL = wwL(IVX) - (testP - wwL(IPRE))/(mL)


        hantei = zR*zL*( vR - vL )/( (zR + zL)*testP )
        testP = testP*(1.0 - hantei)
    enddo

    prest(i) = testP
    velst(i) = (zL*vL + zR*vR)/(zL + zR)
enddo


end Subroutine

end program


!import numpy as np
!import os
!
!import matplotlib.pyplot as plt
!import matplotlib.tri as tri
!
!!preparation for the animation
!from matplotlib.animation import ArtistAnimation
!from IPython.display import HTML
!from IPython.display import Image
!
!def InitialProfile(x): 
!    global x0, Lx
!
!    if x>x0*0.8 and x<x0*1.2: 
!        return 1.0 
!    else: 
!        return 0.0
!
!def main():
!    global Lx, ngrid, ntot, c, x0, flag_legend, nghost, ien, ist, gam
!
!    gam = 1.4
!    tout  = 0.25    ! the final time 
!
!    rhoL = 1.0
!    velL = 0.0
!    preL = 1.0
!
!    rhoR = 0.125
!    velR = 0.0
!    preR = 0.1
!
!
!    Lx    = 1.0     ! the box size [0,Lx]
!    ngrid = 64      ! the total numer of the grids
!    x0    = 0.5*Lx ! position of the initial wave front
!    dt_plot = 0.01  ! time interval of plot
!    CFL = 0.8       ! CFL number
!
!    nghost = 2      ! the number of ghost cells
!    ntot = ngrid + 2*nghost ! total cell number including the ghost cells
!    ist = nghost;   ! the first index of the computational domain
!    ien = nghost+ngrid; ! the last index of the computational domain
!
!    gam1 = gam - 1.0
!
!    !definition of variables
!    xv = np.zeros(ntot)   ! cell center   xv[i] --> x_i
!    xf = np.zeros(ntot+1) ! cell boundary xf[i] --> x_{i+1/2}
!    rho = np.zeros(ntot)    ! u^{n+1} u[i] --> u_i^{n+1}
!    vel = np.zeros(ntot)   ! u^{n+1} u[i] --> u_i^{n+1}
!    pre = np.zeros(ntot) ! u^n     uold[i] --> u_i^n
!    mom = np.zeros(ntot)   ! u^{n+1} u[i] --> u_i^{n+1}
!    ene = np.zeros(ntot) ! u^n     uold[i] --> u_i^n
!    rho0 = np.zeros(ntot)    ! u^{n+1} u[i] --> u_i^{n+1}
!    mom0 = np.zeros(ntot)   ! u^{n+1} u[i] --> u_i^{n+1}
!    ene0 = np.zeros(ntot) ! u^n     uold[i] --> u_i^n
!    Frho = np.zeros(ntot) ! u^n     uold[i] --> u_i^n
!    Fmom = np.zeros(ntot) ! u^n     uold[i] --> u_i^n
!    Fene = np.zeros(ntot) ! u^n     uold[i] --> u_i^n
!    u = np.zeros(ntot) ! u^n     uold[i] --> u_i^n
!
!    dx = Lx/ngrid
!    sod_ana = np.loadtxt("sod_t0.25.dat")
!
!    PreparationAnimation()
!    
!    new_dir_path = 'outputs_hydro_FVM/'
!    if not os.path.exists(new_dir_path):
!        os.mkdir(new_dir_path)
!    
!    !cell boundary position
!    for i in range(ntot+1):
!        xf[i] = (i - ist)*dx
!
!    !cell center position
!    for i in range(ntot):
!        xv[i] = 0.5*(xf[i] + xf[i+1])
!
!    !initial profile 
!    for i in range(ntot):
!        rho[i] = (rhoR - rhoL)*0.5*( np.tanh((xv[i]-0.5*Lx)/(0.5*dx)) + 1.0) + rhoL
!        vel[i] = (velR - velL)*0.5*( np.tanh((xv[i]-0.5*Lx)/(0.5*dx)) + 1.0) + velL
!        pre[i] = (preR - preL)*0.5*( np.tanh((xv[i]-0.5*Lx)/(0.5*dx)) + 1.0) + preL
!
!        mom[i] = rho[i]*vel[i]
!        ene[i] = pre[i]/gam1 + 0.5*rho[i]*vel[i]**2
!    
!    time = 0
!    time_plot = 0
!
!    dt = CFL*np.min(dx/(np.sqrt(gam*pre[ist:ien]/rho[ist:ien]) + abs(vel[ist:ien])))
!    
!    flag_legend=True
!    while time < tout: 
!
!        ! substep
!        for i in range(ist,ien): 
!            rho0[i] = rho[i] 
!            mom0[i] = mom[i] 
!            ene0[i] = ene[i] 
!        ! prediction
!        CalculateRoeNumericalFlux(rho, vel, pre, Frho, Fmom, Fene)
!        UpdateConservedVariables(0.5*dt, dx, rho, rho, mom, mom, ene, ene, Frho, Fmom, Fene)
!        ConvertConserved2Primitive( rho, mom, ene, vel, pre )
!        SetBoundaryCondition(xv,rho,vel,pre)
!    
!        ! correction
!        CalculateRoeNumericalFlux(rho, vel, pre, Frho, Fmom, Fene)
!        UpdateConservedVariables(dt, dx, rho, rho0, mom, mom0, ene, ene0, Frho, Fmom, Fene)
!        ConvertConserved2Primitive( rho, mom, ene, vel, pre )
!        SetBoundaryCondition(xv,rho,vel,pre)
!
!        time = time + dt
!        time_plot = time_plot + dt
!    
!        if time_plot > dt_plot:
!            AddPictureToAnimation(time,xv,rho,vel,pre,sod_ana)
!            time_plot = 0.0
!   
!        if  time + dt > tout:
!            dt = tout - time;
!
!        dt = CFL*np.min(dx/(np.sqrt(gam*pre[ist:ien]/rho[ist:ien]) + abs(vel[ist:ien])))
!    
!    for i in range(ist,ien):
!        print xv[i], rho[i], vel[i], pre[i] 
!    animationA = ArtistAnimation(figA, artistA, interval=100)
!    animationA.save('animate_hydro_FVM_2nd.mp4',writer='ffmpeg')
!    HTML(animationA.to_html5_video())
!
!def CalculateRoeNumericalFlux(rho, vel, pre, Frho, Fmom, Fene): 
!    global ien, ist, gam, ntot
!
!    gam1 = gam - 1.0
!
!    rhol = np.zeros(ntot)   ! rhoL(i+1/2)
!    rhor = np.zeros(ntot)   ! rhoR(i+1/2)
!    vell = np.zeros(ntot)   ! velL(i+1/2)
!    velr = np.zeros(ntot)   ! velR(i+1/2)
!    prel = np.zeros(ntot)   ! preL(i+1/2)
!    prer = np.zeros(ntot)   ! preR(i+1/2)
!
!    ! calculate w_{i+1/2,L}
!    for i in range(ist-1,ien): 
!        rhol[i] = rho[i] 
!        vell[i] = vel[i] 
!        prel[i] = pre[i] 
!
!    ! calculate u_{i+1/2,R}
!    for i in range(ist-1,ien): 
!        rhor[i] = rho[i+1] 
!        velr[i] = vel[i+1]
!        prer[i] = pre[i+1] 
!
!    GodunovNumericalFlux(ist,ien,gam, rhol, vell, prel, rhor, velr, prer, Frho, Fmom, Fene)
!
!def UpdateConservedVariables(dt, dx, rho, rho0, mom, mom0, ene, ene0, Frho, Fmom, Fene):
!    global ien, ist, gam, ntot
!
!    nu = dt/dx
!    for i in range(ist,ien): 
!        rho[i] = rho0[i] - nu*( Frho[i] - Frho[i-1] ) 
!        mom[i] = mom0[i] - nu*( Fmom[i] - Fmom[i-1] ) 
!        ene[i] = ene0[i] - nu*( Fene[i] - Fene[i-1] )
!
!def ConvertConserved2Primitive( rho, mom, ene, vel, pre ): 
!    global gam
!
!    gam1 = gam - 1
!    ! conver the conserved variables to the primitive variables 
!    for i in range(ist,ien): 
!        vel[i] = mom[i]/rho[i] 
!        pre[i] = gam1*(ene[i] - 0.5*rho[i]*vel[i]**2)
!
!def SetBoundaryCondition(x,rho,vel,pre):
!    global nghost, ien, ist
!    
!!    for i in range(nghost): 
!!        rho[i] = rho[ien+i-nghost]
!!        vel[i] = vel[ien+i-nghost]
!!        pre[i] = pre[ien+i-nghost]
!!
!!    for i in range(nghost): 
!!        rho[ien+i] = rho[ist+i]
!!        vel[ien+i] = vel[ist+i]
!!        pre[ien+i] = pre[ist+i]
!
!def ExactRiemanSolver(gam, rhol, vell, prel, rhor, velr, prer ):
!    gamp1 = gam + 1.0
!    gamm1 = gam - 1.0
!
!    cL = np.sqrt(gam*prel*rhol) 
!    cR = np.sqrt(gam*prer*rhor)
!
!    testP = ( cR*prel + cL*prer - cR*cL*(velr - vell) )/(cR + cL)
!    if testP<0: testP = 1e-8
!    hantei = 10
!    zR = 0.0
!    zL = 0.0
!    vR = 0.0
!    vL = 0.0
!    while abs(hantei) > 1e-4: 
!        wsqL = 0.5*cL**2*(gamp1*testP/prel + gamm1)/gam 
!        wL = np.sqrt(wsqL) 
!        zL = 2.0*wsqL*wL/(wsqL + cL**2)
!
!        wsqR = 0.5*cR**2*(gamp1*testP/prer + gamm1)/gam 
!        wR = np.sqrt(wsqR) 
!        zR = 2.0*wsqR*wR/(wsqR + cR**2)
!
!        vR = velr + (testP - prer)/(wR)
!        vL = vell - (testP - prel)/(wL)
!
!        hantei = zR*zL*( vR - vL )/( (zR + zL)*testP )
!        testP = testP*(1.0 - hantei)
!
!    return testP, (zL*vL + zR*vR)/(zL + zR)
!
!def GodunovNumericalFlux(ist,ien,gam, rhol, vell, prel, rhor, velr, prer, Frho, Fmom, Fene):
!
!    for i in range(ist-1,ien): 
!        gamp1 = gam + 1.0
!        gamm1 = gam - 1.0
!    
!        prest, velst = ExactRiemanSolver(gam, rhol[i], vell[i], prel[i], rhor[i], velr[i], prer[i] )
!    
!        aL = np.sqrt(gam*prel[i]/rhol[i])
!        aR = np.sqrt(gam*prer[i]/rhor[i])
!    
!        S = np.zeros(3)
!    
!        if prest >= prel[i]: !shock
!            rhostL = rhol[i]*(prest/prel[i] + gamm1/gamp1)/( gamm1/gamp1*prest/prel[i] + 1.0)
!            S[0] = vell[i] - aL*np.sqrt( ( gamp1*prest/prel[i] + gamm1 )/(2.0*gam) )
!        else: !rarefuction 
!            rhostL = rhol[i]*(prest/prel[i])**(1.0/gam)
!            S[0] = 0.5*( vell[i] - aL + velst - np.sqrt(gam*prest/rhostL) )
!        
!        S[1] = velst !const dicontinuity
!    
!        if prest >= prer[i]: !shock
!            rhostR = rhor[i]*(prest/prer[i] + gamm1/gamp1)/( gamm1/gamp1*prest/prer[i] + 1.0)
!            S[2] = velr[i] + aR*np.sqrt( ( gamp1*prest/prer[i] + gamm1 )/(2.0*gam) )
!        else: !rarefuction
!            rhostR = rhor[i]*(prest/prer[i])**(1.0/gam)
!            S[2] = 0.5*( velr[i] + aR + velst + np.sqrt(gam*prest/rhostR) )
!
!        if S[0] >= 0: ! W_L
!            Frho[i] = rhol[i]*vell[i]
!            Fmom[i] = prel[i] + rhol[i]*vell[i]**2
!            Fene[i] = ( gam*prel[i]/(gam-1.0) + 0.5*rhol[i]*vell[i]**2 )*vell[i]
!        elif S[2] <= 0: ! W_R
!            Frho[i] = rhor[i]*velr[i]
!            Fmom[i] = prer[i] + rhor[i]*velr[i]**2
!            Fene[i] = ( gam*prer[i]/(gam-1.0) + 0.5*rhor[i]*velr[i]**2 )*velr[i]
!        elif S[1] >= 0: ! W_L*
!            Frho[i] = rhostL*velst
!            Fmom[i] = prest + rhostL*velst**2
!            Fene[i] = ( gam*prest/(gam-1.0) + 0.5*rhostL*velst**2 )*velst
!        else:
!            Frho[i] = rhostR*velst
!            Fmom[i] = prest + rhostR*velst**2
!            Fene[i] = ( gam*prest/(gam-1.0) + 0.5*rhostR*velst**2 )*velst
!
!            
!
!         
!        
!    
!
!
!def PreparationAnimation():
!    global figA, artistA, ax1, Lx, c, x0
!    !set up for animation
!    figA = plt.figure(figsize=(8,6))
!    plt.rcParams['font.size']=20
!    artistA = []
!    ax1 = figA.add_subplot(1,1,1)
!
!    !set up for animation
!    figA = plt.figure(figsize=(8,6))
!    plt.rcParams['font.size']=20
!    artistA = []
!    ax1 = figA.add_subplot(1,1,1)
!
!def AddPictureToAnimation(time,xv,rho,vel,pre,sod_ana): 
!    global ax1, ntot, artistA, Lx, x0, flag_legend, gam
!    timetxt='time=%.2f'%(time)
!    ax1.set_xlim(0,Lx)
!    ax1.set_ylim(0,1.1)
!    ax1.set_xlabel(r'$x$')
!    title = ax1.text(0.5*Lx,2.1,timetxt,horizontalalignment="center")
!    im1 = ax1.plot(xv,rho,'o-',color='red',label=r"$\rho$")
!    im2 = ax1.plot(xv,vel,'o-',color='green',label=r"$v$")
!    im3 = ax1.plot(xv,pre,'o-',color='blue',label=r"$P$")
!
!    im4 = ax1.plot((sod_ana[:,0]-0.5)*time/0.25+0.5,sod_ana[:,1],'-',color='black')
!    im5 = ax1.plot((sod_ana[:,0]-0.5)*time/0.25+0.5,sod_ana[:,2],'-',color='black')
!    im6 = ax1.plot((sod_ana[:,0]-0.5)*time/0.25+0.5,sod_ana[:,3],'-',color='black')
!
!!    im2 = ax1.plot(xv,usol,'-',color='blue',label="analytic solution")
!    if flag_legend: 
!        ax1.legend(loc='upper right') 
!        flag_legend = False
!    
!    figA.tight_layout()
!    artistA.append(im1+im2+im3+im4+im5+im6+[title])
!
!if __name__ == "__main__":
!    main()
