module globals
    integer nx, nxtot;
    integer ngh;
    integer ny, nytot;
    integer is, ie, js, je
    integer, parameter :: IDEN=1
    integer, parameter :: IMX=2
    integer, parameter :: IMY=3
    integer, parameter :: IENE=4
    integer, parameter :: NVARS=4

    integer, parameter :: IVX=2
    integer, parameter :: IVY=3
    integer, parameter :: IPRE=4
    real(8), parameter :: CFL = 0.4d0/2.0d0/2.0d0/2.0d0/2.0d0/2.0d0/2.0d0
    real(8), parameter :: gam = 5.0d0/3.0d0
    real(8), parameter :: gamp1 = gam + 1.0d0
    real(8), parameter :: gamm1 = gam - 1.0d0
end module

program Godunov_method
    use globals
    implicit none
    real(8), allocatable :: xf(:), xv(:)
    real(8), allocatable :: yf(:), yv(:)
    real(8), allocatable :: U(:,:,:) ! conserved quantities at x_{i}
    real(8), allocatable :: W(:,:,:) ! primitive quantities at x_{i}
    real(8), allocatable :: Fx(:,:,:) ! numerical flux at x_{i+1/2}
    real(8), allocatable :: Fy(:,:,:) ! numerical flux at x_{i+1/2}
    real(8), allocatable :: Uo(:,:,:) ! conserved quantities at x_{i}
    integer :: i, j
    real(8) :: xmin, xmax, pi, dx, dy
    real(8) :: ymin, ymax
    real(8) :: time, dt, timeout
    integer :: itime
    real(8) :: amp, kx, ky, tana, sina, cosa,denana, error, omega, Lx

    pi = acos(-1.0d0)

    nx = 64*2*2*2/2/2/2/2/2/2
    ny = nx/2
    ngh = 2
    nxtot = nx + 2*ngh
    is = ngh + 1 ! index of the leftmost gird
    ie = ngh + nx  ! index of the rightmost gird

    nytot = ny + 2*ngh
    js = ngh + 1 ! index of the leftmost gird
    je = ngh + ny  ! index of the rightmost gird

    allocate( xf(nxtot+1) ) 
    allocate( xv(nxtot) ) 
    allocate( yf(nytot+1) ) 
    allocate( yv(nytot) ) 
    allocate( U(nxtot,nytot,NVARS) ) 
    allocate( Uo(nxtot,nytot,NVARS) ) 
    allocate( W(nxtot,nytot,NVARS) ) 
    allocate( Fx(nxtot+1,nytot,NVARS) ) 
    allocate( Fy(nxtot,nytot+1,NVARS) ) 

    xmin = 0.0d0
    xmax = 1.0d0
    dx = (xmax - xmin)/dble(nx)
    do i=1, nxtot+1
         xf(i) = xmin + dble(i-is)*dx
    end do
    do i=1, nxtot
         xv(i) = 0.5d0*(xf(i) + xf(i+1))
    enddo

    ymin = 0.0d0
    ymax = 0.5d0
    dy = (ymax - ymin)/dble(ny)
    do j=1, nytot+1
         yf(j) = ymin + dble(j-js)*dy
    end do
    do j=1, nytot
         yv(j) = 0.5d0*(yf(j) + yf(j+1))
    enddo

    call SetInitialCondition( xv, yv, xf, yf, U )
    call Cons2Prim( U, W )
    call ApplyBoundaryCondition( W )

    time = 0.0d0
    itime = 0
    timeout = 1.0/dsqrt(5.0d0)
!
    do while( time < timeout ) 
          Uo(:,:,:) = U(:,:,:)
          dt = NewDt( xf, yf, W )
          dt =  min( dt, timeout - time)

          call CalculateNumericalFlux( xf, yf, xv, yv, W, Fx, Fy )
          call UpdateConservedVariables( 0.5*dt, xf, yf, U, Fx, Fy, U )
          call Cons2Prim( U, W )
          call ApplyBoundaryCondition( W )

          call CalculateNumericalFlux( xf, yf, xv, yv, W, Fx, Fy )
          call UpdateConservedVariables( dt, xf, yf, Uo, Fx, Fy, U )
          call Cons2Prim( U, W )
          call ApplyBoundaryCondition( W )

          time = time + dt
          itime = itime + 1
!          print*,time
          print*, time, W(is,js,IDEN)
!if(itime.gt.10)exit
    enddo

    error = 0.0d0

    tana = 2.0d0
    cosa = 1.0d0/dsqrt(1.0d0 + tana**2)
    sina = tana/dsqrt(1.0d0 + tana**2)

    amp = 1d-5
    kx = 2.0d0*pi
    ky = kx*tana
    omega = dsqrt(kx**2+ky**2)

    error = 0.0d0
    do j=js,je
    do i=is,ie
        denana = 1.0d0 + 1d-5*dsin(kx*xv(i) + ky*yv(j) - omega*time)
        error = error + (W(i,j,IDEN) - denana)**2
!       print*, xv(i), yv(j), W(i,j,IDEN), W(i,j,IVX), W(i,j,IVY), W(i,j,IPRE)
    enddo
    enddo

    print*, time/timeout, nx, dsqrt(error/(nx*ny))


contains
!---------------------------------------------------------
!---------------------------------------------------------
Subroutine SetInitialCondition( xv, yv, xf, yf, U )
use globals
implicit none
real(8), intent(in)  :: xv(:), xf(:) 
real(8), intent(in)  :: yv(:), yf(:) 
real(8), intent(out) :: U(:,:,:)
integer :: i
real(8) :: amp, Lx, kx,ky,pi
real(8) :: den, vel, pre, k_dot_x
real(8) :: tana,cosa, sina


    pi = acos(-1.0d0)

    tana = 2.0d0
    cosa = 1.0/dsqrt(1.0 + tana**2)
    sina = tana/dsqrt(1.0 + tana**2)

    amp = 1d-5
    Lx = xf(ie+1) - xf(is)
    kx = 2.0d0*pi
    ky = kx*tana

    do j=js,je
    do i=is,ie
         k_dot_x = kx*xv(i) + ky*yv(j)
         den = 1.0d0 + amp*dsin(k_dot_x)
         vel = dsqrt(gam)*amp*dsin(k_dot_x)
         pre = 1.0d0/gam + amp*dsin(k_dot_x)

!         print*, i,xv(i), U(i,IDEN), U(i,IMX), U(i,IENE)

!      if( xv(i) <0.5d0 ) then
!      if( yv(j) <0.5d0 ) then
!          den = 1.0d0
!          vel = 0.0d0
!          pre = 1.0d0
!      else 
!          den = 0.125d0
!          vel = 0.0d0
!          pre = 0.1d0
!      endif

         U(i,j,IDEN) = den
!         U(i,j,IMX ) = den*vel
!         U(i,j,IMY ) = 0.0d0
         U(i,j,IMX ) = den*vel*kx/dsqrt(kx**2+ky**2)
         U(i,j,IMY ) = den*vel*ky/dsqrt(kx**2+ky**2)
         U(i,j,IENE) = pre/(gam - 1.0d0) + 0.5d0*den*vel**2
    enddo
    enddo

end Subroutine

!---------------------------------------------------------
!---------------------------------------------------------

Subroutine Cons2Prim( U, W )
use globals
implicit none
real(8), intent(in ) :: U(:,:,:)
real(8), intent(out) :: W(:,:,:)
real(8) :: deni


    do j=js,je
    do i=is,ie
        W(i,j,IDEN) = U(i,j,IDEN)
        deni = 1.0/U(i,j,IDEN)

        W(i,j,IVX) = deni*U(i,j,IMX)
        W(i,j,IVY) = deni*U(i,j,IMY)
        W(i,j,IPRE) = ( U(i,j,IENE) - 0.5d0*( U(i,j,IMX)**2 + U(i,j,IMY)**2 )*deni )*(gam - 1.0d0)
    enddo
    enddo

end Subroutine

!---------------------------------------------------------
!---------------------------------------------------------

real(8) function NewDt( xf, yf, W )
use globals
implicit none
real(8), intent(in) :: xf(:), yf(:), W(:,:,:)

    dt = 1d100
    do j=js, je
    do i=is, ie
         dt = min( dt, min(xf(i+1) - xf(i),yf(j+1) - yf(j))/ & 
             (dsqrt(gam*W(i,j,IPRE)/W(i,j,IDEN)) + dabs(W(i,j,IVX)) + dabs(W(i,j,IVY)) ) )
    enddo
    enddo

    NewDt = CFL*dt

end function

!!---------------------------------------------------------
!!---------------------------------------------------------

Subroutine ApplyBoundaryCondition( W )
use globals
implicit none
integer :: i, j
real(8), intent(out) :: W(:,:,:)
real(8) :: deni


! left boundary 
    do j=js,je
    do i=1,ngh
!        W(is-i,j,IDEN) = W(is+i-1,j,IDEN)
!        W(is-i,j,IVX ) = W(is+i-1,j,IVX )
!        W(is-i,j,IVY ) = W(is+i-1,j,IVY )
!        W(is-i,j,IPRE) = W(is+i-1,j,IPRE)

        W(is-i,j,IDEN) = W(ie-i+1,j,IDEN)
        W(is-i,j,IVX ) = W(ie-i+1,j,IVX )
        W(is-i,j,IVY ) = W(ie-i+1,j,IVY )
        W(is-i,j,IPRE) = W(ie-i+1,j,IPRE)
    enddo
    enddo

! right boundary
    do j=js,je
    do i=1,ngh
!        W(ie+i,j,IDEN) = W(ie-i+1,j,IDEN)
!        W(ie+i,j,IVX ) = W(ie-i+1,j,IVX )
!        W(ie+i,j,IVY ) = W(ie-i+1,j,IVY )
!        W(ie+i,j,IPRE) = W(ie-i+1,j,IPRE)

        W(ie+i,j,IDEN) = W(is+i-1,j,IDEN)
        W(ie+i,j,IVX ) = W(is+i-1,j,IVX )
        W(ie+i,j,IVY ) = W(is+i-1,j,IVY )
        W(ie+i,j,IPRE) = W(is+i-1,j,IPRE)
    enddo
    enddo


! bottom boundary 
    do j=1,ngh
    do i=is,ie
!        W(i,js-j,IDEN) = W(i,js+j-1,IDEN)
!        W(i,js-j,IVX ) = W(i,js+j-1,IVX )
!        W(i,js-j,IVY ) = W(i,js+j-1,IVY )
!        W(i,js-j,IPRE) = W(i,js+j-1,IPRE)

        W(i,js-j,IDEN) = W(i,je-j+1,IDEN)
        W(i,js-j,IVX ) = W(i,je-j+1,IVX )
        W(i,js-j,IVY ) = W(i,je-j+1,IVY )
        W(i,js-j,IPRE) = W(i,je-j+1,IPRE)
    enddo
    enddo

! top boundary
    do j=1,ngh
    do i=is,ie
!        W(i,je+j,IDEN) = W(i,je-j+1,IDEN)
!        W(i,je+j,IVX ) = W(i,je-j+1,IVX )
!        W(i,je+j,IVY ) = W(i,je-j+1,IVY )
!        W(i,je+j,IPRE) = W(i,je-j+1,IPRE)

        W(i,je+j,IDEN) = W(i,js+j-1,IDEN)
        W(i,je+j,IVX ) = W(i,js+j-1,IVX )
        W(i,je+j,IVY ) = W(i,js+j-1,IVY )
        W(i,je+j,IPRE) = W(i,js+j-1,IPRE)
    enddo
    enddo

end Subroutine

!---------------------------------------------------------
!---------------------------------------------------------
Subroutine CalculateNumericalFlux( xf, yf, xv, yv, W, Fx, Fy )
use globals
implicit none
real(8), intent(in ) :: xf(:), xv(:), yf(:), yv(:), W(:,:,:)
real(8), intent(out) :: Fx(:,:,:), Fy(:,:,:)

integer :: i, n
real(8), allocatable :: WL(:,:,:), WR(:,:,:)
real(8) :: dwL, dwR, dW, dwc

allocate( WL(nxtot,nytot,NVARS), WR(nxtot,nytot,NVARS) )

!---( 1st order )---
do n=1, NVARS
   do j=js, je
   do i=is-1, ie+1

        dwL = W(i,j,n) - W(i-1,j,n)
        dwR = W(i+1,j,n) - W(i,j,n)
        dwC = 0.5d0*( dwL + dwR )
        if (dwL*dwR > 0) then
!            dW = sign(dwL)*min( 2.0*dabs(dwL), 2.0*dabs(dwR), dabs(dwC))
            dW = 2.0d0*dwL*dwR/(dwL+dwR+1d-20)
        else
            dW = 0.0d0
        endif
       WL(i,j,n) = W(i,j,n) + 0.5d0*dW
       WR(i-1,j,n) = W(i,j,n) - 0.5d0*dW

!       WL(i,j,n) = W(i,j,n)
!       WR(i-1,j,n) = W(i,j,n)
   enddo
   enddo
enddo

 call HLLNumericalFluxX( is-1, ie, js, je, WL, WR, Fx )

!---( 1st order )---
do n=1, NVARS
   do j=js-1, je+1
   do i=is, ie
!       WL(i,j,n) = W(i,j,n)
!       WR(i,j,n) = W(i,j+1,n)
        dwL = W(i,j,n) - W(i,j-1,n)
        dwR = W(i,j+1,n) - W(i,j,n)
        dwC = 0.5d0*( dwL + dwR )
        if (dwL*dwR > 0) then
!            dW = sign(dwL)*min( 2.0*dabs(dwL), 2.0*dabs(dwR), dabs(dwC))
            dW = 2.0d0*dwL*dwR/(dwL+dwR+1d-20)
!            drho[i] = 2.0*dwL*dwR/(dwL+dwR+1e-8)
        else
            dW = 0.0d0
        endif
        WL(i,j   ,n) = W(i,j,n) + 0.5d0*dW
        WR(i,j-1 ,n) = W(i,j,n) - 0.5d0*dW

!       WL(i,j,n) = W(i,j,n)
!       WR(i,j-1,n) = W(i,j,n)
   enddo
   enddo
enddo

 call HLLNumericalFluxY( is, ie, js-1, je, WL, WR, Fy )

end Subroutine
!
!Subroutine GodunovNumericalFlux( iss, iee, WL, WR, F )
!use globals
!implicit none
!integer, intent(in) :: iss, iee
!real(8), intent(in ) :: WL(:,:), WR(:,:)
!real(8), intent(out) :: F(:,:)
!real(8) :: wwL(NVARS), wwR(NVARS)
!
!integer :: i
!real(8) :: prest(nxtot), velst(nxtot)
!real(8) :: aL, aR, S(3), rhostL, rhostR
!
!     call ExactRiemanSolver(iss, iee, WL, WR, velst, prest)
!
!     do i=iss, iee
!        wwL(:) = WL(i,:)
!        wwR(:) = WR(i,:)
!
!        aL = dsqrt(gam*wwL(IPRE)/wwL(IDEN))
!        aR = dsqrt(gam*wwR(IPRE)/wwR(IDEN))
!
!        if ( prest(i) >= wwL(IPRE) ) then !shock
!            rhostL = wwL(IDEN)*(prest(i)/wwL(IPRE) + gamm1/gamp1)/( gamm1/gamp1*prest(i)/wwL(IPRE) + 1.0d0)
!            S(1) = wwL(IVX) - aL*dsqrt( ( gamp1*prest(i)/wwL(IPRE) + gamm1 )/(2.0d0*gam) )
!        else !rarefuction 
!            rhostL = wwL(IDEN)*(prest(i)/wwL(IPRE))**(1.0d0/gam)
!            S(1) = 0.5d0*( wwL(IVX) - aL + velst(i) - dsqrt(gam*prest(i)/rhostL) )
!        endif
!
!        S(2) = velst(i) !const dicontinuity
!
!        if(prest(i) >= wwR(IPRE)) then !shock
!            rhostR = wwR(IDEN)*(prest(i)/wwR(IPRE) + gamm1/gamp1)/( gamm1/gamp1*prest(i)/wwR(IPRE) + 1.0d0)
!            S(3) = wwR(IVX) + aR*dsqrt( ( gamp1*prest(i)/wwR(IPRE) + gamm1 )/(2.0d0*gam) )
!        else !rarefuction
!            rhostR = wwR(IDEN)*(prest(i)/wwR(IPRE))**(1.0d0/gam)
!            S(3) = 0.5d0*( wwR(IVX) + aR + velst(i) + dsqrt(gam*prest(i)/rhostR) )
!        endif
!
!        if( S(1) >= 0) then ! W_L
!            F(i,IDEN) = wwL(IDEN)*wwL(IVX)
!            F(i,IMX ) = wwL(IPRE) + wwL(IDEN)*wwL(IVX)**2
!            F(i,IENE) = ( gam*wwL(IPRE)/gamm1 + 0.5d0*wwL(IDEN)*wwL(IVX)**2 )*wwL(IVX)
!        else if (S(3) <= 0) then ! W_R
!            F(i,IDEN) = wwR(IDEN)*wwR(IVX)
!            F(i,IMX ) = wwR(IPRE) + wwR(IDEN)*wwR(IVX)**2
!            F(i,IENE) = ( gam*wwR(IPRE)/gamm1 + 0.5d0*wwR(IDEN)*wwR(IVX)**2 )*wwR(IVX)
!        else if (S(2) >= 0) then !W_L*
!            F(i,IDEN) = rhostL*velst(i)
!            F(i,IMX)  = prest(i) + rhostL*velst(i)**2
!            F(i,IENE) = ( gam*prest(i)/gamm1 + 0.5d0*rhostL*velst(i)**2 )*velst(i)
!        else  !W_R*
!            F(i,IDEN) = rhostR*velst(i)
!            F(i,IMX)  = prest(i) + rhostR*velst(i)**2
!            F(i,IENE) = ( gam*prest(i)/gamm1 + 0.5d0*rhostR*velst(i)**2 )*velst(i)
!        endif
!     enddo
!
!
!end Subroutine
!
Subroutine HLLNumericalFluxX( iss, iee, jss, jee, WL, WR, Fx )
use globals
implicit none
integer, intent(in) :: iss, iee, jss, jee
real(8), intent(in ) :: WL(:,:,:), WR(:,:,:)
real(8), intent(out) :: Fx(:,:,:)
real(8) :: wwL(NVARS), wwR(NVARS)

integer :: i
real(8) :: cL, cR, eneL, eneR
real(8) :: FrhoL, FmoxL, FmoyL, FeneL
real(8) :: FrhoR, FmoxR, FmoyR, FeneR, SL, SR

     do j=jss, jee
     do i=iss, iee
        wwL(:) = WL(i,j,:)
        wwR(:) = WR(i,j,:)

        cL = dsqrt(gam*wwL(IPRE)/wwL(IDEN))
        cR = dsqrt(gam*wwL(IPRE)/wwL(IDEN))

        SL = min(wwL(IVX) - cL, wwR(IVX) - cR)
        SR = max(wwL(IVX) + cL, wwR(IVX) + cR)

        eneL = 0.5*wwL(IDEN)*( wwL(IVX)**2 + wwL(IVY)**2 ) + wwL(IPRE)/(gam - 1.0)
        eneR = 0.5*wwR(IDEN)*( wwR(IVX)**2 + wwR(IVY)**2 ) + wwR(IPRE)/(gam - 1.0)

        FrhoL = wwL(IDEN)*wwL(IVX)
        FmoxL = wwL(IPRE) + wwL(IDEN)*wwL(IVX)**2
        FmoyL = wwL(IDEN)*wwL(IVX)*wwL(IVY)
        FeneL = ( wwL(IPRE) + eneL )*wwL(IVX)

        FrhoR = wwR(IDEN)*wwR(IVX)
        FmoxR = wwR(IPRE) + wwR(IDEN)*wwR(IVX)**2
        FmoyR =  wwR(IDEN)*wwR(IVX)*wwR(IVY)
        FeneR = ( wwR(IPRE) + eneR )*wwR(IVX)


        if( SL > 0) then ! W_L
            Fx(i,j,IDEN) = FrhoL
            Fx(i,j,IMX ) = FmoxL
            Fx(i,j,IMY ) = FmoyL
            Fx(i,j,IENE) = FeneL
        else if (SR <= 0) then ! W_R
            Fx(i,j,IDEN) = FrhoR
            Fx(i,j,IMX ) = FmoxR
            Fx(i,j,IMY ) = FmoyR
            Fx(i,j,IENE) = FeneR
        else 
            Fx(i,j,IDEN) = (SR*FrhoL - SL*FrhoR + SR*SL*(wwR(IDEN) - wwL(IDEN)))/(SR - SL)
            Fx(i,j,IVX) = (SR*FmoxL - SL*FmoxR + SR*SL*(wwR(IDEN)*wwR(IVX) - wwL(IDEN)*wwL(IVX)) )/(SR - SL)
            Fx(i,j,IVY) = (SR*FmoyL - SL*FmoyR + SR*SL*(wwR(IDEN)*wwR(IVY) - wwL(IDEN)*wwL(IVY)) )/(SR - SL)
            Fx(i,j,IENE) = (SR*FeneL - SL*FeneR + SR*SL*(eneR - eneL))/(SR - SL)
        endif
     enddo
     enddo


end Subroutine
!
Subroutine HLLNumericalFluxY( iss, iee, jss, jee, WL, WR, Fy )
use globals
implicit none
integer, intent(in) :: iss, iee, jss, jee
real(8), intent(in ) :: WL(:,:,:), WR(:,:,:)
real(8), intent(out) :: Fy(:,:,:)
real(8) :: wwL(NVARS), wwR(NVARS)

integer :: i
real(8) :: cL, cR, eneL, eneR
real(8) :: FrhoL, FmoxL, FmoyL, FeneL
real(8) :: FrhoR, FmoxR, FmoyR, FeneR, SL, SR

     do j=jss, jee
     do i=iss, iee
        wwL(:) = WL(i,j,:)
        wwR(:) = WR(i,j,:)

        cL = dsqrt(gam*wwL(IPRE)/wwL(IDEN))
        cR = dsqrt(gam*wwL(IPRE)/wwL(IDEN))

        SL = min(wwL(IVY) - cL, wwR(IVY) - cR)
        SR = max(wwL(IVY) + cL, wwR(IVY) + cR)

        eneL = 0.5*wwL(IDEN)*( wwL(IVX)**2 + wwL(IVY)**2 ) + wwL(IPRE)/(gam - 1.0)
        eneR = 0.5*wwR(IDEN)*( wwR(IVX)**2 + wwR(IVY)**2 ) + wwR(IPRE)/(gam - 1.0)

        FrhoL = wwL(IDEN)*wwL(IVY)
        FmoxL = wwL(IDEN)*wwL(IVX)*wwL(IVY)
        FmoyL = wwL(IPRE) + wwL(IDEN)*wwL(IVY)**2
        FeneL = ( wwL(IPRE) + eneL )*wwL(IVY)

        FrhoR = wwR(IDEN)*wwR(IVY)
        FmoxR =  wwR(IDEN)*wwR(IVX)*wwR(IVY)
        FmoyR = wwR(IPRE) + wwR(IDEN)*wwR(IVY)**2
        FeneR = ( wwR(IPRE) + eneR )*wwR(IVY)

        if( SL > 0) then ! W_L
            Fy(i,j,IDEN) = FrhoL
            Fy(i,j,IMX ) = FmoxL
            Fy(i,j,IMY ) = FmoyL
            Fy(i,j,IENE) = FeneL
        else if (SR <= 0) then ! W_R
            Fy(i,j,IDEN) = FrhoR
            Fy(i,j,IMX ) = FmoxR
            Fy(i,j,IMY ) = FmoyR
            Fy(i,j,IENE) = FeneR
        else 
            Fy(i,j,IDEN) = (SR*FrhoL - SL*FrhoR + SR*SL*(wwR(IDEN) - wwL(IDEN)))/(SR - SL)
            Fy(i,j,IVX) = (SR*FmoxL - SL*FmoxR + SR*SL*(wwR(IDEN)*wwR(IVX) - wwL(IDEN)*wwL(IVX)) )/(SR - SL)
            Fy(i,j,IVY) = (SR*FmoyL - SL*FmoyR + SR*SL*(wwR(IDEN)*wwR(IVY) - wwL(IDEN)*wwL(IVY)) )/(SR - SL)
            Fy(i,j,IENE) = (SR*FeneL - SL*FeneR + SR*SL*(eneR - eneL))/(SR - SL)
        endif
     enddo
     enddo


end Subroutine
!
!
Subroutine UpdateConservedVariables( dt, xf, yf, Uold, Fx, Fy, U )
use globals
implicit none
integer i, j, n
real(8), intent(in) :: dt, Fx(:,:,:), Fy(:,:,:)
real(8), intent(in) :: xf(:), yf(:)
real(8), intent(inout) :: Uold(:,:,:)
real(8), intent(out) :: U(:,:,:)

    do n=1,NVARS
        do j=js,je
        do i=is,ie
             U(i,j,n) = Uold(i,j,n) - dt*( Fx(i,j,n) - Fx(i-1,j,n) )/(xf(i) - xf(i-1)) &
                                    - dt*( Fy(i,j,n) - Fy(i,j-1,n) )/(yf(j) - yf(j-1)) 
        enddo
        enddo
    enddo

end Subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subroutine ExactRiemanSolver(iss, iee, wL, wR, velst, prest)
!use globals
!implicit none
!integer, intent(in) :: iss, iee
!real(8), intent(in) :: wL(:,:), wR(:,:)
!real(8), intent(out) :: velst(:), prest(:)
!real(8) :: wwL(NVARS), wwR(NVARS)
!real(8) :: cL, cR, zR, zL, vR, vL, mL, mR, msqL, msqR
!real(8) :: testP, hantei
!
!do i=iss, iee
!    wwL(:) = wL(i,:)
!    wwR(:) = wR(i,:)
!
!    wwL(IDEN) = 1.0d0
!    wwL(IVX) = 0.0d0
!    wwL(IPRE) = 1.0d0
!    wwR(IDEN) = 0.125d0
!    wwR(IVX) = 0.0d0
!    wwR(IPRE) = 0.1d0
!
!    cL = dsqrt(gam*wwL(IPRE)*wwL(IDEN)) 
!    cR = dsqrt(gam*wwR(IPRE)*wwR(IDEN))
!
!    testP = ( cR*wwL(IPRE) + cL*wwR(IPRE) - cR*cL*(wwL(IVX) - wwR(IVX)) )/(cR + cL)
!    if(testP<0) then
!        testP = 1d-8
!    endif
!
!    zR = 0.0d0
!    zL = 0.0d0
!    vR = 0.0d0
!    vL = 0.0d0
!    hantei = 100.0d0
!    do while( abs(hantei) > 1d-4 )
!        if( testP >= wwL(IPRE) ) then !shock
!            msqL = 0.5d0*cL**2*(gamp1*testP/wwL(IPRE) + gamm1)/gam 
!            mL = sqrt(msqL) 
!            zL = 2.0d0*msqL*mL/(msqL + cL**2)
!        else  !rarefaction
!            mL = gamm1/(2.0*gam)*(1.0 - testP/wwL(IPRE))/(1.0 - (testP/wwL(IPRE))**(gamm1/(2.0*gam)) )*cL
!            zL = cL*(testP/wwL(IPRE))**(1.0 - gamm1/(2.0*gam))
!        endif
!        if( testP >= wwR(IPRE) ) then !shock
!            msqR = 0.5d0*cR**2*(gamp1*testP/wwR(IPRE) + gamm1)/gam 
!            mR = sqrt(msqR) 
!            zR = 2.0d0*msqR*mR/(msqR + cR**2)
!        else ! rarefaction
!            mR = gamm1/(2.0*gam)*(1.0 - testP/wwR(IPRE))/(1.0 - (testP/wwR(IPRE))**(gamm1/(2.0*gam)) )*cR
!            zR = cR*(testP/wwR(IPRE))**(1.0 - gamm1/(2.0*gam))
!        endif
!
!
!        vR = wwR(IVX) + (testP - wwR(IPRE))/(mR)
!        vL = wwL(IVX) - (testP - wwL(IPRE))/(mL)
!
!
!        hantei = zR*zL*( vR - vL )/( (zR + zL)*testP )
!        testP = testP*(1.0 - hantei)
!    enddo
!
!    prest(i) = testP
!    velst(i) = (zL*vL + zR*vR)/(zL + zR)
!enddo
!
!
!end Subroutine
!
end program
!
!
!!import numpy as np
!!import os
!!
!!import matplotlib.pyplot as plt
!!import matplotlib.tri as tri
!!
!!!preparation for the animation
!!from matplotlib.animation import ArtistAnimation
!!from IPython.display import HTML
!!from IPython.display import Image
!!
!!def InitialProfile(x): 
!!    global x0, Lx
!!
!!    if x>x0*0.8 and x<x0*1.2: 
!!        return 1.0 
!!    else: 
!!        return 0.0
!!
!!def main():
!!    global Lx, ngrid, nxtot, c, x0, flag_legend, nghost, ien, ist, gam
!!
!!    gam = 1.4
!!    tout  = 0.25    ! the final time 
!!
!!    rhoL = 1.0
!!    velL = 0.0
!!    preL = 1.0
!!
!!    rhoR = 0.125
!!    velR = 0.0
!!    preR = 0.1
!!
!!
!!    Lx    = 1.0     ! the box size [0,Lx]
!!    ngrid = 64      ! the total numer of the grids
!!    x0    = 0.5*Lx ! position of the initial wave front
!!    dt_plot = 0.01  ! time interval of plot
!!    CFL = 0.8       ! CFL number
!!
!!    nghost = 2      ! the number of ghost cells
!!    nxtot = ngrid + 2*nghost ! total cell number including the ghost cells
!!    ist = nghost;   ! the first index of the computational domain
!!    ien = nghost+ngrid; ! the last index of the computational domain
!!
!!    gam1 = gam - 1.0
!!
!!    !definition of variables
!!    xv = np.zeros(nxtot)   ! cell center   xv[i] --> x_i
!!    xf = np.zeros(nxtot+1) ! cell boundary xf[i] --> x_{i+1/2}
!!    rho = np.zeros(nxtot)    ! u^{n+1} u[i] --> u_i^{n+1}
!!    vel = np.zeros(nxtot)   ! u^{n+1} u[i] --> u_i^{n+1}
!!    pre = np.zeros(nxtot) ! u^n     uold[i] --> u_i^n
!!    mom = np.zeros(nxtot)   ! u^{n+1} u[i] --> u_i^{n+1}
!!    ene = np.zeros(nxtot) ! u^n     uold[i] --> u_i^n
!!    rho0 = np.zeros(nxtot)    ! u^{n+1} u[i] --> u_i^{n+1}
!!    mom0 = np.zeros(nxtot)   ! u^{n+1} u[i] --> u_i^{n+1}
!!    ene0 = np.zeros(nxtot) ! u^n     uold[i] --> u_i^n
!!    Frho = np.zeros(nxtot) ! u^n     uold[i] --> u_i^n
!!    Fmom = np.zeros(nxtot) ! u^n     uold[i] --> u_i^n
!!    Fene = np.zeros(nxtot) ! u^n     uold[i] --> u_i^n
!!    u = np.zeros(nxtot) ! u^n     uold[i] --> u_i^n
!!
!!    dx = Lx/ngrid
!!    sod_ana = np.loadtxt("sod_t0.25.dat")
!!
!!    PreparationAnimation()
!!    
!!    new_dir_path = 'outputs_hydro_FVM/'
!!    if not os.path.exists(new_dir_path):
!!        os.mkdir(new_dir_path)
!!    
!!    !cell boundary position
!!    for i in range(nxtot+1):
!!        xf[i] = (i - ist)*dx
!!
!!    !cell center position
!!    for i in range(nxtot):
!!        xv[i] = 0.5*(xf[i] + xf[i+1])
!!
!!    !initial profile 
!!    for i in range(nxtot):
!!        rho[i] = (rhoR - rhoL)*0.5*( np.tanh((xv[i]-0.5*Lx)/(0.5*dx)) + 1.0) + rhoL
!!        vel[i] = (velR - velL)*0.5*( np.tanh((xv[i]-0.5*Lx)/(0.5*dx)) + 1.0) + velL
!!        pre[i] = (preR - preL)*0.5*( np.tanh((xv[i]-0.5*Lx)/(0.5*dx)) + 1.0) + preL
!!
!!        mom[i] = rho[i]*vel[i]
!!        ene[i] = pre[i]/gam1 + 0.5*rho[i]*vel[i]**2
!!    
!!    time = 0
!!    time_plot = 0
!!
!!    dt = CFL*np.min(dx/(np.sqrt(gam*pre[ist:ien]/rho[ist:ien]) + abs(vel[ist:ien])))
!!    
!!    flag_legend=True
!!    while time < tout: 
!!
!!        ! substep
!!        for i in range(ist,ien): 
!!            rho0[i] = rho[i] 
!!            mom0[i] = mom[i] 
!!            ene0[i] = ene[i] 
!!        ! prediction
!!        CalculateRoeNumericalFlux(rho, vel, pre, Frho, Fmom, Fene)
!!        UpdateConservedVariables(0.5*dt, dx, rho, rho, mom, mom, ene, ene, Frho, Fmom, Fene)
!!        ConvertConserved2Primitive( rho, mom, ene, vel, pre )
!!        SetBoundaryCondition(xv,rho,vel,pre)
!!    
!!        ! correction
!!        CalculateRoeNumericalFlux(rho, vel, pre, Frho, Fmom, Fene)
!!        UpdateConservedVariables(dt, dx, rho, rho0, mom, mom0, ene, ene0, Frho, Fmom, Fene)
!!        ConvertConserved2Primitive( rho, mom, ene, vel, pre )
!!        SetBoundaryCondition(xv,rho,vel,pre)
!!
!!        time = time + dt
!!        time_plot = time_plot + dt
!!    
!!        if time_plot > dt_plot:
!!            AddPictureToAnimation(time,xv,rho,vel,pre,sod_ana)
!!            time_plot = 0.0
!!   
!!        if  time + dt > tout:
!!            dt = tout - time;
!!
!!        dt = CFL*np.min(dx/(np.sqrt(gam*pre[ist:ien]/rho[ist:ien]) + abs(vel[ist:ien])))
!!    
!!    for i in range(ist,ien):
!!        print xv[i], rho[i], vel[i], pre[i] 
!!    animationA = ArtistAnimation(figA, artistA, interval=100)
!!    animationA.save('animate_hydro_FVM_2nd.mp4',writer='ffmpeg')
!!    HTML(animationA.to_html5_video())
!!
!!def CalculateRoeNumericalFlux(rho, vel, pre, Frho, Fmom, Fene): 
!!    global ien, ist, gam, nxtot
!!
!!    gam1 = gam - 1.0
!!
!!    rhol = np.zeros(nxtot)   ! rhoL(i+1/2)
!!    rhor = np.zeros(nxtot)   ! rhoR(i+1/2)
!!    vell = np.zeros(nxtot)   ! velL(i+1/2)
!!    velr = np.zeros(nxtot)   ! velR(i+1/2)
!!    prel = np.zeros(nxtot)   ! preL(i+1/2)
!!    prer = np.zeros(nxtot)   ! preR(i+1/2)
!!
!!    ! calculate w_{i+1/2,L}
!!    for i in range(ist-1,ien): 
!!        rhol[i] = rho[i] 
!!        vell[i] = vel[i] 
!!        prel[i] = pre[i] 
!!
!!    ! calculate u_{i+1/2,R}
!!    for i in range(ist-1,ien): 
!!        rhor[i] = rho[i+1] 
!!        velr[i] = vel[i+1]
!!        prer[i] = pre[i+1] 
!!
!!    GodunovNumericalFlux(ist,ien,gam, rhol, vell, prel, rhor, velr, prer, Frho, Fmom, Fene)
!!
!!def UpdateConservedVariables(dt, dx, rho, rho0, mom, mom0, ene, ene0, Frho, Fmom, Fene):
!!    global ien, ist, gam, nxtot
!!
!!    nu = dt/dx
!!    for i in range(ist,ien): 
!!        rho[i] = rho0[i] - nu*( Frho[i] - Frho[i-1] ) 
!!        mom[i] = mom0[i] - nu*( Fmom[i] - Fmom[i-1] ) 
!!        ene[i] = ene0[i] - nu*( Fene[i] - Fene[i-1] )
!!
!!def ConvertConserved2Primitive( rho, mom, ene, vel, pre ): 
!!    global gam
!!
!!    gam1 = gam - 1
!!    ! conver the conserved variables to the primitive variables 
!!    for i in range(ist,ien): 
!!        vel[i] = mom[i]/rho[i] 
!!        pre[i] = gam1*(ene[i] - 0.5*rho[i]*vel[i]**2)
!!
!!def SetBoundaryCondition(x,rho,vel,pre):
!!    global nghost, ien, ist
!!    
!!!    for i in range(nghost): 
!!!        rho[i] = rho[ien+i-nghost]
!!!        vel[i] = vel[ien+i-nghost]
!!!        pre[i] = pre[ien+i-nghost]
!!!
!!!    for i in range(nghost): 
!!!        rho[ien+i] = rho[ist+i]
!!!        vel[ien+i] = vel[ist+i]
!!!        pre[ien+i] = pre[ist+i]
!!
!!def ExactRiemanSolver(gam, rhol, vell, prel, rhor, velr, prer ):
!!    gamp1 = gam + 1.0
!!    gamm1 = gam - 1.0
!!
!!    cL = np.sqrt(gam*prel*rhol) 
!!    cR = np.sqrt(gam*prer*rhor)
!!
!!    testP = ( cR*prel + cL*prer - cR*cL*(velr - vell) )/(cR + cL)
!!    if testP<0: testP = 1e-8
!!    hantei = 10
!!    zR = 0.0
!!    zL = 0.0
!!    vR = 0.0
!!    vL = 0.0
!!    while abs(hantei) > 1e-4: 
!!        wsqL = 0.5*cL**2*(gamp1*testP/prel + gamm1)/gam 
!!        wL = np.sqrt(wsqL) 
!!        zL = 2.0*wsqL*wL/(wsqL + cL**2)
!!
!!        wsqR = 0.5*cR**2*(gamp1*testP/prer + gamm1)/gam 
!!        wR = np.sqrt(wsqR) 
!!        zR = 2.0*wsqR*wR/(wsqR + cR**2)
!!
!!        vR = velr + (testP - prer)/(wR)
!!        vL = vell - (testP - prel)/(wL)
!!
!!        hantei = zR*zL*( vR - vL )/( (zR + zL)*testP )
!!        testP = testP*(1.0 - hantei)
!!
!!    return testP, (zL*vL + zR*vR)/(zL + zR)
!!
!!def GodunovNumericalFlux(ist,ien,gam, rhol, vell, prel, rhor, velr, prer, Frho, Fmom, Fene):
!!
!!    for i in range(ist-1,ien): 
!!        gamp1 = gam + 1.0
!!        gamm1 = gam - 1.0
!!    
!!        prest, velst = ExactRiemanSolver(gam, rhol[i], vell[i], prel[i], rhor[i], velr[i], prer[i] )
!!    
!!        aL = np.sqrt(gam*prel[i]/rhol[i])
!!        aR = np.sqrt(gam*prer[i]/rhor[i])
!!    
!!        S = np.zeros(3)
!!    
!!        if prest >= prel[i]: !shock
!!            rhostL = rhol[i]*(prest/prel[i] + gamm1/gamp1)/( gamm1/gamp1*prest/prel[i] + 1.0)
!!            S[0] = vell[i] - aL*np.sqrt( ( gamp1*prest/prel[i] + gamm1 )/(2.0*gam) )
!!        else: !rarefuction 
!!            rhostL = rhol[i]*(prest/prel[i])**(1.0/gam)
!!            S[0] = 0.5*( vell[i] - aL + velst - np.sqrt(gam*prest/rhostL) )
!!        
!!        S[1] = velst !const dicontinuity
!!    
!!        if prest >= prer[i]: !shock
!!            rhostR = rhor[i]*(prest/prer[i] + gamm1/gamp1)/( gamm1/gamp1*prest/prer[i] + 1.0)
!!            S[2] = velr[i] + aR*np.sqrt( ( gamp1*prest/prer[i] + gamm1 )/(2.0*gam) )
!!        else: !rarefuction
!!            rhostR = rhor[i]*(prest/prer[i])**(1.0/gam)
!!            S[2] = 0.5*( velr[i] + aR + velst + np.sqrt(gam*prest/rhostR) )
!!
!!        if S[0] >= 0: ! W_L
!!            Frho[i] = rhol[i]*vell[i]
!!            Fmom[i] = prel[i] + rhol[i]*vell[i]**2
!!            Fene[i] = ( gam*prel[i]/(gam-1.0) + 0.5*rhol[i]*vell[i]**2 )*vell[i]
!!        elif S[2] <= 0: ! W_R
!!            Frho[i] = rhor[i]*velr[i]
!!            Fmom[i] = prer[i] + rhor[i]*velr[i]**2
!!            Fene[i] = ( gam*prer[i]/(gam-1.0) + 0.5*rhor[i]*velr[i]**2 )*velr[i]
!!        elif S[1] >= 0: ! W_L*
!!            Frho[i] = rhostL*velst
!!            Fmom[i] = prest + rhostL*velst**2
!!            Fene[i] = ( gam*prest/(gam-1.0) + 0.5*rhostL*velst**2 )*velst
!!        else:
!!            Frho[i] = rhostR*velst
!!            Fmom[i] = prest + rhostR*velst**2
!!            Fene[i] = ( gam*prest/(gam-1.0) + 0.5*rhostR*velst**2 )*velst
!!
!!            
!!
!!         
!!        
!!    
!!
!!
!!def PreparationAnimation():
!!    global figA, artistA, ax1, Lx, c, x0
!!    !set up for animation
!!    figA = plt.figure(figsize=(8,6))
!!    plt.rcParams['font.size']=20
!!    artistA = []
!!    ax1 = figA.add_subplot(1,1,1)
!!
!!    !set up for animation
!!    figA = plt.figure(figsize=(8,6))
!!    plt.rcParams['font.size']=20
!!    artistA = []
!!    ax1 = figA.add_subplot(1,1,1)
!!
!!def AddPictureToAnimation(time,xv,rho,vel,pre,sod_ana): 
!!    global ax1, nxtot, artistA, Lx, x0, flag_legend, gam
!!    timetxt='time=%.2f'%(time)
!!    ax1.set_xlim(0,Lx)
!!    ax1.set_ylim(0,1.1)
!!    ax1.set_xlabel(r'$x$')
!!    title = ax1.text(0.5*Lx,2.1,timetxt,horizontalalignment="center")
!!    im1 = ax1.plot(xv,rho,'o-',color='red',label=r"$\rho$")
!!    im2 = ax1.plot(xv,vel,'o-',color='green',label=r"$v$")
!!    im3 = ax1.plot(xv,pre,'o-',color='blue',label=r"$P$")
!!
!!    im4 = ax1.plot((sod_ana[:,0]-0.5)*time/0.25+0.5,sod_ana[:,1],'-',color='black')
!!    im5 = ax1.plot((sod_ana[:,0]-0.5)*time/0.25+0.5,sod_ana[:,2],'-',color='black')
!!    im6 = ax1.plot((sod_ana[:,0]-0.5)*time/0.25+0.5,sod_ana[:,3],'-',color='black')
!!
!!!    im2 = ax1.plot(xv,usol,'-',color='blue',label="analytic solution")
!!    if flag_legend: 
!!        ax1.legend(loc='upper right') 
!!        flag_legend = False
!!    
!!    figA.tight_layout()
!!    artistA.append(im1+im2+im3+im4+im5+im6+[title])
!!
!!if __name__ == "__main__":
!!    main()
