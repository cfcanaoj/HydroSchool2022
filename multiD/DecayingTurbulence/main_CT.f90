program main
implicit none

integer::ntime = 0 ! counter of the timestep
integer,parameter::ntimemax=200000     ! the maximum timesteps
real(8)::time  = 0.0d0  ! time 
real(8),parameter:: timemax=20.0d0
real(8),parameter:: dtout=0.2d0

integer, parameter :: flag_flux = 3 ! 1 (HLL), 2 (HLLC), 3 (HLLD)

integer,parameter::nx=128*1*4 ! the number of grids in the simulation box
integer,parameter::ny=128*2*4 ! the number of grids in the simulation box
integer,parameter::nz=1          ! the number of grids in the simulation box
integer,parameter::mgn=2         ! the number of ghost cells
integer,parameter::nxtot=nx+2*mgn+1 ! the total number of grids including ghost cells
integer,parameter::nytot=ny+2*mgn+1 ! the total number of grids including ghost cells
integer,parameter::nztot=2 ! the total number of grids including ghost cells
integer,parameter::is=mgn+1         ! the index of the leftmost grid
integer,parameter::js=mgn+1         ! the index of the leftmost grid
integer,parameter::ks=1         ! the index of the leftmost grid
integer,parameter::ie=nx+mgn     ! the index of the rightmost grid
integer,parameter::je=ny+mgn     ! the index of the rightmost grid
integer,parameter::ke=1 ! the index of the rightmost grid
real(8),parameter::x1min=-0.5d0,x1max=0.5d0
real(8),parameter::x2min=-1.0d0,x2max=1.0d0
real(8),parameter::x3min=0.0d0,x3max=1.0d0

real(8),parameter::Ccfl=0.2d0


integer, parameter :: IDN = 1
integer, parameter :: IM1 = 2
integer, parameter :: IM2 = 3
integer, parameter :: IM3 = 4
integer, parameter :: IPR = 5
integer, parameter :: ISC = 6
integer, parameter :: IB1 = 7
integer, parameter :: IB2 = 8
integer, parameter :: IB3 = 9
integer, parameter :: NVAR = 6
integer, parameter :: NFLX = 9

integer, parameter :: IVX = 2
integer, parameter :: IVY = 3
integer, parameter :: IVZ = 4
integer, parameter :: IEN = 5

real(8),parameter::gam=5.0d0/3.0d0 !! adiabatic index

integer, parameter :: nevo = 3
real(8) :: phys_evo(nevo)

real(8),dimension(nxtot)::xf,xv
real(8),dimension(nytot)::yf,yv
real(8),dimension(nztot)::zf,zv
real(8),dimension(nxtot,nytot,nztot,NVAR) :: Uo
real(8),dimension(nxtot,nytot,nztot,NVAR) :: U
real(8),dimension(nxtot,nytot,nztot,NVAR) :: Q
real(8),dimension(nxtot,nytot,nztot,3) :: Bso 
real(8),dimension(nxtot,nytot,nztot,3) :: Bs 
real(8),dimension(nxtot,nytot,nztot,3) :: Bc
real(8),dimension(nxtot,nytot,nztot,NVAR) :: F
real(8),dimension(nxtot,nytot,nztot,NVAR) :: G
real(8),dimension(nxtot,nytot,nztot,NVAR) :: H
real(8),dimension(nxtot,nytot,nztot,3) :: E 
real(8) :: dt
integer :: i,j,k


!      write(6,*) "setup grids and initial condition"
    call GenerateGrid(xf, xv, yf, yv, zf, zv)
    call GenerateProblem(xv, yv, zv, Q, Bs, Bc)
    call ConsvVariable(Q, Bc, U)
    call BoundaryCondition( Q, Bs, Bc )
!    call Output( .TRUE., xf, xv, yf, yv, Q, Bc )



    open(1,file="vy_evo_B0.8_ct_hll1.dat",action="write")
! main loop
    mloop: do 
        call TimestepControl(xf, yf, zf, Q, Bc, dt)
        if( time + dt > timemax ) dt = timemax - time

        Uo(:,:,:,:) = U(:,:,:,:)
        Bso(:,:,:,:) = Bs(:,:,:,:)


        call NumericalFlux( xf, yf, zf, Q, Bc, F, G, H, E )
        call UpdateConsv( 0.5d0*dt, xf, yf, zf, F, G, H, E, Q, U, Bs, U, Bs )
        call PrimVariable( U, Bs, Q, Bc )
        call BoundaryCondition( Q, Bs, Bc )

        call NumericalFlux( xf, yf, zf, Q, Bc, F, G, H, E )
        call UpdateConsv( dt, xf, yf, zf, F, G, H, E, Q, Uo, Bso, U, Bs )
        call PrimVariable( U, Bs, Q, Bc )
        call BoundaryCondition( Q, Bs, Bc )

         time=time+dt
         ntime = ntime + 1

         if( mod(ntime,10) .eq. 0 ) then
             call Analysis(xv,yv,Q,Bc,Bs,phys_evo)
             write(1,*) time, phys_evo(1:nevo)
             call flush(1)
         endif
!         call Output( .FALSE., xf, xv, yf, yv, Q, Bc)
!         call Output( .true., xf, xv, yf, yv, Q, Bc)

         print*, "time = ",time, "dt = ",dt

         if(time >= timemax) exit mloop
         if( ntime == 10) exit
      enddo mloop
      close(1)
!      call Output( .TRUE., xf, xv, yf, yv, Q, Bc)


!      write(6,*) "program has been finished"
contains 

    subroutine GenerateGrid(xf, xv, yf, yv, zf, zv)
    implicit none
    real(8), intent(out) :: xf(:), xv(:)
    real(8), intent(out) :: yf(:), yv(:)
    real(8), intent(out) :: zf(:), zv(:)
    real(8) :: dx,dy
    integer::i,j


    dx=(x1max-x1min)/dble(nx)
    do i=1,nxtot
         xf(i) = dx*(i-(mgn+1))+x1min
    enddo
    do i=1,nxtot-1
         xv(i) = 0.5d0*(xf(i+1)+xf(i))
    enddo

    dy=(x2max-x2min)/dble(ny)
    do j=1,nytot
         yf(j) = dx*(j-(mgn+1))+x2min
    enddo
    do j=1,nytot-1
         yv(j) = 0.5d0*(yf(j+1)+yf(j))
    enddo

    return
    end subroutine GenerateGrid

    subroutine GenerateProblem(xv, yv, zv, Q, Bs, Bc )
    implicit none
    integer::i, j, k
    real(8), intent(in ) :: xv(:), yv(:), zv(:)
    real(8), intent(out) :: Q(:,:,:,:)
    real(8), intent(out) :: Bs(:,:,:,:)
    real(8), intent(out) :: Bc(:,:,:,:)
    real(8) :: pi, den, B0, rho1, rho2, dv, wid, v1, v2, sig

    pi = dacos(-1.0d0)

    rho1 = 1.0d0
    rho2 = 1.0d0
    dv   = 2.00d0
    wid  = 0.05d0
    sig  = 0.2d0
!   B0  = 0.0d0*dsqrt( dv**2*0.5d0*rho1*rho2/(rho1+rho2))
    B0  = sqrt(2.0d0/3.0d0)

    do k=ks,ke
    do j=js,je
    do i=is,ie
        Q(i,j,k,IDN) = 1.0d0 !+ 0.5d0*( dtanh( (yv(j)+0.25d0)/wid ) - tanh( (yv(j)-0.25d0)/wid) )
        Q(i,j,k,IVX)  = 0.5*dv*( dtanh( (yv(j)+0.5d0)/wid ) - dtanh( (yv(j) - 0.5d0)/wid ) - 1.0d0 )
        Q(i,j,k,IVY)  = 0.001d0*dsin(2.0d0*pi*xv(i))* &
             ( dexp( - (yv(j) + 0.5d0)**2/sig**2 ) +  &
               dexp( - (yv(j) - 0.5d0)**2/sig**2 ) )
        Q(i,j,k,IPR) = 1.0d0

        Q(i,j,k,ISC) = 0.5d0*( dtanh( (yv(j)+0.5d0)/wid ) - tanh( (yv(j)-0.5d0)/wid) )
    enddo
    enddo
    enddo

    do k=ks,ke
    do j=js,je
    do i=is,ie+1
        Bs(i,j,k,1) = B0
    enddo
    enddo
    enddo

    do k=ks,ke
    do j=js,je+1
    do i=is,ie
        Bs(i,j,k,2) = 0.0d0
    enddo
    enddo
    enddo

    do k=ks,ke+1
    do j=js,je
    do i=is,ie
        Bs(i,j,k,3) = 0.0d0
    enddo
    enddo
    enddo

    call CellCenterMagneticField(is, ie, js, je, ks, ke, Bs, Bc)


    return
    end subroutine GenerateProblem

    subroutine BoundaryCondition( Q, Bs, Bc )
    implicit none
    real(8), intent(inout) :: Q(:,:,:,:)
    real(8), intent(inout) :: Bs(:,:,:,:)
    real(8), intent(inout) :: Bc(:,:,:,:)
    integer::i,j,k

    ! x inner boundary
    do k=ks,ke
    do j=js-mgn,je+mgn
    do i=1,mgn
        Q(is-i,j,k,:)  = Q(ie+1-i,j,k,:)
    enddo
    enddo
    enddo

    do k=ks,ke
    do j=js-mgn,je+mgn
    do i=1,mgn
          Bs(is-i,j,k,1) = Bs(ie+1-i,j,k,1)
    enddo
    enddo
    enddo

    do k=ks,ke
    do j=js-mgn,je+mgn+1
    do i=1,mgn
          Bs(is-i,j,k,2) = Bs(ie+1-i,j,k,2)
    enddo
    enddo
    enddo

    do k=ks,ke+1
    do j=js-mgn,je+mgn
    do i=1,mgn
          Bs(is-i,j,k,3) = Bs(ie+1-i,j,k,3)
    enddo
    enddo
    enddo

    ! x outer boundary
    do k=ks,ke
    do j=js-mgn,je+mgn
    do i=1,mgn
        Q(ie+i,j,k,:) = Q(is+i-1,j,k,:)
    enddo
    enddo
    enddo

    do k=ks,ke
    do j=js-mgn,je+mgn
    do i=1,mgn
        Bs(ie+i+1,j,k,1) = Bs(is+i,j,k,1)
    enddo
    enddo
    enddo

    do k=ks,ke
    do j=js-mgn,je+mgn+1
    do i=1,mgn
        Bs(ie+i,j,k,2) = Bs(is+i-1,j,k,2)
    enddo
    enddo
    enddo

    do k=ks,ke+1
    do j=js-mgn,je+mgn
    do i=1,mgn
        Bs(ie+i,j,k,3) = Bs(is+i-1,j,k,3)
    enddo
    enddo
    enddo

    ! y inner boundary
    do k=ks,ke
    do j=1,mgn
    do i=is-mgn,ie+mgn
        Q(i,js-j,k,:)  = Q(i,je+1-j,k,:)
    enddo
    enddo
    enddo

    do k=ks,ke
    do j=1,mgn
    do i=is-mgn,ie+mgn+1
        Bs(i,js-j,k,1) = Bs(i,je+1-j,k,1)
    enddo
    enddo
    enddo

    do k=ks,ke
    do j=1,mgn
    do i=is-mgn,ie+mgn
        Bs(i,js-j,k,2) = Bs(i,je+1-j,k,2)
    enddo
    enddo
    enddo


    do k=ks,ke+1
    do j=1,mgn
    do i=is-mgn,ie+mgn
        Bs(i,js-j,k,3) = Bs(i,je+1-j,k,3)
    enddo
    enddo
    enddo

    ! y outer boundary
    do k=ks,ke
    do j=1,mgn
    do i=is-mgn,ie+mgn
       Q(i,je+j,k,:)  = Q(i,js+j-1,k,:)
    enddo
    enddo
    enddo

    do k=ks,ke
    do j=1,mgn
    do i=is-mgn,ie+mgn+1
          Bs(i,je+j,k,1) = Bs(i,js+j-1,k,1)
    enddo
    enddo
    enddo

    do k=ks,ke
    do j=1,mgn
    do i=is-mgn,ie+mgn
        Bs(i,je+j+1,k,2) = Bs(i,js+j,k,2)
    enddo
    enddo
    enddo

    do k=ks,ke+1
    do j=1,mgn
    do i=is-mgn,ie+mgn
        Bs(i,je+j,k,3) = Bs(i,js+j-1,k,3)
    enddo
    enddo
    enddo


    ! boundary condition for the cell centered B field
    call CellCenterMagneticField(is-mgn, is-1, js-mgn, je+mgn, ks, ke, Bs, Bc)
    call CellCenterMagneticField(ie+1,   ie+mgn, js-mgn, je+mgn, ks, ke, Bs, Bc)
    call CellCenterMagneticField(is-mgn, ie+mgn, js-mgn, js-1, ks, ke, Bs, Bc)
    call CellCenterMagneticField(is-mgn, ie+mgn, je+1, je+mgn, ks, ke, Bs, Bc)
    return
    end subroutine BoundaryCondition
!
    subroutine ConsvVariable(Q, Bc, U)
    implicit none
    real(8), intent(in) :: Q(:,:,:,:)
    real(8), intent(in) :: Bc(:,:,:,:)
    real(8), intent(out) :: U(:,:,:,:)
    integer::i,j,k

        do k=ks,ke
        do j=js,je
        do i=is,ie
            U(i,j,k,IDN) = Q(i,j,k,IDN)
            U(i,j,k,IM1) = Q(i,j,k,IDN)*Q(i,j,k,IVX)
            U(i,j,k,IM2) = Q(i,j,k,IDN)*Q(i,j,k,IVY)
            U(i,j,k,IM3) = Q(i,j,k,IDN)*Q(i,j,k,IVZ)
            U(i,j,k,IEN) = 0.5d0*Q(i,j,k,IDN)*( Q(i,j,k,IVX)**2 + Q(i,j,k,IVY)**2 + Q(i,j,k,IVZ)**2 ) &
                       + 0.5d0*( Bc(i,j,k,1)**2 + Bc(i,j,k,2)**2 + Bc(i,j,k,3)**2 ) &
                       + Q(i,j,k,IPR)/(gam - 1.0d0)
            U(i,j,k,ISC) = Q(i,j,k,IDN)*Q(i,j,k,ISC)
        enddo
        enddo
        enddo
      
    return
    end subroutine Consvvariable

    subroutine PrimVariable( U, Bs, Q, Bc )
    implicit none
    real(8), intent(in) :: U(:,:,:,:),Bs(:,:,:,:)
    real(8), intent(out) :: Q(:,:,:,:),Bc(:,:,:,:)
    integer::i,j,k
    real(8) :: inv_d;

        call CellCenterMagneticField(is, ie, js, je, ks, ke, Bs, Bc)

        do k=ks,ke
        do j=js,je
        do i=is,ie
            Q(i,j,k,IDN) = U(i,j,k,IDN)
            inv_d = 1.0d0/U(i,j,k,IDN)
            Q(i,j,k,IVX) = U(i,j,k,IM1)*inv_d
            Q(i,j,k,IVY) = U(i,j,k,IM2)*inv_d
            Q(i,j,k,IVZ) = U(i,j,k,IM3)*inv_d
            Q(i,j,k,IPR) = ( U(i,j,k,IEN) &
                        - 0.5d0*(U(i,j,k,IM1)**2 + U(i,j,k,IM2)**2 + U(i,j,k,IM3)**2)*inv_d  &
                        - 0.5d0*(Bc(i,j,k,1)**2 + Bc(i,j,k,2)**2 + Bc(i,j,k,3)**2) )*(gam-1.0d0)
            Q(i,j,k,ISC) = U(i,j,k,ISC)*inv_d
        enddo
        enddo
        enddo

    return
    end subroutine PrimVariable

    subroutine CellCenterMagneticField(ibeg, ifin, jbeg, jfin, kbeg, kfin, Bs, Bc )
    implicit none
    integer, intent(in) :: ibeg, ifin, jbeg, jfin, kbeg, kfin
    real(8), intent(in) :: Bs(:,:,:,:)
    real(8), intent(out) :: Bc(:,:,:,:)
    integer::i,j,k
    real(8) :: inv_d;

        do k=kbeg,kfin
        do j=jbeg,jfin
        do i=ibeg,ifin
            Bc(i,j,k,1) = 0.5d0*( Bs(i+1,j,k,1) + Bs(i,j,k,1) )
            Bc(i,j,k,2) = 0.5d0*( Bs(i,j+1,k,2) + Bs(i,j,k,2) )
            Bc(i,j,k,3) = 0.5d0*( Bs(i,j,k+1,3) + Bs(i,j,k,3) )
        enddo
        enddo
        enddo

    return
    end subroutine CellCenterMagneticField

    subroutine TimestepControl(xf, yf, zf, Q, Bc, dt1)
    implicit none
    real(8), intent(in) :: xf(:), yf(:), zf(:), Q(:,:,:,:), Bc(:,:,:,:)
    real(8), intent(out) :: dt1
    real(8)::dtl1,dtl2,dtl3
    real(8)::dtlocal,dtmin,cf
    integer::i,j,k

        dtmin=1.0d90

        do k=ks,ke
        do j=js,je
        do i=is,ie
            cf = dsqrt( (gam*Q(i,j,k,IPR) + Bc(i,j,k,1)**2 + Bc(i,j,k,2)**2 + Bc(i,j,k,3)**2)/Q(i,j,k,IDN))
         
            dtl1 =(xf(i+1)-xf(i))/(abs(Q(i,j,k,IVX)) + cf)
            dtl2 =(yf(j+1)-yf(j))/(abs(Q(i,j,k,IVY)) + cf)
            dtlocal = min(dtl1,dtl2)
            if(dtlocal .lt. dtmin) dtmin = dtlocal
        enddo
        enddo
        enddo

        dt1 = Ccfl* dtmin*2.0d0

    return
    end subroutine TimestepControl

!---------------------------------------------------------------------
!     van Leer monotonicity limiter 
!---------------------------------------------------------------------
    subroutine vanLeer(n,dvp,dvm,dv)
    implicit none
    real(8),intent(in)::dvp(:),dvm(:)
    integer,intent(in) :: n
    real(8),intent(out)::dv(:)
    real(8) :: sgn
    integer :: i

        do i=1,n
            if(dvp(i)*dvm(i) .gt. 0.0d0) then
                dv(i) = 2.0d0*dvp(i)*dvm(i)/(dvp(i)+dvm(i))
            else
                dv(i) = 0.0d0
            endif
        enddo

    return
    end subroutine vanLeer

!---------------------------------------------------------------------
!     NumericalFlux
!---------------------------------------------------------------------
!     computes the numerical flux at the cell boundary 
!
!     Input: Q: primitive variables at the cell center
!
!     Input: B: magnetic fields
!
!     Output: flux : the numerical flux estimated at the cell boundary
!---------------------------------------------------------------------
    subroutine NumericalFlux( xf, yf, zf, Q, Bc, F, G, H, E)
    implicit none
    real(8), intent(in) :: xf(:), yf(:), zf(:)
    real(8), intent(in) :: Q(:,:,:,:)
    real(8), intent(in) :: Bc(:,:,:,:)
    real(8), intent(out) :: F(:,:,:,:)
    real(8), intent(out) :: G(:,:,:,:)
    real(8), intent(out) :: H(:,:,:,:)
    real(8), intent(out) :: E(:,:,:,:)
    
    integer::i,j,k
    real(8),dimension(nxtot,nytot,nztot,NFLX):: Ql,Qr
    real(8),dimension(NFLX):: flx
    real(8) :: dQm(NVAR), dQp(NVAR), dQmon(NVAR)
    real(8) :: ddmon, dvmon, dpmon
    real(8) :: Qltest(NFLX), Qrtest(NFLX);
    
    real(8),dimension(nxtot,nytot,nztot) :: e2_xf, e3_xf
    real(8),dimension(nxtot,nytot,nztot) :: e1_yf, e3_yf
    real(8),dimension(nxtot,nytot,nztot) :: weight1, weight2, weight3
    real(8) :: wghtCT
    
    
    ! numerical flux in the x direction
    ! hydro part
        do k=ks,ke
        do j=js-1,je+1
        do i=is-1,ie+1
            dQp(1:NVAR) = Q(i+1,j,k,1:NVAR) - Q(i  ,j,k,1:NVAR)
            dQm(1:NVAR) = Q(i  ,j,k,1:NVAR) - Q(i-1,j,k,1:NVAR)
    
            call vanLeer(NVAR, dQp, dQm, dQmon)
    
             ! Ql(i,j,k) --> W_(i-1/2,j,k)
             ! Qr(i,j,k) --> W_(i-1/2,j,k)
            Ql(i+1,j,k,1:NVAR) = Q(i,j,k,1:NVAR) + 0.5d0*dQmon(1:NVAR)
            Qr(i  ,j,k,1:NVAR) = Q(i,j,k,1:NVAR) - 0.5d0*dQmon(1:NVAR)
        enddo
        enddo
        enddo
    
        ! B field part
        do k=ks,ke
        do j=js-1,je+1
        do i=is-1,ie+1
            dQp(1:3) = Bc(i+1,j,k,1:3) - Bc(i  ,j,k,1:3)
            dQm(1:3) = Bc(i  ,j,k,1:3) - Bc(i-1,j,k,1:3)
    
            call vanLeer(3, dQp, dQm, dQmon)
    
             ! Ql(i,j,k) --> W_(i-1/2,j,k)
             ! Qr(i,j,k) --> W_(i-1/2,j,k)
            Ql(i+1,j,k,NVAR+1:NFLX) = Bc(i,j,k,1:3) + 0.5d0*dQmon(1:3)
            Qr(i  ,j,k,NVAR+1:NFLX) = Bc(i,j,k,1:3) - 0.5d0*dQmon(1:3)
        enddo
        enddo
        enddo
    
        if (flag_flux == 1 ) then
            do k=ks,ke
            do j=js-1,je+1
            do i=is,ie+1
                call HLL(1,Ql(i,j,k,:),Qr(i,j,k,:),Bs(i,j,k,1),xf(i+1)-xf(i),flx,wghtCT)
                 if (flx(IDN) >= 0.0) then 
                     flx(ISC) = flx(IDN)*Ql(i,j,k,ISC);
                 else 
                     flx(ISC) = flx(IDN)*Qr(i,j,k,ISC);
                 endif
        
                 F(i,j,k,1:NVAR)  = flx(1:NVAR)
                 e3_xf(i,j,k) =  -flx(IB2)
                 e2_xf(i,j,k) =  +flx(IB3)
        
                 if (flx(IDN) >= 0.0) then 
                     flx(ISC) = flx(IDN)*Ql(i,j,k,ISC);
                 else 
                     flx(ISC) = flx(IDN)*Qr(i,j,k,ISC);
                 endif

!                 if(flx(IDN).ne.flx(IDN)) then
!                     print*,i,j,k,flx(IDN),Ql(i,j,k,IDN),Qr(i,j,k,IDN)
!                 endif
        
                 weight1(i,j,k) = wghtCT
            enddo
            enddo
            enddo
        else if (flag_flux == 3 ) then
            do k=ks,ke
            do j=js-1,je+1
            do i=is,ie+1
                call HLLD(1,Ql(i,j,k,:),Qr(i,j,k,:),Bs(i,j,k,1),xf(i+1)-xf(i),flx,wghtCT)
                 if (flx(IDN) >= 0.0) then 
                     flx(ISC) = flx(IDN)*Ql(i,j,k,ISC);
                 else 
                     flx(ISC) = flx(IDN)*Qr(i,j,k,ISC);
                 endif
        
                 F(i,j,k,1:NVAR)  = flx(1:NVAR)
                 e3_xf(i,j,k) =  -flx(IB2)
                 e2_xf(i,j,k) =  +flx(IB3)

!                 if(flx(IDN).ne.flx(IDN)) then
!                     print*,i,j,k,flx(IDN),Ql(i,j,k,IDN),Qr(i,j,k,IDN)
!                 endif
        
        
                 weight1(i,j,k) = wghtCT
            enddo
            enddo
            enddo
        end if 
    
    
          ! numerical flux in the y direction
          do k=ks,ke
          do j=js-1,je+1
          do i=is-1,ie+1
             dQp(1:NVAR) = Q(i,j+1,k,1:NVAR) - Q(i,j  ,k,1:NVAR)
             dQm(1:NVAR) = Q(i,j  ,k,1:NVAR) - Q(i,j-1,k,1:NVAR)
    
             call vanLeer(NVAR, dQp, dQm, dQmon)
    
             ! Ql(i,j,k) --> W_(i-1/2,j,k)
             ! Qr(i,j,k) --> W_(i-1/2,j,k)
             Ql(i,j+1,k,1:NVAR) = Q(i,j,k,1:NVAR) + 0.5d0*dQmon(1:NVAR)
             Qr(i,j  ,k,1:NVAR) = Q(i,j,k,1:NVAR) - 0.5d0*dQmon(1:NVAR)
          enddo
          enddo
          enddo
    
          ! B field part
          do k=ks,ke
          do j=js-1,je+1
          do i=is-1,ie+1
             dQp(1:3) = Bc(i,j+1,k,1:3) - Bc(i,j  ,k,1:3)
             dQm(1:3) = Bc(i,j  ,k,1:3) - Bc(i,j-1,k,1:3)
    
             call vanLeer(3, dQp, dQm, dQmon)
    
             ! Ql(i,j,k) --> W_(i-1/2,j,k)
             ! Qr(i,j,k) --> W_(i-1/2,j,k)
             Ql(i,j+1,k,NVAR+1:NFLX) = Bc(i,j,k,1:3) + 0.5d0*dQmon(1:3)
             Qr(i,j  ,k,NVAR+1:NFLX) = Bc(i,j,k,1:3) - 0.5d0*dQmon(1:3)
          enddo
          enddo
          enddo
    
        if( flag_flux == 1 ) then
          do k=ks,ke
          do j=js,je+1
          do i=is-1,ie+1
             call HLL(2,Ql(i,j,k,:),Qr(i,j,k,:),Bs(i,j,k,2),yf(j+1) - yf(j), flx,wghtCT)
             if (flx(IDN) >= 0.0) then 
                 flx(ISC) = flx(IDN)*Ql(i,j,k,ISC);
             else 
                 flx(ISC) = flx(IDN)*Qr(i,j,k,ISC);
             endif
    
             G(i,j,k,1:NVAR) = flx(1:NVAR)
             e1_yf(i,j,k) =  - flx(IB3)
             e3_yf(i,j,k) =    flx(IB1)
    
    
             weight2(i,j,k) = wghtCT
          enddo
          enddo
          enddo
        else if (flag_flux == 3 ) then
          do k=ks,ke
          do j=js,je+1
          do i=is-1,ie+1
             call HLLD(2,Ql(i,j,k,:),Qr(i,j,k,:),Bs(i,j,k,2),yf(j+1) - yf(j), flx,wghtCT)

             if (flx(IDN) >= 0.0) then 
                 flx(ISC) = flx(IDN)*Ql(i,j,k,ISC);
             else 
                 flx(ISC) = flx(IDN)*Qr(i,j,k,ISC);
             endif
    
             G(i,j,k,1:NVAR) = flx(1:NVAR)
             e1_yf(i,j,k) =  - flx(IB3)
             e3_yf(i,j,k) =    flx(IB1)
    
    
             weight2(i,j,k) = wghtCT
          enddo
          enddo
          enddo
        endif
    
          call ElectricField( Q, Bc, e2_xf, e3_xf, e3_yf, e1_yf, weight1, weight2, E )
    
        return
        end subroutine Numericalflux
    
!---------------------------------------------------------------------
!     HLL Riemann Solver
!---------------------------------------------------------------------
!     solve the HLL Riemann solver 
!
!     Input: Ql, Qr: primitive variables containing the perpendicular B fields 
!                    at the left and right states
!            1D array (IDN, IVX, IVY, IVZ, IPR, IBperp1, IBperp2)
!                                 |
!                                 |
!                           Ql    |    Qr
!                                 |
!                                -->
!                                flx
!
!     Input: b1    : magnetic field perpendicular to the initial discontinuity
!
!     Output: flx  : flux estimated at the initial discontinuity
!            index: (IDN, IVX, IVY, IVZ, IPR, IBperp1, IBperp2)
!---------------------------------------------------------------------
      subroutine HLL(idir,Ql,Qr,b1,dx,flx,wghtCT)
      implicit none
      integer, intent(in) :: idir
      real(8),intent(in)  :: Ql(:), Qr(:)
      real(8),intent(in)  :: dx, b1
      real(8),intent(out) :: flx(:), wghtCT
      integer :: IVpara, IVperp1, IVperp2
      integer :: IBpara, IBperp1, IBperp2
      real(8):: Ul(NFLX), Ur(NFLX)
      real(8):: Fl(NFLX), Fr(NFLX)
      real(8):: Ust(NFLX)
      real(8):: Fst(NFLX)
      real(8):: cfl,cfr
      real(8):: sl, sr
      real(8):: pbl, pbr, ptotl, ptotr, v_over_c
      integer :: i, n

      if( idir == 1 ) then
           IVpara  = IVX
           IVperp1 = IVY
           IVperp2 = IVZ
           IBpara  = IB1
           IBperp1 = IB2
           IBperp2 = IB3
      else if (idir == 2 ) then
           IVpara  = IVY
           IVperp1 = IVZ
           IVperp2 = IVX
           IBpara  = IB2
           IBperp1 = IB3
           IBperp2 = IB1
      endif
          
          pbl = 0.5d0*(b1**2 + Ql(IBperp1)**2 + Ql(IBperp2)**2)
          pbr = 0.5d0*(b1**2 + Qr(IBperp1)**2 + Qr(IBperp2)**2)
          ptotl = Ql(IPR) + pbl
          ptotr = Qr(IPR) + pbr

          ! conserved variables in the left and right states
          Ul(IDN) = Ql(IDN)
          Ul(IVpara) = Ql(IDN)*Ql(IVpara)
          Ul(IVperp1) = Ql(IDN)*Ql(IVperp1)
          Ul(IVperp2) = Ql(IDN)*Ql(IVperp2)
          Ul(IEN) = 0.5d0*Ql(IDN)*( Ql(IVpara)**2 + Ql(IVperp1)**2 + Ql(IVperp2)**2) & 
                  + pbl + Ql(IPR)/(gam - 1.0d0)
          Ul(IBperp1) = Ql(IBperp1)
          Ul(IBperp2) = Ql(IBperp2)

          Ur(IDN) = Qr(IDN)
          Ur(IVpara) = Qr(IDN)*Qr(IVpara)
          Ur(IVperp1) = Qr(IDN)*Qr(IVperp1)
          Ur(IVperp2) = Qr(IDN)*Qr(IVperp2)
          Ur(IEN) = 0.5d0*Qr(IDN)*( Qr(IVpara)**2 + Qr(IVperp1)**2 + Qr(IVperp2)**2) & 
                  + pbr + Qr(IPR)/(gam - 1.0d0)
          Ur(IBperp1) = Qr(IBperp1)
          Ur(IBperp2) = Qr(IBperp2)

    !--- Step 3.  Compute L/R fluxes
          Fl(IDN) = Ul(IVpara)
          Fl(IVpara) = Ul(IVpara)*Ql(IVpara) + ptotl - b1**2
          Fl(IVperp1) = Ul(IVperp1)*Ql(IVpara) - b1*Ql(IBperp1)
          Fl(IVperp2) = Ul(IVperp2)*Ql(IVpara) - b1*Ql(IBperp2)
          Fl(IEN) = ( Ul(IEN) + ptotl - b1**2 )*Ql(IVpara) &
                  - b1*( Ql(IBperp1)*Ql(IVperp1) + Ql(IBperp2)*Ql(IVperp2) )
          Fl(IBperp1) = Ql(IBperp1)*Ql(IVpara) - b1*Ql(IVperp1)
          Fl(IBperp2) = Ql(IBperp2)*Ql(IVpara) - b1*Ql(IVperp2)

          Fr(IDN) = Ur(IVpara)
          Fr(IVpara) = Ur(IVpara)*Qr(IVpara) + ptotr - b1**2
          Fr(IVperp1) = Ur(IVperp1)*Qr(IVpara) - b1*Qr(IBperp1)
          Fr(IVperp2) = Ur(IVperp2)*Qr(IVpara) - b1*Qr(IBperp2)
          Fr(IEN) = ( Ur(IEN) + ptotr - b1**2 )*Qr(IVpara) &
                  - b1*( Qr(IBperp1)*Qr(IVperp1) + Qr(IBperp2)*Qr(IVperp2) )
          Fr(IBperp1) = Qr(IBperp1)*Qr(IVpara) - b1*Qr(IVperp1)
          Fr(IBperp2) = Qr(IBperp2)*Qr(IVpara) - b1*Qr(IVperp2)

!          cfl = dsqrt( (gam*Ql(IPR) + Ql(IBperp1)**2 + Ql(IBperp2)**2 + b1**2)/Ql(IDN))
!          cfr = dsqrt( (gam*Qr(IPR) + Qr(IBperp1)**2 + Qr(IBperp2)**2 + b1**2)/Qr(IDN))
          cfl = dsqrt( 0.5d0*( 2.0d0*pbl + gam*Ql(IPR) &
                     + dsqrt( (2.0d0*pbl - gam*Ql(IPR))**2 &
                     + 4.0d0*gam*Ql(IPR)*( Ql(IBperp1)**2 + Ql(IBperp2)**2 ) ) )/Ql(IDN) )
          cfr = dsqrt( 0.5d0*( 2.0d0*pbr + gam*Qr(IPR) &
                    + dsqrt( (2.0d0*pbr - gam*Qr(IPR))**2 &
                     + 4.0d0*gam*Qr(IPR)*( Qr(IBperp1)**2 + Qr(IBperp2)**2 ) ) )/Qr(IDN) )

!          sl = min(Ql(IVX),Qr(IVX)) - max(cfl,cfr)
!          sr = max(Ql(IVX),Qr(IVX)) + max(cfl,cfr)
          sl = min(Ql(IVpara) - cfl,Qr(IVpara) - cfr)
          sr = max(Ql(IVpara) + cfl,Qr(IVpara) + cfr)
!          Fst(:)  = (sr*Fl(:) - sl*Fr(:) + sl*sr*( Ur(:) - Ul(:) ))/(sr - sl)
!          Ust(:) = ( sr*Ur(:) - sl*Ul(:) - Fr(:) + Fl(:) )/(sr - sl)

          if( sl > 0.0d0 ) then
               flx(:) = Fl(:)
          else if (sr <= 0.0d0 ) then
               flx(:) = Fr(:)
          else 
               flx(:)  = (sr*Fl(:) - sl*Fr(:) + sl*sr*( Ur(:) - Ul(:) ))/(sr - sl)
          endif

           v_over_c = 1024.0d0*dt*flx(IDN)/( dx*( Ql(IDN) + Qr(IDN) ) )
           wghtCT = 0.5d0 + max( -0.5d0, min(0.5d0, v_over_c) )

!          flx(IPS) = 0.0d0

!         do i=1,NFLX1D 
!         if( flx(i) .ne. flx(i) ) then 
!             print*,(sr*Fl(:) - sl*Fr(:) + sl*sr*( Ur(:) - Ul(:) ))/(sr - sl)
!         stop
!         endif
!         enddo
!
      return
      end subroutine HLL

!---------------------------------------------------------------------
!     HLLD Riemann Solver
!---------------------------------------------------------------------
!     solve the HLL Riemann solver 
!
!     Input: Ql, Qr: primitive variables containing the perpendicular B fields 
!                    at the left and right states
!            1D array (IDN, IVX, IVY, IVZ, IPR, IBperp1, IBperp2)
!                                 |
!                                 |
!                           Ql    |    Qr
!                                 |
!                                -->
!                                flx
!
!     Input: b1    : magnetic field perpendicular to the initial discontinuity
!
!     Output: flx  : flux estimated at the initial discontinuity
!            index: (IDN, IVX, IVY, IVZ, IPR, IBperp1, IBperp2)
!---------------------------------------------------------------------
      subroutine HLLD(idir,Ql,Qr,b1,dx,flx,wghtCT)
      implicit none
      integer, intent(in) :: idir
      real(8),intent(in)  :: Ql(:), Qr(:)
      real(8),intent(in)  :: dx
      real(8),intent(in) :: b1
      real(8),intent(out) :: flx(:)
      real(8),intent(out) :: wghtCT
      integer :: IVpara, IVperp1, IVperp2
      integer :: IBpara, IBperp1, IBperp2, id
      real(8):: Ul(NFLX), Ur(NFLX)
      real(8):: Ulst(NFLX), Urst(NFLX)
      real(8):: Uldst(NFLX), Urdst(NFLX)
      real(8):: Fl(NFLX), Fr(NFLX)
      real(8):: test(NFLX)
      real(8):: cfl,cfr
      real(8):: S0, S1, S2, S3, S4
      real(8):: pbl, pbr, ptotl, ptotr
      real(8) :: sqrtdl, sqrtdr, v_dot_B_stl, v_dot_B_str
      real(8) :: Ulst_d_inv, Urst_d_inv, sum_sqrtd_inv, tmp
      real(8) :: ptot_stl, ptot_str,ptot_st, Cl, Cr, Cml, Cmr, Cml_inv, Cmr_inv, bxsgn
      real(8) :: v_over_c
      integer :: i, n

      if( idir == 1 ) then
           IVpara  = IVX
           IVperp1 = IVY
           IVperp2 = IVZ
           IBpara  = IB1
           IBperp1 = IB2
           IBperp2 = IB3
      else if (idir == 2 ) then
           IVpara  = IVY
           IVperp1 = IVZ
           IVperp2 = IVX
           IBpara  = IB2
           IBperp1 = IB3
           IBperp2 = IB1
      endif
          
          pbl = 0.5d0*(b1**2 + Ql(IBperp1)**2 + Ql(IBperp2)**2)
          pbr = 0.5d0*(b1**2 + Qr(IBperp1)**2 + Qr(IBperp2)**2)
          ptotl = Ql(IPR) + pbl
          ptotr = Qr(IPR) + pbr

          cfl = dsqrt( 0.5d0*( 2.0d0*pbl + gam*Ql(IPR) &
                     + dsqrt( (2.0d0*pbl - gam*Ql(IPR))**2 &
                     + 4.0d0*gam*Ql(IPR)*( Ql(IBperp1)**2 + Ql(IBperp2)**2 ) ) )/Ql(IDN) )
          cfr = dsqrt( 0.5d0*( 2.0d0*pbr + gam*Qr(IPR) &
                    + dsqrt( (2.0d0*pbr - gam*Qr(IPR))**2 &
                     + 4.0d0*gam*Qr(IPR)*( Qr(IBperp1)**2 + Qr(IBperp2)**2 ) ) )/Qr(IDN) )
!          cfl = dsqrt( (gam*Ql(IPR) + Ql(IBperp1)**2 + Ql(IBperp2)**2 + b1**2)/Ql(IDN))
!          cfr = dsqrt( (gam*Qr(IPR) + Qr(IBperp1)**2 + Qr(IBperp2)**2 + b1**2)/Qr(IDN))

          S0 = min( Ql(IVpara) - cfl, Qr(IVpara) - cfr)
          S4 = max( Ql(IVpara) + cfl, Qr(IVpara) + cfr)

          ! conserved variables in the left and right states
          Ul(IDN) = Ql(IDN)
          Ul(IVpara) = Ql(IDN)*Ql(IVpara)
          Ul(IVperp1) = Ql(IDN)*Ql(IVperp1)
          Ul(IVperp2) = Ql(IDN)*Ql(IVperp2)
          Ul(IEN) = 0.5d0*Ql(IDN)*( Ql(IVpara)**2 + Ql(IVperp1)**2 + Ql(IVperp2)**2) & 
                  + pbl + Ql(IPR)/(gam - 1.0d0)
          Ul(IBperp1) = Ql(IBperp1)
          Ul(IBperp2) = Ql(IBperp2)

          Ur(IDN) = Qr(IDN)
          Ur(IVpara) = Qr(IDN)*Qr(IVpara)
          Ur(IVperp1) = Qr(IDN)*Qr(IVperp1)
          Ur(IVperp2) = Qr(IDN)*Qr(IVperp2)
          Ur(IEN) = 0.5d0*Qr(IDN)*( Qr(IVpara)**2 + Qr(IVperp1)**2 + Qr(IVperp2)**2) & 
                  + pbr + Qr(IPR)/(gam - 1.0d0)
          Ur(IBperp1) = Qr(IBperp1)
          Ur(IBperp2) = Qr(IBperp2)

    !--- Step 3.  Compute L/R fluxes
          Fl(IDN) = Ul(IVpara)
          Fl(IVpara) = Ul(IVpara)*Ql(IVpara) + ptotl - b1**2
          Fl(IVperp1) = Ul(IVperp1)*Ql(IVpara) - b1*Ql(IBperp1)
          Fl(IVperp2) = Ul(IVperp2)*Ql(IVpara) - b1*Ql(IBperp2)
          Fl(IEN) = ( Ul(IEN) + ptotl - b1**2 )*Ql(IVpara) &
                  - b1*( Ql(IBperp1)*Ql(IVperp1) + Ql(IBperp2)*Ql(IVperp2) )
          Fl(IBperp1) = Ql(IBperp1)*Ql(IVpara) - b1*Ql(IVperp1)
          Fl(IBperp2) = Ql(IBperp2)*Ql(IVpara) - b1*Ql(IVperp2)

          Fr(IDN) = Ur(IVpara)
          Fr(IVpara) = Ur(IVpara)*Qr(IVpara) + ptotr - b1**2
          Fr(IVperp1) = Ur(IVperp1)*Qr(IVpara) - b1*Qr(IBperp1)
          Fr(IVperp2) = Ur(IVperp2)*Qr(IVpara) - b1*Qr(IBperp2)
          Fr(IEN) = ( Ur(IEN) + ptotr - b1**2 )*Qr(IVpara) &
                  - b1*( Qr(IBperp1)*Qr(IVperp1) + Qr(IBperp2)*Qr(IVperp2) )
          Fr(IBperp1) = Qr(IBperp1)*Qr(IVpara) - b1*Qr(IVperp1)
          Fr(IBperp2) = Qr(IBperp2)*Qr(IVpara) - b1*Qr(IVperp2)

    !--- Step 4.  Compute middle and Alfven wave speeds
          Cl = S0 - Ql(IVpara)
          Cr = S4 - Qr(IVpara)

          S2 = ( Cr*Ur(IVpara) - Cl*Ul(IVpara) + (ptotl - ptotr) ) &
                 /( Cr*Ur(IDN) - Cl*Ul(IDN) )

          Cml = S0 - S2
          Cmr = S4 - S2
          Cml_inv = 1.0d0/Cml
          Cmr_inv = 1.0d0/Cmr

          Ulst(IDN) = Ul(IDN)*Cl*Cml_inv
          Urst(IDN) = Ur(IDN)*Cr*Cmr_inv
          Ulst_d_inv = 1.0d0/Ulst(IDN)
          Urst_d_inv = 1.0d0/Urst(IDN)
          sqrtdl = dsqrt(Ulst(IDN))
          sqrtdr = dsqrt(Urst(IDN))


          S1 = S2 - dabs(b1)/sqrtdl
          S3 = S2 + dabs(b1)/sqrtdr

    !--- Step 5.  Compute intermediate states
         ptot_stl = ptotl + Ul(IDN)*Cl*(S2 - Ql(IVpara))
         ptot_str = ptotr + Ur(IDN)*Cr*(S2 - Qr(IVpara))

         ptot_st = 0.5d0*(ptot_stl + ptot_str)

         Ulst(IVpara) = Ulst(IDN)*S2
         if( dabs( Ul(IDN)*Cl*Cml-b1**2) < 1.0d-8*ptot_st ) then
             Ulst(IVperp1) = Ulst(IDN)*Ql(IVperp1)
             Ulst(IVperp2) = Ulst(IDN)*Ql(IVperp2)

             Ulst(IBperp1) = Ul(IBperp1)
             Ulst(IBperp2) = Ul(IBperp2)
         else 
             tmp = b1*( Cl - Cml )/(Ul(IDN)*Cl*Cml - b1**2)
             Ulst(IVperp1) = Ulst(IDN)*( Ql(IVperp1) - Ul(IBperp1)*tmp )
             Ulst(IVperp2) = Ulst(IDN)*( Ql(IVperp2) - Ul(IBperp2)*tmp )

             tmp = (Ul(IDN)*Cl**2 - b1**2)/( Ul(IDN)*Cl*Cml - b1**2)
             Ulst(IBperp1) = Ul(IBperp1)*tmp
             Ulst(IBperp2) = Ul(IBperp2)*tmp
         endif

         v_dot_B_stl = ( Ulst(IVpara)*b1 + Ulst(IVperp1)*Ulst(IBperp1) + Ulst(IVperp2)*Ulst(IBperp2) )*Ulst_d_inv
         Ulst(IEN) = ( Cl*Ul(IEN) - ptotl*Ql(IVpara) + ptot_st*S2 &
                     + b1*( Ql(IVpara)*b1 + Ql(IVperp1)*Ul(IBperp1) + Ql(IVperp2)*Ul(IBperp2) - v_dot_B_stl) )*Cml_inv

         Urst(IVpara) = Urst(IDN)*S2
         if( dabs( Ur(IDN)*Cr*Cmr-b1**2) < 1.0d-8*ptot_st ) then
             Urst(IVperp1) = Urst(IDN)*Qr(IVperp1)
             Urst(IVperp2) = Urst(IDN)*Qr(IVperp2)

             Urst(IBperp1) = Ur(IBperp1)
             Urst(IBperp2) = Ur(IBperp2)
         else 
             tmp = b1*( Cr - Cmr )/(Ur(IDN)*Cr*Cmr - b1**2)
             Urst(IVperp1) = Urst(IDN)*( Qr(IVperp1) - Ur(IBperp1)*tmp )
             Urst(IVperp2) = Urst(IDN)*( Qr(IVperp2) - Ur(IBperp2)*tmp )

             tmp = (Ur(IDN)*Cr**2 - b1**2)/( Ur(IDN)*Cr*Cmr - b1**2)
             Urst(IBperp1) = Ur(IBperp1)*tmp
             Urst(IBperp2) = Ur(IBperp2)*tmp
         endif

         v_dot_B_str = ( Urst(IVpara)*b1 + Urst(IVperp1)*Urst(IBperp1) + Urst(IVperp2)*Urst(IBperp2) )*Urst_d_inv
         Urst(IEN) = ( Cr*Ur(IEN) - ptotr*Qr(IVpara) + ptot_st*S2 &
                     + b1*( Qr(IVpara)*b1 + Qr(IVperp1)*Ur(IBperp1) + Qr(IVperp2)*Ur(IBperp2) - v_dot_B_str) )*Cmr_inv
       
         if( 0.5d0*b1**2 < 1.0d-8*ptot_st )then
             Uldst(:) = Ulst(:)
             Urdst(:) = Urst(:)
         else 
             sum_sqrtd_inv = 1.0d0/(sqrtdl + sqrtdr) 
             if (b1 > 0.0d0 ) then 
                 bxsgn = 1.0d0
             else 
                 bxsgn = -1.0d0
             endif

             Uldst(IDN) = Ulst(IDN)
             Urdst(IDN) = Urst(IDN)

             Uldst(IVpara) = Ulst(IVpara)
             Urdst(IVpara) = Urst(IVpara)

             tmp = sum_sqrtd_inv*(  sqrtdl*(Ulst(IVperp1)*Ulst_d_inv) + sqrtdr*(Urst(IVperp1)*Urst_d_inv) &
                                  + bxsgn*(Urst(IBperp1) - Ulst(IBperp1)) )
             Uldst(IVperp1) = Uldst(IDN)*tmp
             Urdst(IVperp1) = Urdst(IDN)*tmp
!

             tmp = sum_sqrtd_inv*(  sqrtdl*(Ulst(IVperp2)*Ulst_d_inv) + sqrtdr*(Urst(IVperp2)*Urst_d_inv) &
                                  + bxsgn*(Urst(IBperp2) - Ulst(IBperp2)) )
             Uldst(IVperp2) = Uldst(IDN)*tmp
             Urdst(IVperp2) = Urdst(IDN)*tmp

             tmp = sum_sqrtd_inv*(  sqrtdl*Urst(IBperp1) + sqrtdr*Ulst(IBperp1) &
                       + bxsgn*sqrtdl*sqrtdr*( (Urst(IVperp1)*Urst_d_inv) - (Ulst(IVperp1)*Ulst_d_inv) ) )
             Uldst(IBperp1) = tmp
             Urdst(IBperp1) = tmp

             tmp = sum_sqrtd_inv*(  sqrtdl*Urst(IBperp2) + sqrtdr*Ulst(IBperp2) &
                       + bxsgn*sqrtdl*sqrtdr*( (Urst(IVperp2)*Urst_d_inv) - (Ulst(IVperp2)*Ulst_d_inv) ) )
             Uldst(IBperp2) = tmp
             Urdst(IBperp2) = tmp
!
             tmp = S2*b1 + (Uldst(IVperp1)*Uldst(IBperp1) + Uldst(IVperp2)*Uldst(IBperp2))/Uldst(IDN)
             Uldst(IEN) = Ulst(IEN) - sqrtdl*bxsgn*(v_dot_B_stl - tmp)
             Urdst(IEN) = Urst(IEN) + sqrtdr*bxsgn*(v_dot_B_str - tmp)
         endif

!         test = (S4 - S3)*Urst + (S3 - S2)*Urdst + (S2 - S1)*Uldst + (S1 - S0)*Ulst - S4*Ur + S0*Ul + Fr - Fl
!         do id=1,NFLX
!         if( abs(test(id)) > 1.0d-10 ) then
!         print*,test(id)
!         stop
!         endif
!         enddo
         

    !--- Step 6.  Compute flux
          if( S0 >= 0.0d0 ) then
               flx(:) = Fl(:)
          else if (S4 <= 0.0d0 ) then
               flx(:) = Fr(:)
          else  if  (S1 >= 0.0d0) then
               flx(:) = Fl(:) + S0*(Ulst(:) - Ul(:))
           else if (S2 >= 0.0d0) then
               flx(:) = Fl(:) + S0*(Ulst(:) - Ul(:)) + S1*(Uldst(:) - Ulst(:))
           else if (S3 > 0.0d0 ) then
               flx(:) = Fr(:) + S4*(Urst(:) - Ur(:)) + S3*(Urdst(:) - Urst(:))
           else
               flx(:) = Fr(:) + S4*(Urst(:) - Ur(:)) 
           endif
           flx(IBpara) = 0.0d0

           v_over_c = 1024.0d0*dt*flx(IDN)/( dx*( Ql(IDN) + Qr(IDN) ) )
           wghtCT = 0.5d0 + max( -0.5d0, min(0.5d0, v_over_c) )

      return
      end subroutine HLLD

!---------------------------------------------------------------------
!     ElectricField
!---------------------------------------------------------------------
!     computes the numerical flux at the cell boundary 
!
!     Input: Q: primitive variables at the cell center
!
!     Input: B: magnetic fields
!
!     Output: flux : the numerical flux estimated at the cell boundary
!---------------------------------------------------------------------
      subroutine ElectricField( Q, Bc, e2_xf, e3_xf, e3_yf, e1_yf, weight1, weight2, E )
      implicit none
      integer::i,j,k
      real(8), intent(in)  :: Q(:,:,:,:), Bc(:,:,:,:)
      real(8), intent(in)  :: e2_xf(:,:,:)
      real(8), intent(in)  :: e3_xf(:,:,:)
      real(8), intent(in)  :: e1_yf(:,:,:)
      real(8), intent(in)  :: e3_yf(:,:,:)
      real(8), intent(in)  :: weight1(:,:,:)
      real(8), intent(in)  :: weight2(:,:,:)
      real(8), intent(out) :: E(:,:,:,:)
      real(8) :: Etmp(nxtot,nytot,nztot) 
      real(8) :: de3_l1, de3_r1, de3_l2, de3_r2

      do k=ks,ke
      do j=js-1, je+1
      do i=is-1, ie+1
           Etmp(i,j,k) = Q(i,j,k,IVY)*Bc(i,j,k,1) - Q(i,j,k,IVX)*Bc(i,j,k,2)
      enddo
      enddo
      enddo

      do j=js, je
      do i=is, ie+1
           E(i,j,ke+1,2) = e2_xf(i,j,ks)
           E(i,j,ks  ,2) = e2_xf(i,j,ks)
      enddo
      enddo

      do j=js, je+1
      do i=is, ie
           E(i,j,ke+1,1) = e1_yf(i,j,ks)
           E(i,j,ks  ,1) = e1_yf(i,j,ks)
      enddo
      enddo


      do k=ks,ke
      do j=js, je+1
      do i=is, ie+1
          de3_l2 = (1.0d0-weight1(i,j-1,k))*(e3_yf(i  ,j,k) - Etmp(i  ,j-1,k)) + &
                   (      weight1(i,j-1,k))*(e3_yf(i-1,j,k) - Etmp(i-1,j-1,k))

          de3_r2 = (1.0d0-weight1(i,j  ,k))*(e3_yf(i  ,j,k) - Etmp(i  ,j  ,k)) + &
                   (      weight1(i,j  ,k))*(e3_yf(i-1,j,k) - Etmp(i-1,j  ,k))

          de3_l1 = (1.0d0-weight2(i-1,j,k))*(e3_xf(i,j  ,k) - Etmp(i-1,j  ,k)) + &
                   (      weight2(i-1,j,k))*(e3_xf(i,j-1,k) - Etmp(i-1,j-1,k))

          de3_r1 = (1.0d0-weight2(i  ,j,k))*(e3_xf(i,j  ,k) - Etmp(i  ,j  ,k)) + &
                   (      weight2(i  ,j,k))*(e3_xf(i,j-1,k) - Etmp(i  ,j-1,k))

         E(i,j,k,3) = 0.25d0*( de3_l1 + de3_r1 + de3_l2 + de3_r2 + &
                          e3_yf(i-1,j,k) + e3_yf(i,j,k) + e3_xf(i,j-1,k) + e3_xf(i,j,k))
      enddo
      enddo
      enddo


      return
      end subroutine ElectricField 
!!=====================================================================

      subroutine UpdateConsv( dt1, xf, yf, zf, F, G, H, E, Q, Uo, Bso, U, Bs)
      implicit none
      real(8), intent(in) :: dt1
      real(8), intent(in)  :: xf(:), yf(:), zf(:)
      real(8), intent(in)  :: F(:,:,:,:), G(:,:,:,:), H(:,:,:,:)
      real(8), intent(in)  :: Uo(:,:,:,:), Q(:,:,:,:), Bso(:,:,:,:)
      real(8), intent(out) :: U(:,:,:,:), Bs(:,:,:,:), E(:,:,:,:)

      real(8) :: divB 
      integer::i,n,j,k

      do n=1,NVAR
      do k=ks,ke
      do j=js,je
      do i=is,ie
         U(i,j,k,n) = Uo(i,j,k,n) + dt1*(- F(i+1,j,k,n) + F(i,j,k,n))/(xf(i+1)-xf(i)) &
                                  + dt1*(- G(i,j+1,k,n) + G(i,j,k,n))/(yf(j+1)-yf(j)) 
      enddo
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=js,je
      do i=is,ie+1
           Bs(i,j,k,1) = Bso(i,j,k,1) &
                       - dt1*(E(i,j+1,k,3) - E(i,j,k,3))/(yv(j+1) - yv(j))
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=js,je+1
      do i=is,ie
           Bs(i,j,k,2) = Bso(i,j,k,2) &
                       + dt1*(E(i+1,j,k,3) - E(i,j,k,3))/(xv(i+1) - xv(i))
      enddo
      enddo
      enddo

      do k=ks,ke+1
      do j=js,je
      do i=is,ie
           Bs(i,j,k,3) = Bso(i,j,k,3) &
                       - dt1*(E(i+1,j,k,2) - E(i,j,k,2))/(xf(i+1) - xf(i)) &
                       + dt1*(E(i,j+1,k,1) - E(i,j,k,1))/(yf(j+1) - yf(j))
      enddo
      enddo
      enddo



      return
      end subroutine UpdateConsv

      subroutine Output( flag_output, xf, xv, yf, yv, Q, Bc )
      implicit none
      logical, intent(in) :: flag_output ! false --> output per dtout, true --> force to output
      real(8), intent(in) :: xf(:), xv(:), yf(:), yv(:)
      real(8), intent(in) :: Q(:,:,:,:), Bc(:,:,:,:)
      integer::i,j,k
      character(20),parameter::dirname="snap_B0.8_ct_hll"
      character(20),parameter::base="kh"
      character(20),parameter::suffix=".dat"
      character(40)::filename
      real(8), save :: tout = - dtout
      integer::nout = 0
      integer,parameter:: unitout=17
      integer,parameter:: unitbin=13

      logical, save:: is_inited
      data is_inited /.false./

      if (.not. is_inited) then
         call makedirs(dirname)
         is_inited =.true.
      endif

      if( .not.flag_output) then
          if( time + 1.0d-14.lt. tout+dtout) return
      endif

      write(filename,'(i5.5)') nout
      filename = trim(dirname)//"/"//trim(base)//trim(filename)//trim(suffix)
      open(unitbin,file=filename,form='formatted',action="write")
      write(unitbin,*) "# time = ",time
      write(unitbin,*) "#nx, ny = ", nx, ny
      do k=ks,ke
      do j=js,je
      do i=is,ie
          write(unitbin,*) xv(i), yv(j), Q(i,j,k,IDN), Q(i,j,k,IVX), Q(i,j,k,IVY), Q(i,j,k,IVZ), Q(i,j,k,IPR), &
              Bc(i,j,k,1), Bc(i,j,k,2), Bc(i,j,k,3), Q(i,j,k,ISC)
      enddo
      enddo
      enddo

     close(unitbin)

     write(6,*) "output:  ",filename,time

      nout=nout+1
      tout=tout + dtout

      return
      end subroutine Output

      subroutine makedirs(outdir)
      implicit none
      character(len=*), intent(in) :: outdir
      character(len=256) command
      write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
!      write(*, *) trim(command)
      call system(command)
      end subroutine makedirs

      real(8) function errorBpara(Q)
      implicit none
      real(8), intent(in) :: Q(:,:,:,:)
      integer::i,j,k
      real(8) :: error

      do k=ks,ke
      do j=js,je
      do i=is,ie
         error = dabs( ( Q(i,j,k,IB1) + 2.0d0*Q(i,j,k,IB2) )/sqrt(5.0d0) - 5.0d0/dsqrt(4.0d0*dacos(-1.0d0)))

         errorBpara = errorBpara + error
      enddo
      enddo
      enddo
      errorBpara = errorBpara/dble((ie-is+1)*(je-js+1)*(ke-ks+1))
      
      return
      end function

      real(8) function divergenceB(xf, yf, Bc)
      implicit none
      real(8), intent(in) :: xf(:), yf(:), Bc(:,:,:,:)
      integer::i,j,k
      real(8) :: error

      do k=ks,ke
      do j=js,je
      do i=is,ie
         error = dabs( ( Bc(i+1,j,k,1) - Bc(i-1,j,k,1) )/(xf(i+1)-xf(i-1))  &
                    + ( Bc(i,j+1,k,2) - Bc(i,j-1,k,2) )/(yf(j+1)-yf(j-1)) )/ &
                    dsqrt( Bc(i,j,k,IB1)**2 + Bc(i,j,k,2)**2 )*min(xf(i+1) - xf(i),yf(j+1)-yf(j))

         divergenceB = divergenceB + error
      enddo
      enddo
      enddo
      divergenceB = divergenceB/dble((ie-is+1)*(je-js+1)*(ke-ks+1))
      
      return
      end function


      subroutine Analysis(xv,yv,Q,Bc,Bs,phys_evo)
      real(8), intent(in) :: xv(:), yv(:), Q(:,:,:,:), Bc(:,:,:,:), Bs(:,:,:,:)
      real(8), intent(out) :: phys_evo(:)
      integer::i,j,k
      real(8) :: dvy, er_divBc, er_divBs

      dvy = 0.0d0
      er_divBc = 0.0d0
      er_divBs = 0.0d0
      do k=ks,ke
      do j=js,je
      do i=is,ie
           dvy = dvy + Q(i,j,k,IVY)**2
           er_divBs = er_divBs + ( Bs(i+1,j,k,1) - Bs(i,j,k,1) + Bs(i,j+1,k,2) - Bs(i,j,k,2) )**2 &
                       /( Bc(i,j,k,1)**2 + Bc(i,j,k,2)**2 )
           er_divBc = er_divBc + 0.5d0*( Bc(i+1,j,k,1) - Bc(i-1,j,k,1) + Bc(i,j+1,k,2) - Bc(i,j-1,k,2) )**2 &
                                       /( Bc(i,j,k,1)**2 + Bc(i,j,k,2)**2 )
      enddo
      enddo
      enddo
      phys_evo(1) = sqrt(dvy/dble(nx*ny))
      phys_evo(2) = sqrt(er_divBc/dble(nx*ny))
      phys_evo(3) = sqrt(er_divBs/dble(nx*ny))
      
      return
      end subroutine

end program main
