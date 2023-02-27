program main
!$ use omp_lib
implicit none

! time evolution
integer :: ntime = 0    ! counter of the timestep
real(8) :: time = 0.0d0  ! time 
real(8) :: dt   = 0.0d0  ! time width
real(8),parameter:: timemax=25.0d0 ! simulation end time

! option
integer, parameter :: flag_HDC = 0 ! 1 --> HDC on , 0 --> HDC off
integer, parameter :: flag_flux = 2 ! 1 (HLL), 2 (HLLD)

! coordinate 
integer,parameter::nx=50        ! the number of grids in the simulation box
integer,parameter::ny=150 ! the number of grids in the simulation box
integer,parameter::nz=1          ! the number of grids in the simulation box
integer,parameter::ngh=2         ! the number of ghost cells
integer,parameter::nxtot=nx+2*ngh+1 ! the total number of grids including ghost cells
integer,parameter::nytot=ny+2*ngh+1 ! the total number of grids including ghost cells
integer,parameter::nztot=1 ! the total number of grids including ghost cells
integer,parameter::is=ngh+1         ! the index of the leftmost grid
integer,parameter::js=ngh+1         ! the index of the leftmost grid
integer,parameter::ks=1         ! the index of the leftmost grid
integer,parameter::ie=nx+ngh     ! the index of the rightmost grid
integer,parameter::je=ny+ngh     ! the index of the rightmost grid
integer,parameter::ke=1     ! the index of the rightmost grid
real(8),parameter::xmin=-0.25d0,xmax=0.25d0
real(8),parameter::ymin=-0.75d0,ymax=0.75d0
real(8),parameter::zmin=0.0d0,zmax=1.0d0

real(8),parameter::Ccfl=0.4d0
real(8),parameter::gam=5.0d0/3.0d0 !! adiabatic index

real(8) ch       ! advection speed of divergence B
real(8), parameter :: alpha = 0.1d0    ! decay timescale of divergence B

real(8),parameter::grav_accy=-0.1d0  ! gravitaional acceleration


! indices of the conservative variables
integer, parameter :: IDN = 1
integer, parameter :: IMX = 2
integer, parameter :: IMY = 3
integer, parameter :: IMZ = 4
integer, parameter :: IPR = 5
integer, parameter :: IBX = 6
integer, parameter :: IBY = 7
integer, parameter :: IBZ = 8
integer, parameter :: IPS = 9
integer, parameter :: NVAR = 9
integer, parameter :: NFLX = 9

! indices of the primitive variables
integer, parameter :: IVX = 2
integer, parameter :: IVY = 3
integer, parameter :: IVZ = 4
integer, parameter :: IEN = 5



! definition of arrays 
real(8),dimension(nxtot)::xf,xv
real(8),dimension(nytot)::yf,yv
real(8),dimension(nztot)::zf,zv
real(8),dimension(NVAR,nxtot,nytot,nztot) :: Uo
real(8),dimension(NVAR,nxtot,nytot,nztot) :: U
real(8),dimension(NVAR,nxtot,nytot,nztot) :: Q
real(8),dimension(NVAR,nxtot,nytot,nztot) :: F
real(8),dimension(NVAR,nxtot,nytot,nztot) :: G
real(8),dimension(NVAR,nxtot,nytot,nztot) :: H

! output 
character(20),parameter::dirname="nohdc" ! directory name

! snapshot
integer, parameter :: unitsnap = 17
real(8), parameter:: dtsnap=2.0d-1

! realtime analysis 
integer, parameter :: nevo = 2
integer, parameter :: unitevo =11
real(8) :: phys_evo(nevo)

!logical :: flag_binary = .False.
logical :: flag_binary = .true.


integer :: i, j,k
real(8) :: stime

      ! make the directory for output
      call makedirs(trim(dirname))

      write(6,*) "setup grids and initial condition"
      call GenerateGrid(xf, xv, yf, yv, zf, zv)
      call GenerateProblem(xv, yv, zv, Q)
      call PrIMYConsrv(Q, U)
      call BoundaryCondition(xf,yf,Q)
      call Output( .TRUE., flag_binary, dirname, xf, xv, yf, yv, Q )


  write(6,*) "Start the simulation"
  open(unitevo,file=trim(dirname)//'/'//'ana.dat', action="write")
! main loop
      ntime = 1
      do !ntime=1,ntimemax
         dt = TimestepControl(xf, yf, zf, Q)
         if( time + dt > timemax ) dt = timemax - time

         Uo(:,:,:,:) = U(:,:,:,:)

         call NumericalFlux( Q, F, G, H )
         call UpdateConsv( 0.5d0*dt, xf, yf, zf, F, G, H, Q, U, U )
         call Consrv2Prim( U, Q )
         call BoundaryCondition(xf, yf,Q )

         call NumericalFlux( Q, F, G, H )
         call UpdateConsv( dt, xf, yf, zf, F, G, H, Q, Uo, U )
         call Consrv2Prim( U, Q )
         call BoundaryCondition( xf, yf, Q )

         time=time+dt
         ntime = ntime+1
         call Output( .FALSE., flag_binary, dirname, xf, xv, yf, yv, Q)
!         call Output( .true., xf, xv, yf, yv, Q)

         print*, "ntime = ",ntime, "time = ",time, dt

         if( mod(ntime,10) .eq. 0 ) then
             call RealtimeAnalysis(xv,yv,Q,phys_evo)
             write(unitevo,*) time, phys_evo(1:nevo)
         endif

         if(time >= timemax) exit 
      enddo 

      close(unitevo)


!      write(6,*) "program has been finished"
contains
!-------------------------------------------------------------------
!       Generate coordiantes
!       xf,yf,zf --> cell boundary xf(i) <==> x_{i-1/2}
!       xv,yv,zv --> cell center   xv(i) <==> x_{i}
!-------------------------------------------------------------------
subroutine GenerateGrid(xf, xv, yf, yv, zf, zv)
implicit none
real(8), intent(out) :: xf(:), xv(:)
real(8), intent(out) :: yf(:), yv(:)
real(8), intent(out) :: zf(:), zv(:)
real(8) :: dx,dy
integer::i,j


      dx=(xmax-xmin)/dble(nx)
      do i=1,nxtot
         xf(i) = dx*(i-(ngh+1))+xmin
      enddo
      do i=1,nxtot-1
         xv(i) = 0.5d0*(xf(i+1)+xf(i))
      enddo

      dy=(ymax-ymin)/dble(ny)
      do j=1,nytot
         yf(j) = dx*(j-(ngh+1))+ymin
      enddo
      do j=1,nytot-1
         yv(j) = 0.5d0*(yf(j+1)+yf(j))
      enddo

return
end subroutine GenerateGrid
!-------------------------------------------------------------------
!       Generate initial condition of the primitive variables
!-------------------------------------------------------------------
subroutine GenerateProblem(xv, yv, zv, Q )
implicit none
integer::i, j, k
real(8), intent(in ) :: xv(:), yv(:), zv(:)
real(8), intent(out) :: Q(:,:,:,:)
real(8) :: pi, den, B0

      pi = dacos(-1.0d0)
      B0 = 0.5d0*sqrt( abs(grav_accy)/(2.0*2.0d0*pi) )

      do k=ks,ke
      do j=js,je
      do i=is,ie
           if( yv(j) .lt. 0.0d0 ) then
               den = 1.0d0
           else 
               den = 2.0d0
           endif
           Q(IDN,i,j,k) = den
           Q(IVX,i,j,k) = 0.0d0
           Q(IVY,i,j,k) = 0.0d0
           Q(IVZ,i,j,k) = 0.0d0
           Q(IBX,i,j,k) = B0 
           Q(IBY,i,j,k) = 0.0d0
           Q(IBZ,i,j,k) = 0.0d0
           Q(IPR,i,j,k) = 2.5d0 + grav_accy*den*yv(j)

           Q(IVY,i,j,k)= 0.01d0/4.0d0 &
                     & *(-dcos(2.0d0*pi*(xv(i)-(xmax+xmin)/2.0d0)/(xmax-xmin))) &
                     & *(1.0+cos(2.0d0*pi*(yv(j)-(ymax+ymin)/2.0d0)/(ymax-ymin)))
!                     & *dexp( - (xv(i) - (x1max+x1min)/2.0d0)**2/(0.1**2) )
!                     & *(+cos(2.0d0*pi*(xv(i)-(x1max+x1min)/2.0d0)/(x1max-x1min)))
      enddo
      enddo
      enddo

return
end subroutine GenerateProblem

!-------------------------------------------------------------------
!       Boundary Condition of the primitive variables
!-------------------------------------------------------------------
subroutine BoundaryCondition(xf, yf, Q)
implicit none
real(8), intent(inout) :: xf(:), yf(:), Q(:,:,:,:)
integer::i,j,k,ish

      do k=ks,ke
      do j=1,nytot-1
      do i=1,ngh
          Q(:,is-i,j,k)  = Q(:,ie+1-i,j,k)
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=1,nytot-1
      do i=1,ngh

          Q(:,ie+i,j,k)= Q(:,is+i-1,j,k)
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=1,ngh
      do i=1,nxtot-1
          Q(IDN,i,js-j,k)  = Q(IDN,i,js-1+j,k)
          Q(IVX,i,js-j,k)  = Q(IVX,i,js-1+j,k)
          Q(IVY,i,js-j,k)  = -Q(IVY,i,js-1+j,k)
          Q(IVZ,i,js-j,k)  = Q(IVZ,i,js-1+j,k)
          Q(IPR,i,js-j,k)  = Q(IPR,i,js-1+j,k) &
                           - Q(IDN,i,js-1+j,k)*grav_accy*(2*j-1)*(yf(j+1)-yf(j))
          Q(IBX,i,js-j,k)  = Q(IBX,i,js-1+j,k)
          Q(IBY,i,js-j,k)  = Q(IBY,i,js-1+j,k)
          Q(IBZ,i,js-j,k)  = Q(IBZ,i,js-1+j,k)
          Q(IPS,i,js-j,k)  = Q(IPS,i,js-1+j,k)
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=1,ngh
      do i=1,nxtot-1
          Q(IDN,i,je+j,k) = Q(IDN,i,je-j+1,k)
          Q(IVX,i,je+j,k) = Q(IVX,i,je-j+1,k)
          Q(IVY,i,je+j,k) = Q(IVY,i,je-j+1,k)
          Q(IVZ,i,je+j,k) = Q(IVZ,i,je-j+1,k)
          Q(IPR,i,je+j,k) = Q(IPR,i,je-j+1,k) &
                          + Q(IDN,i,je-j+1,k)*grav_accy*(2*j-1)*(yf(j+1)-yf(j))
          Q(IBX,i,je+j,k) = Q(IBX,i,je-j+1,k)
          Q(IBY,i,je+j,k) = Q(IBY,i,je-j+1,k)
          Q(IBZ,i,je+j,k) = Q(IBZ,i,je-j+1,k)
          Q(IPS,i,je+j,k) = Q(IPS,i,je-j+1,k)
      enddo
      enddo
      enddo

return
end subroutine BoundaryCondition
!-------------------------------------------------------------------
!       Primitive variables ===> Conservative variables
!       Input  : Q
!       Output : U
!-------------------------------------------------------------------
subroutine PrIMYConsrv(Q, U)
implicit none
real(8), intent(in) :: Q(:,:,:,:)
real(8), intent(out) :: U(:,:,:,:)
integer::i,j,k

      do k=ks,ke
     !$omp parallel do private(i,j)
      do j=js,je
      do i=is,ie
          U(IDN,i,j,k) = Q(IDN,i,j,k)
          U(IMX,i,j,k) = Q(IDN,i,j,k)*Q(IVX,i,j,k)
          U(IMY,i,j,k) = Q(IDN,i,j,k)*Q(IVY,i,j,k)
          U(IMZ,i,j,k) = Q(IDN,i,j,k)*Q(IVZ,i,j,k)
          U(IEN,i,j,k) = 0.5d0*Q(IDN,i,j,k)*( Q(IVX,i,j,k)**2 + Q(IVY,i,j,k)**2 + Q(IVZ,i,j,k)**2 ) &
                       + 0.5d0*( Q(IBX,i,j,k)**2 + Q(IBY,i,j,k)**2 + Q(IBZ,i,j,k)**2 ) &
                       + Q(IPR,i,j,k)/(gam - 1.0d0)
          U(IBX,i,j,k) = Q(IBX,i,j,k)
          U(IBY,i,j,k) = Q(IBY,i,j,k)
          U(IBZ,i,j,k) = Q(IBZ,i,j,k)
          U(IPS,i,j,k) = Q(IPS,i,j,k)
      enddo
      enddo
      !$omp end parallel do
      enddo
      
return
end subroutine PrIMYConsrv
!-------------------------------------------------------------------
!       Conservative variables ===> Primitive variables
!       Input  : U
!       Output : Q
!-------------------------------------------------------------------
subroutine Consrv2Prim( U, Q )
implicit none
real(8), intent(in) :: U(:,:,:,:)
real(8), intent(out) :: Q(:,:,:,:)
integer::i,j,k
real(8) :: inv_d;
real(8) :: stime
!  stime = omp_get_wtime()

      do k=ks,ke
     !$omp parallel do private( i, j, inv_d )
      do j=js,je
      do i=is,ie
           Q(IDN,i,j,k) = U(IDN,i,j,k)
           inv_d = 1.0d0/U(IDN,i,j,k)
           Q(IVX,i,j,k) = U(IMX,i,j,k)*inv_d
           Q(IVY,i,j,k) = U(IMY,i,j,k)*inv_d
           Q(IVZ,i,j,k) = U(IMZ,i,j,k)*inv_d
           Q(IPR,i,j,k) = ( U(IEN,i,j,k) &
                        - 0.5d0*(U(IMX,i,j,k)**2 + U(IMY,i,j,k)**2 + U(IMZ,i,j,k)**2)*inv_d  &
                        - 0.5d0*(U(IBX,i,j,k)**2 + Q(IBY,i,j,k)**2 + Q(IBZ,i,j,k)**2) )*(gam-1.0d0)
           Q(IBX,i,j,k) = U(IBX,i,j,k)
           Q(IBY,i,j,k) = U(IBY,i,j,k)
           Q(IBZ,i,j,k) = U(IBZ,i,j,k)
           Q(IPS,i,j,k) = U(IPS,i,j,k)
      enddo
      enddo
      !$omp end parallel do
      enddo
!      print*,omp_get_wtime() - stime

return
end subroutine Consrv2Prim
!-------------------------------------------------------------------
!       determine dt 
!-------------------------------------------------------------------
real(8) function TimestepControl(xf, yf, zf, Q)
implicit none
real(8), intent(in) :: xf(:), yf(:), zf(:), Q(:,:,:,:)
real(8)::dtl1
real(8)::dtl2
real(8)::dtl3
real(8)::dtlocal
real(8)::dtmin,cf
integer::i,j,k

      dtmin=1.0d90

      do k=ks,ke
     !$omp parallel do private(i,dtl1,dtl2,cf) reduction (min: dtmin)
      do j=js,je
      do i=is,ie
         cf = dsqrt( (gam*Q(IPR,i,j,k) + Q(IBX,i,j,k)**2 + Q(IBY,i,j,k)**2 + Q(IBZ,i,j,k)**2)/Q(IDN,i,j,k))
         dtl1 =(xf(i+1)-xf(i))/(abs(Q(IVX,i,j,k)) + cf)
         dtl2 =(yf(j+1)-yf(j))/(abs(Q(IVY,i,j,k)) + cf)
!         dtl3 =(zf(j+1)-zf(j))/(abs(Q(IVZ,i,j,k)) + cf)
         dtmin = min(dtl1,dtl2,dtmin)
!         dtlocal = min(dtl1,dtl2)
!         if(dtlocal .lt. dtmin) dtmin = dtlocal
      enddo
      enddo
      !$omp end parallel do
      enddo

      TimestepControl = Ccfl* dtmin
!      write(6,*)"dt",dt

return
end function TimestepControl

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
subroutine NumericalFlux( Q, F, G, H )
implicit none
integer::i,j,k
real(8), intent(in) :: Q(:,:,:,:)
real(8), intent(out) :: F(:,:,:,:)
real(8), intent(out) :: G(:,:,:,:)
real(8), intent(out) :: H(:,:,:,:)
real(8),dimension(NFLX,nxtot,nytot,nztot):: Ql,Qr
real(8),dimension(NFLX):: flx
real(8) :: dQm(NFLX), dQp(NFLX), dQmon(NFLX)
real(8) :: ddmon, dvmon, dpmon

      ch = 1.0d0*Ccfl*min( xf(is+1) - xf(is), yf(js+1) - yf(js ) )/dt

      ! numerical flux in the x direction
!$omp parallel 
      do k=ks,ke
      !$omp do private( i, dQp, dQm, dQmon )
      do j=js,je
      do i=is-1,ie+1
         dQp(1:NVAR) = Q(1:NVAR,i+1,j,k) - Q(1:NVAR,i  ,j,k)
         dQm(1:NVAR) = Q(1:NVAR,i  ,j,k) - Q(1:NVAR,i-1,j,k)

         call vanLeer(NFLX, dQp, dQm, dQmon)

         ! Ql(i,j,k) --> W_(i-1/2,j,k)
         ! Qr(i,j,k) --> W_(i-1/2,j,k)
         Ql(1:NVAR,i+1,j,k) = Q(1:NVAR,i,j,k) + 0.5d0*dQmon(1:NVAR)
         Qr(1:NVAR,i  ,j,k) = Q(1:NVAR,i,j,k) - 0.5d0*dQmon(1:NVAR)

      enddo
      enddo
      !$omp end do
      enddo

      if( flag_flux == 1 ) then
          do k=ks,ke
         !$omp do private( i, flx )
          do j=js,je
          do i=is,ie+1
            call HLL(1,Ql(:,i,j,k),Qr(:,i,j,k),flx)
    
!             flx(IBX) = flag_HDC*0.5d0*(Ql(IPS,i,j,k) + Qr(IPS,i,j,k) - ch*(Qr(IBX,i,j,k) - Ql(IBX,i,j,k)) )
!             flx(IPS) = flag_HDC*0.5d0*(ch*ch*(Ql(IBX,i,j,k) + Qr(IBX,i,j,k)) - ch*(Qr(IPS,i,j,k) - Ql(IPS,i,j,k)) )
!    
!!             if (flx(IDN) >= 0.0) then 
!                 flx(ISC) = flx(IDN)*Ql(ISC,i,j,k);
!             else 
!                 flx(ISC) = flx(IDN)*Qr(ISC,i,j,k);
!             endif
    
             F(:,i,j,k)  = flx(:)
          enddo
          enddo
          !$omp end do
          enddo
      else if ( flag_flux == 2 ) then
          do k=ks,ke
         !$omp do private( i, flx )
          do j=js,je
          do i=is,ie+1
             call HLLD(1,Ql(:,i,j,k),Qr(:,i,j,k),flx)
    
!             flx(IBX) = flag_HDC*0.5d0*(Ql(IPS,i,j,k) + Qr(IPS,i,j,k) - ch*(Qr(IBX,i,j,k) - Ql(IBX,i,j,k)) )
!             flx(IPS) = flag_HDC*0.5d0*(ch*ch*(Ql(IBX,i,j,k) + Qr(IBX,i,j,k)) - ch*(Qr(IPS,i,j,k) - Ql(IPS,i,j,k)) )
!    
!             if (flx(IDN) >= 0.0) then 
!                 flx(ISC) = flx(IDN)*Ql(ISC,i,j,k);
!             else 
!                 flx(ISC) = flx(IDN)*Qr(ISC,i,j,k);
!             endif
    
             F(:,i,j,k)  = flx(:)
          enddo
          enddo
          !$omp end do
          enddo
      endif

      ! numerical flux in the y direction
      do k=ks,ke
      !$omp do private( i, dQp, dQm, dQmon )
      do j=js-1,je+1
      do i=is,ie
         dQp(1:NVAR) = Q(1:NVAR,i,j+1,k) - Q(1:NVAR,i,j  ,k)
         dQm(1:NVAR) = Q(1:NVAR,i,j  ,k) - Q(1:NVAR,i,j-1,k)

         call vanLeer(NFLX, dQp, dQm, dQmon)

         ! Ql(i,j,k) --> W_(i-1/2,j,k)
         ! Qr(i,j,k) --> W_(i-1/2,j,k)
         Ql(1:NVAR,i,j+1,k) = Q(1:NVAR,i,j,k) + 0.5d0*dQmon(1:NVAR)
         Qr(1:NVAR,i,j  ,k) = Q(1:NVAR,i,j,k) - 0.5d0*dQmon(1:NVAR)

      enddo
      enddo
      !$omp end do
      enddo


      if( flag_flux == 1 ) then
          do k=ks,ke
         !$omp do private( i, flx )
          do j=js,je+1
          do i=is,ie
              call HLL(2,Ql(:,i,j,k),Qr(:,i,j,k),flx)

              G(:,i,j,k)  = flx(:)
          enddo
          enddo
          !$omp end do
          enddo
      else if( flag_flux == 2 ) then
          do k=ks,ke
         !$omp do private( i, flx )
          do j=js,je+1
          do i=is,ie
              call HLLD(2,Ql(:,i,j,k),Qr(:,i,j,k),flx)

              G(:,i,j,k)  = flx(:)
          enddo
          enddo
          !$omp end do
          enddo
      endif

!$omp end parallel
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
subroutine HLL(idir,Ql,Qr,flx)
implicit none
integer, intent(in) :: idir
real(8),intent(in)  ::Ql(:), Qr(:)
real(8),intent(out) :: flx(:)
integer :: IVpara, IVperp1, IVperp2
integer :: IBpara, IBperp1, IBperp2
real(8):: b1
real(8):: Ul(NFLX), Ur(NFLX)
real(8):: Fl(NFLX), Fr(NFLX)
real(8):: Ust(NFLX)
real(8):: Fst(NFLX)
real(8):: cfl,cfr
real(8):: sl, sr
real(8):: pbl, pbr, ptotl, ptotr
integer :: i, n

      if( idir == 1 ) then
           IVpara  = IVX
           IVperp1 = IVY
           IVperp2 = IVZ
           IBpara  = IBX
           IBperp1 = IBY
           IBperp2 = IBZ
      else if (idir == 2 ) then
           IVpara  = IVY
           IVperp1 = IVZ
           IVperp2 = IVX
           IBpara  = IBY
           IBperp1 = IBZ
           IBperp2 = IBX
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

          sl = min(Ql(IVpara) - cfl,Qr(IVpara) - cfr)
          sr = max(Ql(IVpara) + cfl,Qr(IVpara) + cfr)

          if( sl > 0.0d0 ) then
               flx(:) = Fl(:)
          else if (sr <= 0.0d0 ) then
               flx(:) = Fr(:)
          else 
               flx(:)  = (sr*Fl(:) - sl*Fr(:) + sl*sr*( Ur(:) - Ul(:) ))/(sr - sl)
          endif

          flx(IBpara) = flag_HDC*0.5d0*(Ql(IPS) + Qr(IPS) - ch*(Qr(IBpara) - Ql(IBpara)) )
          flx(IPS) = flag_HDC*0.5d0*(ch*ch*(Ql(IBpara) + Qr(IBpara)) - ch*(Qr(IPS) - Ql(IPS)) )
    

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
subroutine HLLD(idir,Ql,Qr,flx)
implicit none
integer, intent(in) :: idir
real(8),intent(in)  ::Ql(:), Qr(:)
real(8),intent(out) :: flx(:)
integer :: IVpara, IVperp1, IVperp2
integer :: IBpara, IBperp1, IBperp2
real(8):: b1
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
integer :: i, n

      if( idir == 1 ) then
           IVpara  = IVX
           IVperp1 = IVY
           IVperp2 = IVZ
           IBpara  = IBX
           IBperp1 = IBY
           IBperp2 = IBZ
      else if (idir == 2 ) then
           IVpara  = IVY
           IVperp1 = IVZ
           IVperp2 = IVX
           IBpara  = IBY
           IBperp1 = IBZ
           IBperp2 = IBX
      endif
          b1 = 0.5d0*( Ql(IBpara) + Qr(IBpara) )
          
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
!
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

!          if( sqrtdr .ne. sqrtdr ) then
!              print*, "sqrtdr",sqrtdr, Cr, Cmr
!             print*,"S", S0,S2,S4
!              stop
!          endif

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
!         print*,test(IVperp2)
         

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

             flx(IBpara) = flag_HDC*0.5d0*(Ql(IPS) + Qr(IPS) - ch*(Qr(IBpara) - Ql(IBpara)) )
             flx(IPS) = flag_HDC*0.5d0*(ch*ch*(Ql(IBpara) + Qr(IBpara)) - ch*(Qr(IPS) - Ql(IPS)) )
    
return
end subroutine HLLD
!-------------------------------------------------------------------
!       Update consevative variables U using numerical flux F
!-------------------------------------------------------------------
subroutine UpdateConsv( dt1, xf, yf, zf, F, G, H, Q, Uo, U)
real(8), intent(in) :: dt1
real(8), intent(in)  :: xf(:), yf(:), zf(:)
real(8), intent(in)  :: F(:,:,:,:), G(:,:,:,:), H(:,:,:,:)
real(8), intent(in)  :: Uo(:,:,:,:), Q(:,:,:,:)
real(8), intent(out) :: U(:,:,:,:)
real(8) :: divB , src
integer::i,n,j,k

!$omp parallel 
      do k=ks,ke
      !$omp do private( i,j )
      do j=js,je
      do i=is,ie
         U(:,i,j,k) = Uo(:,i,j,k) + dt1*(- F(:,i+1,j,k) + F(:,i,j,k))/(xf(i+1)-xf(i)) &
                                  + dt1*(- G(:,i,j+1,k) + G(:,i,j,k))/(yf(j+1)-yf(j)) 
      enddo
      enddo
      !$omp end do
      enddo

      ! Source term
      do k=ks,ke
      !$omp do private( i, src )
      do j=js,je
      do i=is,ie
         src = dt1*Q(IDN,i,j,k)*grav_accy
         U(IMY,i,j,k) = U(IMY,i,j,k) + src
         U(IEN,i,j,k) = U(IEN,i,j,k) + src*Q(IVY,i,j,k)

         U(IPS,i,j,k) = U(IPS,i,j,k)*dexp(-alpha*Ccfl*dt1/dt)
      enddo
      enddo
      !$omp end do
      enddo

!$omp end parallel

return
end subroutine UpdateConsv
!-------------------------------------------------------------------
!       Output snapshot files 
!       Input  : flag, dirname, xf, xv, Q
!
!       flag = .true.  --> output snapshot when calling this subroutine
!       flag = .false. --> output snapshot every dtsnap
!-------------------------------------------------------------------
subroutine Output( flag, flag_binary, dirname, xf, xv, yf, yv, Q )
logical, intent(in) :: flag ! false --> output per dtsnap, true --> force to output
logical, intent(in) :: flag_binary ! false --> ascii, true --> binary
character(20), intent(in) :: dirname 
real(8), intent(in) :: xf(:), xv(:), yf(:), yv(:), Q(:,:,:,:)
integer::i,j,k
character(100)::filename
real(8), save :: tsnap = - dtsnap
integer, save :: nsnap

    if( .not.flag) then
        if( time + 1.0d-14.lt. tsnap+dtsnap) return
    endif

        write(filename,'(i5.5)') nsnap
    if( flag_binary ) then
        filename = trim(dirname)//"/snap"//trim(filename)//".bin"
        open(unitsnap,file=filename,form='unformatted',access="stream",action="write")
        write(unitsnap) time
        write(unitsnap) nx
        write(unitsnap) ny
        write(unitsnap) 5
        write(unitsnap) NVAR - 5
        write(unitsnap) xv(is:ie)
        write(unitsnap) yv(js:je)
        write(unitsnap) real(Q(1:5,is:ie,js:je,ks:ke)) ! single precision
        write(unitsnap) real(Q(6:NVAR,is:ie,js:je,ks:ke)) ! single precision
        close(unitsnap)
    else 
        filename = trim(dirname)//"/snap"//trim(filename)//".dat"
        open(unitsnap,file=filename,form='formatted',action="write")
        write(unitsnap,*) "# time = ",time
        write(unitsnap,*) "#nx, ny = ", nx, ny
          do k=ks,ke
          do j=js,je
          do i=is,ie
              write(unitsnap,*) xv(i), yv(j), Q(IDN,i,j,k), Q(IVX,i,j,k), Q(IVY,i,j,k), Q(IVZ,i,j,k), &
                                Q(IPR,i,j,k), Q(IBX,i,j,k), Q(IBY,i,j,k), Q(IBZ,i,j,k), Q(IPS,i,j,k) 
          enddo
          enddo
          enddo
          close(unitsnap)
      endif

    write(6,*) "output binary file:  ",filename,time

    nsnap=nsnap+1
    tsnap=tsnap + dtsnap

return
end subroutine Output
!-------------------------------------------------------------------
!       create directory
!       Input  : the directory to be created
!-------------------------------------------------------------------
subroutine makedirs(outdir)
implicit none
character(len=*), intent(in) :: outdir
character(len=256) command
write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
!      write(*, *) trim(command)
call system(command)
end subroutine makedirs
!-------------------------------------------------------------------
!       Realtime Analysis
!       Input  : xf, xv
!       Output : phys_evo(nevo)
!-------------------------------------------------------------------
subroutine RealtimeAnalysis(xv,yv,Q,phys_evo)
implicit none
real(8), intent(in) :: xv(:), yv(:), Q(:,:,:,:)
real(8), intent(out) :: phys_evo(:)
integer::i,j,k
real(8) :: dvx, er_divB 

      dvx = 0.0d0
      er_divB = 0.0d0
      do k=ks,ke
      do j=js,je
      do i=is,ie
           dvx = dvx + Q(IVX,i,j,k)**2
           er_divB = er_divB + ( Q(IBX,i+1,j,k) - Q(IBX,i-1,j,k) &
                               + Q(IBY,i,j+1,k) - Q(IBY,i,j-1,k) )**2 &
                       /( Q(IBX,i,j,k)**2 + Q(IBY,i,j,k)**2 )
      enddo
      enddo
      enddo
      phys_evo(1) = sqrt(dvx/dble(nx*ny))
      phys_evo(2) = sqrt(er_divB/dble(nx*ny))
      
return
end subroutine

end program main
