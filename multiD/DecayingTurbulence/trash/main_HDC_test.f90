program main
implicit none

integer::ntime                        ! counter of the timestep
integer,parameter::ntimemax=200000     ! the maximum timesteps
real(8)::time,dt                    ! time, timewidth
data time / 0.0d0 /
real(8),parameter:: timemax=20.0d0   
real(8),parameter:: dtout=0.1d0

integer, parameter :: flag_HDC = 1 ! 1 --> HDC on , 0 --> HDC off
integer, parameter :: flag_flux = 3 ! 1 (HLL), 2 (HLLC), 3 (HLLD)

integer,parameter::nx=128*1*4        ! the number of grids in the simulation box
integer,parameter::ny=128*2 *4 ! the number of grids in the simulation box
integer,parameter::nz=1          ! the number of grids in the simulation box
integer,parameter::mgn=2         ! the number of ghost cells
integer,parameter::in=nx+2*mgn+1 ! the total number of grids including ghost cells
integer,parameter::jn=ny+2*mgn+1 ! the total number of grids including ghost cells
integer,parameter::kn=1 ! the total number of grids including ghost cells
integer,parameter::is=mgn+1         ! the index of the leftmost grid
integer,parameter::js=mgn+1         ! the index of the leftmost grid
integer,parameter::ks=1         ! the index of the leftmost grid
integer,parameter::ie=nx+mgn     ! the index of the rightmost grid
integer,parameter::je=ny+mgn     ! the index of the rightmost grid
integer,parameter::ke=1     ! the index of the rightmost grid
real(8),parameter::x1min=-0.5d0,x1max=0.5d0
real(8),parameter::x2min=-1.0d0,x2max=1.0d0
real(8),parameter::x3min=0.0d0,x3max=1.0d0

real(8),parameter::Ccfl=0.4d0
real(8),parameter::gam=5.0d0/3.0d0 !! adiabatic index

real(8) ch       ! advection speed of divergence B
real(8), parameter :: alpha = 0.1d0    ! decay timescale of divergence B

integer, parameter :: IDN = 1
integer, parameter :: IM1 = 2
integer, parameter :: IM2 = 3
integer, parameter :: IM3 = 4
integer, parameter :: IPR = 5
integer, parameter :: IB1 = 6
integer, parameter :: IB2 = 7
integer, parameter :: IB3 = 8
integer, parameter :: IPS = 9
integer, parameter :: ISC = 10
integer, parameter :: NVAR = 10
integer, parameter :: NFLX = 10

integer, parameter :: IV1 = 2
integer, parameter :: IV2 = 3
integer, parameter :: IV3 = 4
integer, parameter :: IEN = 5


integer, parameter :: nevo = 2
integer, parameter :: unitevo =11
real(8) :: phys_evo(nevo)

real(8),dimension(in)::xf,xv
real(8),dimension(jn)::yf,yv
real(8),dimension(kn)::zf,zv
real(8),dimension(NVAR,in,jn,kn) :: Uo
real(8),dimension(NVAR,in,jn,kn) :: U
real(8),dimension(NVAR,in,jn,kn) :: Q
real(8),dimension(NVAR,in,jn,kn) :: F
real(8),dimension(NVAR,in,jn,kn) :: G
real(8),dimension(NVAR,in,jn,kn) :: H
integer :: i, j,k

!      write(6,*) "setup grids and initial condition"
      call GenerateGrid(xf, xv, yf, yv, zf, zv)
      call GenerateProblem(xv, yv, zv, Q)
      call ConsvVariable(Q, U)
      call BoundaryCondition(xf,yf,Q)
!      call Output( .TRUE., xf, xv, yf, yv, Q )




      open(unitevo,file="vy_evo_B0.8_dc1.dat",action="write")
! main loop
      ntime = 1
      mloop: do !ntime=1,ntimemax
         call TimestepControl(xf, yf, zf, Q)
         if( time + dt > timemax ) dt = timemax - time

         Uo(:,:,:,:) = U(:,:,:,:)

         call NumericalFlux( Q, F, G, H )
         call UpdateConsv( 0.5d0*dt, xf, yf, zf, F, G, H, Q, U, U )
         call PrimVariable( U, Q )
         call BoundaryCondition(xf, yf,Q )

         call NumericalFlux( Q, F, G, H )
         call UpdateConsv( dt, xf, yf, zf, F, G, H, Q, Uo, U )
         call PrimVariable( U, Q )
         call BoundaryCondition( xf, yf, Q )

         time=time+dt
         ntime = ntime+1
!         call Output( .FALSE., xf, xv, yf, yv, Q)
!         call Output( .true., xf, xv, yf, yv, Q)

         print*, "ntime = ",ntime, "time = ",time, dt

         if( mod(ntime,10) .eq. 0 ) then
             call Analysis(xv,yv,Q,phys_evo)
             write(unitevo,*) time, phys_evo(1:nevo)
             call flush(unitevo)
         endif

         if(time >= timemax) exit mloop
         if(ntime .eq. 10 )exit
      enddo mloop
      close(unitevo)
!      call Output( .TRUE., xf, xv, yf, yv, Q)

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
      do i=1,in
         xf(i) = dx*(i-(mgn+1))+x1min
      enddo
      do i=1,in-1
         xv(i) = 0.5d0*(xf(i+1)+xf(i))
      enddo

      dy=(x2max-x2min)/dble(ny)
      do j=1,jn
         yf(j) = dx*(j-(mgn+1))+x2min
      enddo
      do j=1,jn-1
         yv(j) = 0.5d0*(yf(j+1)+yf(j))
      enddo

      return
      end subroutine GenerateGrid

      subroutine GenerateProblem(xv, yv, zv, Q )
      implicit none
      integer::i, j, k
      real(8), intent(in ) :: xv(:), yv(:), zv(:)
      real(8), intent(out) :: Q(:,:,:,:)
      real(8) :: pi, den, B0, rho1, rho2, dv, wid, sig

      pi = dacos(-1.0d0)

      rho1 = 1.0d0
      rho2 = 1.0d0
      dv   = 2.00d0
      wid  = 0.05d0
      sig  = 0.2d0
      B0  = 0.0d0*dsqrt( dv**2*0.5d0*rho1*rho2/(rho1+rho2))
!      B0  = sqrt(2.0d0/3.0d0)
!      B0  = 1.0d0

      do k=ks,ke
      do j=js,je
      do i=is,ie
         Q(IDN,i,j,k) = 1.0d0 !+ 0.5d0*( dtanh( (yv(j)+0.25d0)/wid ) - tanh( (yv(j)-0.25d0)/wid) )
         Q(IV1,i,j,k)  = 0.5*dv*( dtanh( (yv(j)+0.5d0)/wid ) - dtanh( (yv(j) - 0.5d0)/wid ) - 1.0d0 )
         Q(IV2,i,j,k)  = 0.001d0*dsin(2.0d0*pi*xv(i))* &
             ( dexp( - (yv(j) + 0.5d0)**2/sig**2 ) +  &
               dexp( - (yv(j) - 0.5d0)**2/sig**2 ) )
         Q(IPR,i,j,k) = 1.0d0
         Q(IB1,i,j,k) = B0
         Q(IB2,i,j,k) = 0.0d0
         Q(IB3,i,j,k) = 0.0d0
         Q(IPS,i,j,k) = 0.0d0

         Q(ISC,i,j,k) = 0.5d0*( dtanh( (yv(j)+0.5d0)/wid ) - tanh( (yv(j)-0.5d0)/wid) )
      enddo
      enddo
      enddo

      return
      end subroutine GenerateProblem

      subroutine BoundaryCondition(xf, yf, Q)
      implicit none
      real(8), intent(inout) :: xf(:), yf(:), Q(:,:,:,:)
      integer::i,j,k,ish


      do k=ks,ke
      do j=1,jn-1
      do i=1,mgn

          Q(:,is-i,j,k)  = Q(:,ie+1-i,j,k)
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=1,jn-1
      do i=1,mgn

          Q(:,ie+i,j,k)= Q(:,is+i-1,j,k)
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=1,mgn
      do i=1,in-1
          Q(:,i,js-j,k)  = Q(:,i,je+1-j,k)
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=1,mgn
      do i=1,in-1
          Q(:,i,je+j,k)  = Q(:,i,js+j-1,k)
      enddo
      enddo
      enddo

      return
      end subroutine BoundaryCondition
!
      subroutine ConsvVariable(Q, U)
      implicit none
      real(8), intent(in) :: Q(:,:,:,:)
      real(8), intent(out) :: U(:,:,:,:)
      integer::i,j,k

      do k=ks,ke
      do j=js,je
      do i=is,ie
          U(IDN,i,j,k) = Q(IDN,i,j,k)
          U(IM1,i,j,k) = Q(IDN,i,j,k)*Q(IV1,i,j,k)
          U(IM2,i,j,k) = Q(IDN,i,j,k)*Q(IV2,i,j,k)
          U(IM3,i,j,k) = Q(IDN,i,j,k)*Q(IV3,i,j,k)
          U(IEN,i,j,k) = 0.5d0*Q(IDN,i,j,k)*( Q(IV1,i,j,k)**2 + Q(IV2,i,j,k)**2 + Q(IV3,i,j,k)**2 ) &
                       + 0.5d0*( Q(IB1,i,j,k)**2 + Q(IB2,i,j,k)**2 + Q(IB3,i,j,k)**2 ) &
                       + Q(IPR,i,j,k)/(gam - 1.0d0)
          U(IB1,i,j,k) = Q(IB1,i,j,k)
          U(IB2,i,j,k) = Q(IB2,i,j,k)
          U(IB3,i,j,k) = Q(IB3,i,j,k)
          U(IPS,i,j,k) = Q(IPS,i,j,k)
          U(ISC,i,j,k) = Q(IDN,i,j,k)*Q(ISC,i,j,k)
      enddo
      enddo
      enddo
      
      return
      end subroutine Consvvariable

      subroutine PrimVariable( U, Q )
      implicit none
      real(8), intent(in) :: U(:,:,:,:)
      real(8), intent(out) :: Q(:,:,:,:)
      integer::i,j,k
      real(8) :: inv_d;

      do k=ks,ke
      do j=js,je
      do i=is,ie
           Q(IDN,i,j,k) = U(IDN,i,j,k)
           inv_d = 1.0d0/U(IDN,i,j,k)
           Q(IV1,i,j,k) = U(IM1,i,j,k)*inv_d
           Q(IV2,i,j,k) = U(IM2,i,j,k)*inv_d
           Q(IV3,i,j,k) = U(IM3,i,j,k)*inv_d
           Q(IPR,i,j,k) = ( U(IEN,i,j,k) &
                        - 0.5d0*(U(IM1,i,j,k)**2 + U(IM2,i,j,k)**2 + U(IM3,i,j,k)**2)*inv_d  &
                        - 0.5d0*(U(IB1,i,j,k)**2 + Q(IB2,i,j,k)**2 + Q(IB3,i,j,k)**2) )*(gam-1.0d0)
           Q(IB1,i,j,k) = U(IB1,i,j,k)
           Q(IB2,i,j,k) = U(IB2,i,j,k)
           Q(IB3,i,j,k) = U(IB3,i,j,k)
           Q(IPS,i,j,k) = U(IPS,i,j,k)
           Q(ISC,i,j,k) = U(ISC,i,j,k)*inv_d
      enddo
      enddo
      enddo

      return
      end subroutine PrimVariable

      subroutine TimestepControl(xf, yf, zf, Q)
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
      do j=js,je
      do i=is,ie
         cf = dsqrt( (gam*Q(IPR,i,j,k) + Q(IB1,i,j,k)**2 + Q(IB2,i,j,k)**2 + Q(IB3,i,j,k)**2)/Q(IDN,i,j,k))
         dtl1 =(xf(i+1)-xf(i))/(abs(Q(IV1,i,j,k)) + cf)
         dtl2 =(yf(j+1)-yf(j))/(abs(Q(IV2,i,j,k)) + cf)
!         dtl3 =(zf(j+1)-zf(j))/(abs(Q(IV3,i,j,k)) + cf)
         dtlocal = min(dtl1,dtl2)
         if(dtlocal .lt. dtmin) dtmin = dtlocal
      enddo
      enddo
      enddo

      dt = Ccfl* dtmin
!      write(6,*)"dt",dt
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
      subroutine NumericalFlux( Q, F, G, H )
      implicit none
      integer::i,j,k
      real(8), intent(in) :: Q(:,:,:,:)
      real(8), intent(out) :: F(:,:,:,:)
      real(8), intent(out) :: G(:,:,:,:)
      real(8), intent(out) :: H(:,:,:,:)
      real(8),dimension(NFLX,in,jn,kn):: Ql,Qr
      real(8),dimension(NFLX):: flx
      real(8) :: dQm(NFLX), dQp(NFLX), dQmon(NFLX)
      real(8) :: ddmon, dvmon, dpmon

      ch = 1.0d0*Ccfl*min( xf(is+1) - xf(is), yf(js+1) - yf(js ) )/dt

      ! numerical flux in the x direction
      do k=ks,ke
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
      enddo

      if( flag_flux == 1 ) then
          do k=ks,ke
          do j=js,je
          do i=is,ie+1
            call HLL(1,Ql(:,i,j,k),Qr(:,i,j,k),flx)
    
             flx(IB1) = flag_HDC*0.5d0*(Ql(IPS,i,j,k) + Qr(IPS,i,j,k) - ch*(Qr(IB1,i,j,k) - Ql(IB1,i,j,k)) )
             flx(IPS) = flag_HDC*0.5d0*(ch*ch*(Ql(IB1,i,j,k) + Qr(IB1,i,j,k)) - ch*(Qr(IPS,i,j,k) - Ql(IPS,i,j,k)) )
    
             if (flx(IDN) >= 0.0) then 
                 flx(ISC) = flx(IDN)*Ql(ISC,i,j,k);
             else 
                 flx(ISC) = flx(IDN)*Qr(ISC,i,j,k);
             endif
    
             F(:,i,j,k)  = flx(:)
          enddo
          enddo
          enddo
      else if ( flag_flux == 2 ) then

      else if ( flag_flux == 3 ) then
          do k=ks,ke
          do j=js,je
          do i=is,ie+1
             call HLLD(1,Ql(:,i,j,k),Qr(:,i,j,k),flx)
    
             flx(IB1) = flag_HDC*0.5d0*(Ql(IPS,i,j,k) + Qr(IPS,i,j,k) - ch*(Qr(IB1,i,j,k) - Ql(IB1,i,j,k)) )
             flx(IPS) = flag_HDC*0.5d0*(ch*ch*(Ql(IB1,i,j,k) + Qr(IB1,i,j,k)) - ch*(Qr(IPS,i,j,k) - Ql(IPS,i,j,k)) )
    
             if (flx(IDN) >= 0.0) then 
                 flx(ISC) = flx(IDN)*Ql(ISC,i,j,k);
             else 
                 flx(ISC) = flx(IDN)*Qr(ISC,i,j,k);
             endif
    
             F(:,i,j,k)  = flx(:)
          enddo
          enddo
          enddo
      endif

      ! numerical flux in the y direction
      do k=ks,ke
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
      enddo

      if( flag_flux == 1 ) then
          do k=ks,ke
          do j=js,je+1
          do i=is,ie
              call HLL(2,Ql(:,i,j,k),Qr(:,i,j,k),flx)

              flx(IB2) = flag_HDC*0.5d0*(Ql(IPS,i,j,k) + Qr(IPS,i,j,k) - ch*(Qr(IB2,i,j,k) - Ql(IB2,i,j,k)) )
              flx(IPS) = flag_HDC*0.5d0*(ch*ch*(Ql(IB2,i,j,k) + Qr(IB2,i,j,k)) - ch*(Qr(IPS,i,j,k) - Ql(IPS,i,j,k)) )

              if (flx(IDN) >= 0.0) then 
                  flx(ISC) = flx(IDN)*Ql(ISC,i,j,k); 
              else 
                  flx(ISC) = flx(IDN)*Qr(ISC,i,j,k); 
              endif

              G(:,i,j,k)  = flx(:)
          enddo
          enddo
          enddo
      else if( flag_flux == 2 ) then

      else if( flag_flux == 3 ) then
          do k=ks,ke
          do j=js,je+1
          do i=is,ie
              call HLLD(2,Ql(:,i,j,k),Qr(:,i,j,k),flx)

              flx(IB2) = flag_HDC*0.5d0*(Ql(IPS,i,j,k) + Qr(IPS,i,j,k) - ch*(Qr(IB2,i,j,k) - Ql(IB2,i,j,k)) )
              flx(IPS) = flag_HDC*0.5d0*(ch*ch*(Ql(IB2,i,j,k) + Qr(IB2,i,j,k)) - ch*(Qr(IPS,i,j,k) - Ql(IPS,i,j,k)) )

              if (flx(IDN) >= 0.0) then 
                  flx(ISC) = flx(IDN)*Ql(ISC,i,j,k); 
              else 
                  flx(ISC) = flx(IDN)*Qr(ISC,i,j,k); 
              endif

              G(:,i,j,k)  = flx(:)
          enddo
          enddo
          enddo
      endif

      return
      end subroutine Numericalflux

!---------------------------------------------------------------------
!     HLL Riemann Solver
!---------------------------------------------------------------------
!     solve the HLL Riemann solver 
!
!     Input: Ql, Qr: primitive variables containing the perpendicular B fields 
!                    at the left and right states
!            1D array (IDN, IV1, IV2, IV3, IPR, IBperp1, IBperp2)
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
!            index: (IDN, IV1, IV2, IV3, IPR, IBperp1, IBperp2)
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
           IVpara  = IV1
           IVperp1 = IV2
           IVperp2 = IV3
           IBpara  = IB1
           IBperp1 = IB2
           IBperp2 = IB3
      else if (idir == 2 ) then
           IVpara  = IV2
           IVperp1 = IV3
           IVperp2 = IV1
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

          sl = min(Ql(IVpara) - cfl,Qr(IVpara) - cfr)
          sr = max(Ql(IVpara) + cfl,Qr(IVpara) + cfr)

          if( sl > 0.0d0 ) then
               flx(:) = Fl(:)
          else if (sr <= 0.0d0 ) then
               flx(:) = Fr(:)
          else 
               flx(:)  = (sr*Fl(:) - sl*Fr(:) + sl*sr*( Ur(:) - Ul(:) ))/(sr - sl)
          endif

      return
      end subroutine HLL

!---------------------------------------------------------------------
!     HLLD Riemann Solver
!---------------------------------------------------------------------
!     solve the HLL Riemann solver 
!
!     Input: Ql, Qr: primitive variables containing the perpendicular B fields 
!                    at the left and right states
!            1D array (IDN, IV1, IV2, IV3, IPR, IBperp1, IBperp2)
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
!            index: (IDN, IV1, IV2, IV3, IPR, IBperp1, IBperp2)
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
           IVpara  = IV1
           IVperp1 = IV2
           IVperp2 = IV3
           IBpara  = IB1
           IBperp1 = IB2
           IBperp2 = IB3
      else if (idir == 2 ) then
           IVpara  = IV2
           IVperp1 = IV3
           IVperp2 = IV1
           IBpara  = IB2
           IBperp1 = IB3
           IBperp2 = IB1
      endif
          b1 = 0.5d0*( Ql(IBpara) + Qr(IBpara) )
          
          pbl = 0.5d0*(b1**2 + Ql(IBperp1)**2 + Ql(IBperp2)**2)
          pbr = 0.5d0*(b1**2 + Qr(IBperp1)**2 + Qr(IBperp2)**2)
          ptotl = Ql(IPR) + pbl
          ptotr = Qr(IPR) + pbr

!          cfl = dsqrt( 0.5d0*( 2.0d0*pbl + gam*Ql(IPR) &
!                     + dsqrt( (2.0d0*pbl - gam*Ql(IPR))**2 &
!                     + 4.0d0*gam*Ql(IPR)*( Ql(IBperp1)**2 + Ql(IBperp2)**2 ) ) )/Ql(IDN) )
!          cfr = dsqrt( 0.5d0*( 2.0d0*pbr + gam*Qr(IPR) &
!                     + dsqrt( (2.0d0*pbr - gam*Qr(IPR))**2 &
!                     + 4.0d0*gam*Qr(IPR)*( Qr(IBperp1)**2 + Qr(IBperp2)**2 ) ) )/Qr(IDN) )
          cfl = dsqrt( (gam*Ql(IPR) + Ql(IBperp1)**2 + Ql(IBperp2)**2 + b1**2)/Ql(IDN))
          cfr = dsqrt( (gam*Qr(IPR) + Qr(IBperp1)**2 + Qr(IBperp2)**2 + b1**2)/Qr(IDN))

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

      return
      end subroutine HLLD
!!=====================================================================

      subroutine UpdateConsv( dt1, xf, yf, zf, F, G, H, Q, Uo, U)
      real(8), intent(in) :: dt1
      real(8), intent(in)  :: xf(:), yf(:), zf(:)
      real(8), intent(in)  :: F(:,:,:,:), G(:,:,:,:), H(:,:,:,:)
      real(8), intent(in)  :: Uo(:,:,:,:), Q(:,:,:,:)
      real(8), intent(out) :: U(:,:,:,:)
      real(8) :: divB , src
      integer::i,n,j,k

      do k=ks,ke
      do j=js,je
      do i=is,ie
         U(:,i,j,k) = Uo(:,i,j,k) + dt1*(- F(:,i+1,j,k) + F(:,i,j,k))/(xf(i+1)-xf(i)) &
                                  + dt1*(- G(:,i,j+1,k) + G(:,i,j,k))/(yf(j+1)-yf(j)) 
      enddo
      enddo
      enddo

      ! Source term
      do k=ks,ke
      do j=js,je
      do i=is,ie
         U(IPS,i,j,k) = U(IPS,i,j,k)*dexp(-alpha*Ccfl*dt1/dt)
      enddo
      enddo
      enddo


      return
      end subroutine UpdateConsv

      subroutine Output( flag_output, xf, xv, yf, yv, Q )
      real(8), intent(in) :: xf(:), xv(:), yf(:), yv(:), Q(:,:,:,:)
      integer::i,j,k
      logical, intent(in) :: flag_output ! false --> output per dtout, true --> force to output
      character(20),parameter::dirname="snap_B0.8_dc1"
      character(20),parameter::base="kh"
      character(20),parameter::suffix=".dat"
      character(100)::filename
      real(8), save :: tout = - dtout
      integer::nout = 0
      integer,parameter:: unitout=17
      integer,parameter:: unitbin=13
      integer,parameter:: gs=1

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
          write(unitbin,*) xv(i), yv(j), Q(IDN,i,j,k), Q(IV1,i,j,k), Q(IV2,i,j,k), Q(IV3,i,j,k), &
                          Q(IPR,i,j,k), Q(IB1,i,j,k), Q(IB2,i,j,k), Q(IB3,i,j,k),Q(ISC,i,j,k), Q(IPS,i,j,k) 
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
         error = dabs( ( Q(IB1,i,j,k) + 2.0d0*Q(IB2,i,j,k) )/sqrt(5.0d0) - 5.0d0/dsqrt(4.0d0*dacos(-1.0d0)))

         errorBpara = errorBpara + error
      enddo
      enddo
      enddo
      errorBpara = errorBpara/dble((ie-is+1)*(je-js+1)*(ke-ks+1))
      
      return
      end function

      real(8) function divergenceB(xf, yf, Q)
      implicit none
      real(8), intent(in) :: xf(:), yf(:), Q(:,:,:,:)
      integer::i,j,k
      real(8) :: error

      do k=ks,ke
      do j=js,je
      do i=is,ie
         error = dabs( ( Q(IB1,i+1,j,k) - Q(IB1,i-1,j,k) )/(xf(i+1)-xf(i-1))  &
                    + ( Q(IB2,i,j+1,k) - Q(IB2,i,j-1,k) )/(yf(j+1)-yf(j-1)) )/ &
                    dsqrt( Q(IB1,i,j,k)**2 + Q(IB2,i,j,k)**2 )*min(xf(i+1) - xf(i),yf(j+1)-yf(j))

         divergenceB = divergenceB + error
      enddo
      enddo
      enddo
      divergenceB = divergenceB/dble((ie-is+1)*(je-js+1)*(ke-ks+1))
      
      return
      end function


      subroutine Analysis(xv,yv,Q,phys_evo)
      implicit none
      real(8), intent(in) :: xv(:), yv(:), Q(:,:,:,:)
      real(8), intent(out) :: phys_evo(:)
      integer::i,j,k
      real(8) :: dvy, er_divB

      dvy = 0.0d0
      er_divB = 0.0d0
      do k=ks,ke
      do j=js,je
      do i=is,ie
           dvy = dvy + Q(IV2,i,j,k)**2
           er_divB = er_divB + ( Q(IB1,i+1,j,k) - Q(IB1,i-1,j,k) + Q(IB2,i,j+1,k) - Q(IB2,i,j-1,k) )**2 &
                       /( Q(IB1,i,j,k)**2 + Q(IB2,i,j,k)**2 )
      enddo
      enddo
      enddo
      phys_evo(1) = sqrt(dvy/dble(nx*ny))
      phys_evo(2) = sqrt(er_divB/dble(nx*ny))
      
      return
      end subroutine


end program main
