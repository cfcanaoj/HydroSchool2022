module commons
implicit none
integer::ntime                        ! counter of the timestep
integer,parameter::ntimemax=200000     ! the maximum timesteps
real(8)::time,dt                    ! time, timewidth
data time / 0.0d0 /
real(8),parameter:: timemax=0.5d0   
real(8),parameter:: dtout=5.0d-3

integer,parameter::nx=128        ! the number of grids in the simulation box
integer,parameter::ny=128        ! the number of grids in the simulation box
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
real(8),parameter::x2min=-0.5d0,x2max=0.5d0
real(8),parameter::x3min=-0.5d0,x3max=0.5d0

integer, parameter :: IDN = 1
integer, parameter :: IM1 = 2
integer, parameter :: IM2 = 3
integer, parameter :: IM3 = 4
integer, parameter :: IPR = 5
integer, parameter :: NVAR = 5

integer, parameter :: IB1 = NVAR+1
integer, parameter :: IB2 = NVAR+2
integer, parameter :: IB3 = NVAR+3
integer, parameter :: NFLX = NVAR+3

integer, parameter :: IBperp1 = NVAR+1
integer, parameter :: IBperp2 = NVAR+2
integer, parameter :: NFLX1D = NVAR+2

integer, parameter :: IV1 = 2
integer, parameter :: IV2 = 3
integer, parameter :: IV3 = 4
integer, parameter :: IEN = 5

end module commons
     
      module eosmod
      implicit none
! adiabatic
      real(8),parameter::gam=5.0d0/3.0d0 !! adiabatic index
!      real(8),parameter::gam=1.4d0 !! adiabatic index
! isothermal
!      real(8)::csiso  !! isothemal sound speed
      end module eosmod

      program main
      use commons
      implicit none

      real(8),dimension(in)::x1a,x1b
      real(8),dimension(jn)::x2a,x2b
      real(8),dimension(kn)::x3a,x3b
      real(8),dimension(in,jn,kn,NVAR) :: U
      real(8),dimension(in,jn,kn,NVAR) :: W
      real(8),dimension(in,jn,kn,3) :: B
      real(8),dimension(in,jn,kn,NFLX) :: flux1
      real(8),dimension(in,jn,kn,NFLX) :: flux2
      real(8),dimension(in,jn,kn,NFLX) :: flux3

!      write(6,*) "setup grids and initial condition"
      call GenerateGrid(x1a, x1b, x2a, x2b, x3a, x3b)
      call GenerateProblem(x1b, x2b, x3b, W, B)
      call ConsvVariable(W, B, U)
      call BoundaryCondition( W, B )
      call Output( x1a, x1b, x2a, x2b, W, B )

! main loop
      mloop: do ntime=1,ntimemax
         call TimestepControl(x1a, x2a, x3a, W, B)
         if( time + dt > timemax ) dt = timemax - time
         call BoundaryCondition( W, B )
         call NumericalFlux( W, B, flux1, flux2, flux3 )
         call UpdateConsv( x1a, x2a, x3a, flux1, flux2, flux3, U, B )
         call PrimVariable( U, B, W )
         time=time+dt
         call Output( x1a, x1b, x2a, x2b, W, B)


         if(time >= timemax) exit mloop
      enddo mloop
      call Output( x1a, x1b, x2a, x2b, W, B)

!      write(6,*) "program has been finished"
contains

      subroutine GenerateGrid(x1a, x1b, x2a, x2b, x3a, x3b)
      use commons
      use eosmod
      implicit none
      real(8), intent(out) :: x1a(:), x1b(:)
      real(8), intent(out) :: x2a(:), x2b(:)
      real(8), intent(out) :: x3a(:), x3b(:)
      real(8) :: dx,dy
      integer::i,j


      dx=(x1max-x1min)/dble(nx)
      do i=1,in
         x1a(i) = dx*(i-(mgn+1))+x1min
      enddo
      do i=1,in-1
         x1b(i) = 0.5d0*(x1a(i+1)+x1a(i))
      enddo

      dy=(x2max-x2min)/dble(ny)
      do j=1,jn
         x2a(j) = dx*(j-(mgn+1))+x1min
      enddo
      do j=1,jn-1
         x2b(j) = 0.5d0*(x2a(j+1)+x2a(j))
      enddo

      return
      end subroutine GenerateGrid

      subroutine GenerateProblem(x1b, x2b, x3b, W, B)
      use commons
      use eosmod
      implicit none
      integer::i, j, k
      real(8), intent(in ) :: x1b(:), x2b(:), x3b(:)
      real(8), intent(out) :: W(:,:,:,:), B(:,:,:,:)
      real(8) :: rho1,rho2,Lsm,u1,u2

      real(8)::pi, B0
      pi=acos(-1.0d0)

      B0 = 1.0d0/gam
      

      do k=ks,ke
      do j=js,je
      do i=is,ie

         W(i,j,k,IDN) = 1.0d0
         W(i,j,k,IPR) = 1.0d0/gam !0.95d0
         W(i,j,k,IV1) = - dsin(2.0d0*pi*x2b(j))
         W(i,j,k,IV2) =   dsin(2.0d0*pi*x1b(i))
         W(i,j,k,IV3) = 0.0d0
         B(i,j,k,1) = -B0*dsin(2.0d0*pi*x2b(j))
         B(i,j,k,2) =  B0*dsin(4.0d0*pi*x1b(i))
         B(i,j,k,3) =  0.0d0
      enddo
      enddo
      enddo

      
!      call BoundaryCondition

      return
      end subroutine GenerateProblem

      subroutine BoundaryCondition(W,B)
      use commons
      implicit none
      real(8), intent(inout) :: W(:,:,:,:), B(:,:,:,:)
      integer::i,j,k

      do k=ks,ke
      do j=1,jn-1
      do i=1,mgn
!          W(is-i,j,k,IDN)  = W(is-1+i,j,k,IDN)
!          W(is-i,j,k,IV1)  = W(is-1+i,j,k,IV1)
!          W(is-i,j,k,IV2)  = W(is-1+i,j,k,IV2)
!          W(is-i,j,k,IV3)  = W(is-1+i,j,k,IV3)
!          W(is-i,j,k,IPR)  = W(is-1+i,j,k,IPR)
!          B(is-i,j,k,1)  = B(is-1+i,j,k,1)
!          B(is-i,j,k,2)  = B(is-1+i,j,k,2)
!          B(is-i,j,k,3)  = B(is-1+i,j,k,3)

          W(is-i,j,k,IDN)  = W(ie+1-i,j,k,IDN)
          W(is-i,j,k,IV1)  = W(ie+1-i,j,k,IV1)
          W(is-i,j,k,IV2)  = W(ie+1-i,j,k,IV2)
          W(is-i,j,k,IV3)  = W(ie+1-i,j,k,IV3)
         W(is-i,j,k,IPR)  = W(ie+1-i,j,k,IPR)
          B(is-i,j,k,1)  = B(ie+1-i,j,k,1)
          B(is-i,j,k,2)  = B(ie+1-i,j,k,2)
          B(is-i,j,k,3)  = B(ie+1-i,j,k,3)
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=1,jn-1
      do i=1,mgn
!          W(ie+i,j,k,IDN) = W(ie-i+1,j,k,IDN)
!          W(ie+i,j,k,IV1) = W(ie-i+1,j,k,IV1)
!          W(ie+i,j,k,IV2) = W(ie-i+1,j,k,IV2)
!          W(ie+i,j,k,IV3) = W(ie-i+1,j,k,IV3)
!          W(ie+i,j,k,IPR) = W(ie-i+1,j,k,IPR)
!          B(ie+i,j,k,1) = B(ie-i+1,j,k,1)
!          B(ie+i,j,k,2) = B(ie-i+1,j,k,2)
!          B(ie+i,j,k,3) = B(ie-i+1,j,k,3)

          W(ie+i,j,k,IDN) = W(is+i-1,j,k,IDN)
          W(ie+i,j,k,IV1) = W(is+i-1,j,k,IV1)
          W(ie+i,j,k,IV2) = W(is+i-1,j,k,IV2)
          W(ie+i,j,k,IV3) = W(is+i-1,j,k,IV3)
          W(ie+i,j,k,IPR) = W(is+i-1,j,k,IPR)
          B(ie+i,j,k,1) = B(is+i-1,j,k,1)
          B(ie+i,j,k,2) = B(is+i-1,j,k,2)
          B(ie+i,j,k,3) = B(is+i-1,j,k,3)
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=1,mgn
      do i=1,in-1
!          W(i,js-j,k,IDN)  = W(i,js-1+j,k,IDN)
!          W(i,js-j,k,IV1)  = W(i,js-1+j,k,IV1)
!          W(i,js-j,k,IV2)  = W(i,js-1+j,k,IV2)
!          W(i,js-j,k,IV3)  = W(i,js-1+j,k,IV3)
!          W(i,js-j,k,IPR)  = W(i,js-1+j,k,IPR)
!          B(i,js-j,k,1)  = B(i,js-1+j,k,1)
!          B(i,js-j,k,2)  = B(i,js-1+j,k,2)
!          B(i,js-j,k,3)  = B(i,js-1+j,k,3)

          W(i,js-j,k,IDN)  = W(i,je+1-j,k,IDN)
          W(i,js-j,k,IV1)  = W(i,je+1-j,k,IV1)
          W(i,js-j,k,IV2)  = W(i,je+1-j,k,IV2)
          W(i,js-j,k,IV3)  = W(i,je+1-j,k,IV3)
          W(i,js-j,k,IPR)  = W(i,je+1-j,k,IPR)
          B(i,js-j,k,1)  = B(i,je+1-j,k,1)
          B(i,js-j,k,2)  = B(i,je+1-j,k,2)
          B(i,js-j,k,3)  = B(i,je+1-j,k,3)
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=1,mgn
      do i=1,in-1
!          W(i,je+j,k,IDN)  = W(i,je-j+1,k,IDN)
!          W(i,je+j,k,IV1)  = W(i,je-j+1,k,IV1)
!          W(i,je+j,k,IV2)  = W(i,je-j+1,k,IV2)
!          W(i,je+j,k,IV3)  = W(i,je-j+1,k,IV3)
!          W(i,je+j,k,IPR)  = W(i,je-j+1,k,IPR)
!          B(i,je+j,k,1)  = B(i,je-j+1,k,1)
!          B(i,je+j,k,2)  = B(i,je-j+1,k,2)
!          B(i,je+j,k,3)  = B(i,je-j+1,k,3)

          W(i,je+j,k,IDN)  = W(i,js+j-1,k,IDN)
          W(i,je+j,k,IV1)  = W(i,js+j-1,k,IV1)
          W(i,je+j,k,IV2)  = W(i,js+j-1,k,IV2)
          W(i,je+j,k,IV3)  = W(i,js+j-1,k,IV3)
          W(i,je+j,k,IPR)  = W(i,js+j-1,k,IPR)
          B(i,je+j,k,1)  = B(i,js+j-1,k,1)
          B(i,je+j,k,2)  = B(i,js+j-1,k,2)
          B(i,je+j,k,3)  = B(i,js+j-1,k,3)
      enddo
      enddo
      enddo



      return
      end subroutine BoundaryCondition
!
      subroutine ConsvVariable(W, B, U)
      use commons
      use eosmod
      implicit none
      real(8), intent(in) :: W(:,:,:,:), B(:,:,:,:)
      real(8), intent(out) :: U(:,:,:,:)
      integer::i,j,k

      do k=ks,ke
      do j=js,je
      do i=is,ie
          U(i,j,k,IDN) = W(i,j,k,IDN)
          U(i,j,k,IM1) = W(i,j,k,IDN)*W(i,j,k,IV1)
          U(i,j,k,IM2) = W(i,j,k,IDN)*W(i,j,k,IV2)
          U(i,j,k,IM3) = W(i,j,k,IDN)*W(i,j,k,IV3)
          U(i,j,k,IEN) = 0.5d0*W(i,j,k,IDN)*( W(i,j,k,IV1)**2 + W(i,j,k,IV2)**2 + W(i,j,k,IV3)**2 ) &
                       + 0.5d0*( B(i,j,k,1)**2 + B(i,j,k,2)**2 + B(i,j,k,3)**2 ) &
                       + W(i,j,k,IPR)/(gam - 1.0d0)
      enddo
      enddo
      enddo
      
      return
      end subroutine Consvvariable

      subroutine PrimVariable( U, B, W )
      use commons
      use eosmod
      implicit none
      real(8), intent(in) :: U(:,:,:,:), B(:,:,:,:)
      real(8), intent(out) :: W(:,:,:,:)
      integer::i,j,k
      real(8) :: inv_d;

      do k=ks,ke
      do j=js,je
      do i=is,ie
           W(i,j,k,IDN) = U(i,j,k,IDN)
           inv_d = 1.0d0/U(i,j,k,IDN)
           W(i,j,k,IV1) = U(i,j,k,IM1)*inv_d
           W(i,j,k,IV2) = U(i,j,k,IM2)*inv_d
           W(i,j,k,IV3) = U(i,j,k,IM3)*inv_d
           W(i,j,k,IPR) = ( U(i,j,k,IEN) &
                        - 0.5d0*(U(i,j,k,IM1)**2 + U(i,j,k,IM2)**2 + U(i,j,k,IM3)**2)*inv_d  &
                        - 0.5d0*(B(i,j,k,1)**2 + B(i,j,k,2)**2 + B(i,j,k,3)**2) )*(gam-1.0d0)
      enddo
      enddo
      enddo

      return
      end subroutine PrimVariable

      subroutine TimestepControl(x1a, x2a, x3a, W, B)
      use commons
      use eosmod
      implicit none
      real(8), intent(in) :: x1a(:), x2a(:), x3a(:), W(:,:,:,:), B(:,:,:,:)
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
         cf = dsqrt( (gam*W(i,j,k,IPR) + B(i,j,k,1)**2 + B(i,j,k,2)**2 + B(i,j,k,3)**2)/W(i,j,k,IDN))
         dtl1 =(x1a(i+1)-x1a(i))/(abs(W(i,j,k,IV1)) + cf)
         dtl2 =(x2a(j+1)-x2a(j))/(abs(W(i,j,k,IV2)) + cf)
!         dtl3 =(x3a(j+1)-x3a(j))/(abs(W(i,j,k,IV3)) + cf)
         dtlocal = min(dtl1,dtl2)
         if(dtlocal .lt. dtmin) dtmin = dtlocal
      enddo
      enddo
      enddo

      dt = 0.2d0 * dtmin
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
      integer :: i

      do i=1,n
         if(dvp(i)*dvm(i) .gt. 0.0d0)then
            dv(i) =2.0d0*dvp(i)*dvm(i)/(dvp(i)+dvm(i))
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
!     Input: W: primitive variables at the cell center
!
!     Input: B: magnetic fields
!
!     Output: flux : the numerical flux estimated at the cell boundary
!---------------------------------------------------------------------
      subroutine NumericalFlux( W, B, flux1, flux2, flux3 )
      use commons !, only: is, ie, in
      implicit none
      integer::i,j,k
      real(8), intent(in) :: W(:,:,:,:), B(:,:,:,:)
      real(8), intent(out) :: flux1(:,:,:,:)
      real(8), intent(out) :: flux2(:,:,:,:)
      real(8), intent(out) :: flux3(:,:,:,:)
      real(8),dimension(in,jn,kn,NFLX1D):: Wl,Wr
      real(8),dimension(NFLX1D):: flx
      real(8) :: dWm(NFLX1D), dWp(NFLX1D), dWmon(NFLX1D)
      real(8) :: ddmon, dvmon, dpmon

      ! numerical flux in the x direction
      do k=ks,ke
      do j=js,je
      do i=is-1,ie+1
         dWp(1:NVAR) = W(i+1,j,k,1:NVAR) - W(i  ,j,k,1:NVAR)
         dWm(1:NVAR) = W(i  ,j,k,1:NVAR) - W(i-1,j,k,1:NVAR)

         dWp(NVAR+1) = B(i+1,j,k,2) - B(i,  j,k,2)
         dWm(NVAR+1) = B(i  ,j,k,2) - B(i-1,j,k,2)

         dWp(NVAR+2) = B(i+1,j,k,3) - B(i,  j,k,3)
         dWm(NVAR+2) = B(i  ,j,k,3) - B(i-1,j,k,3)

         call vanLeer(NFLX1D, dWp, dWm, dWmon)

         ! Wl(i,j,k) --> W_(i-1/2,j,k)
         ! Wr(i,j,k) --> W_(i-1/2,j,k)
         Wl(i+1,j,k,1:NVAR) = W(i,j,k,1:NVAR) + 0.5d0*dWmon(1:NVAR)*1.0d0
         Wr(i  ,j,k,1:NVAR) = W(i,j,k,1:NVAR) - 0.5d0*dWmon(1:NVAR)*1.0d0

         Wl(i+1,j,k,NVAR+1) = B(i,j,k,2) + 0.5d0*dWmon(NVAR+1)*1.0d0
         Wr(i  ,j,k,NVAR+1) = B(i,j,k,2) - 0.5d0*dWmon(NVAR+1)*1.0d0

         Wl(i+1,j,k,NVAR+2) = B(i,j,k,3) + 0.5d0*dWmon(NVAR+2)*1.0d0
         Wr(i  ,j,k,NVAR+2) = B(i,j,k,3) - 0.5d0*dWmon(NVAR+2)*1.0d0
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=js,je
      do i=is,ie+1
        call HLL(1,Wl(i,j,k,:),Wr(i,j,k,:),0.5d0*(B(i-1,j,k,1) + B(i,j,k,1)),flx)

         flux1(i,j,k,IDN)  = flx(IDN)
         flux1(i,j,k,IM1)  = flx(IM1)
         flux1(i,j,k,IM2)  = flx(IM2)
         flux1(i,j,k,IM3)  = flx(IM3)
         flux1(i,j,k,IEN)  = flx(IEN)
         flux1(i,j,k,IB1)  = 0.0
         flux1(i,j,k,IB2)  = flx(NVAR+1)
         flux1(i,j,k,IB3)  = flx(NVAR+2)
      enddo
      enddo
      enddo

      ! numerical flux in the y direction
      do k=ks,ke
      do j=js-1,je+1
      do i=is,ie
         dWp(1:NVAR) = W(i,j+1,k,1:NVAR) - W(i,j  ,k,1:NVAR)
         dWm(1:NVAR) = W(i,j  ,k,1:NVAR) - W(i,j-1,k,1:NVAR)

         dWp(NVAR+1) = B(i,j+1,k,3) - B(i,j  ,k,3)
         dWm(NVAR+1) = B(i,j  ,k,3) - B(i,j-1,k,3)

         dWp(NVAR+2) = B(i,j+1,k,1) - B(i,j  ,k,1)
         dWm(NVAR+2) = B(i,j,  k,1) - B(i,j-1,k,1)

         call vanLeer(NFLX1D, dWp, dWm, dWmon)

         ! Wl(i,j,k) --> W_(i-1/2,j,k)
         ! Wr(i,j,k) --> W_(i-1/2,j,k)
         Wl(i,j+1,k,1:NVAR) = W(i,j,k,1:NVAR) + 0.5d0*dWmon(1:NVAR)*1.0d0
         Wr(i,j  ,k,1:NVAR) = W(i,j,k,1:NVAR) - 0.5d0*dWmon(1:NVAR)*1.0d0

         Wl(i,j+1,k,NVAR+1) = B(i,j,k,3) + 0.5d0*dWmon(NVAR+1)*1.0d0
         Wr(i,j  ,k,NVAR+1) = B(i,j,k,3) - 0.5d0*dWmon(NVAR+1)*1.0d0

         Wl(i,j+1,k,NVAR+2) = B(i,j,k,1) + 0.5d0*dWmon(NVAR+2)*1.0d0
         Wr(i,j  ,k,NVAR+2) = B(i,j,k,1) - 0.5d0*dWmon(NVAR+2)*1.0d0

      enddo
      enddo
      enddo

      do k=ks,ke
      do j=js,je+1
      do i=is,ie
         call HLL(2,Wl(i,j,k,:),Wr(i,j,k,:),0.5d0*(B(i-1,j,k,2) + B(i,j,k,2)),flx)

         flux2(i,j,k,IDN)  = flx(IDN)
         flux2(i,j,k,IM1)  = flx(IM1)
         flux2(i,j,k,IM2)  = flx(IM2)
         flux2(i,j,k,IM3)  = flx(IM3)
         flux2(i,j,k,IEN)  = flx(IEN)
         flux2(i,j,k,IB1)  = flx(NVAR+2)
         flux2(i,j,k,IB2)  = 0.0d0
         flux2(i,j,k,IB3)  = flx(NVAR+1)
!         if( i == is+10 ) print*,x2b(j),flx(IM3)
      enddo
      enddo
      enddo


      return
      end subroutine Numericalflux

!---------------------------------------------------------------------
!     HLL Riemann Solver
!---------------------------------------------------------------------
!     solve the HLL Riemann solver 
!
!     Input: Wl, Wr: primitive variables containing the perpendicular B fields 
!                    at the left and right states
!            1D array (IDN, IV1, IV2, IV3, IPR, IBperp1, IBperp2)
!                                 |
!                                 |
!                           Wl    |    Wr
!                                 |
!                                -->
!                                flx
!
!     Input: b1    : magnetic field perpendicular to the initial discontinuity
!
!     Output: flx  : flux estimated at the initial discontinuity
!            index: (IDN, IV1, IV2, IV3, IPR, IBperp1, IBperp2)
!---------------------------------------------------------------------
      subroutine HLL(idir,Wl,Wr,b1,flx)
      use commons !, only : is, ie, NVAR
      use eosmod
      implicit none
      integer, intent(in) :: idir
      real(8),intent(in)::Wl(:), Wr(:)
      real(8),intent(in) :: b1
      real(8),intent(out) :: flx(:)
      integer :: IVpara, IVperp1, IVperp2
      real(8):: Ul(NFLX1D), Ur(NFLX1D)
      real(8):: Fl(NFLX1D), Fr(NFLX1D)
      real(8):: Ust(NFLX1D)
      real(8):: Fst(NFLX1D)
      real(8):: cfl,cfr
      real(8):: sl, sr
      real(8):: pbl, pbr, ptotl, ptotr
      integer :: i, n

      if( idir == 1 ) then
           IVpara = IV1
           IVperp1 = IV2
           IVperp2 = IV3
      else if (idir == 2 ) then
           IVpara = IV2
           IVperp1 = IV3
           IVperp2 = IV1
      else 
           IVpara = IV3
           IVperp1 = IV1
           IVperp2 = IV2
      endif

          
          pbl = 0.5d0*(b1**2 + Wl(IBperp1)**2 + Wl(IBperp2)**2)
          pbr = 0.5d0*(b1**2 + Wr(IBperp1)**2 + Wr(IBperp2)**2)
          ptotl = Wl(IPR) + pbl
          ptotr = Wr(IPR) + pbr

          ! conserved variables in the left and right states
          Ul(IDN) = Wl(IDN)
          Ul(IM1) = Wl(IDN)*Wl(IVpara)
          Ul(IM2) = Wl(IDN)*Wl(IVperp1)
          Ul(IM3) = Wl(IDN)*Wl(IVperp2)
          Ul(IEN) = 0.5d0*Wl(IDN)*( Wl(IVpara)**2 + Wl(IVperp1)**2 + Wl(IVperp2)**2) & 
                  + pbl + Wl(IPR)/(gam - 1.0d0)
          Ul(IBperp1) = Wl(IBperp1)
          Ul(IBperp2) = Wl(IBperp2)

          Ur(IDN) = Wr(IDN)
          Ur(IM1) = Wr(IDN)*Wr(IVpara)
          Ur(IM2) = Wr(IDN)*Wr(IVperp1)
          Ur(IM3) = Wr(IDN)*Wr(IVperp2)
          Ur(IEN) = 0.5d0*Wr(IDN)*( Wr(IVpara)**2 + Wr(IVperp1)**2 + Wr(IVperp2)**2) & 
                  + pbr + Wr(IPR)/(gam - 1.0d0)
          Ur(IBperp1) = Wr(IBperp1)
          Ur(IBperp2) = Wr(IBperp2)

          ! flux in the left and right states
          Fl(IDN) = Ul(IM1)
          Fl(IM1) = Wl(IDN)*Wl(IVpara)**2 + ptotl - b1**2
          Fl(IM2) = Wl(IDN)*Wl(IVpara)*Wl(IVperp1) - b1*Wl(IBperp1)
          Fl(IM3) = Wl(IDN)*Wl(IVpara)*Wl(IVperp2) - b1*Wl(IBperp2)
          Fl(IEN) = ( Ul(IEN) + ptotl )*Wl(IVpara) &
                  - b1*( b1*Wl(IVpara) + Wl(IBperp1)*Wl(IVperp1) + Wl(IBperp2)*Wl(IVperp2) )
          Fl(IBperp1) = Wl(IBperp1)*Wl(IVpara) - b1*Wl(IVperp1)
          Fl(IBperp2) = Wl(IBperp2)*Wl(IVpara) - b1*Wl(IVperp2)



          Fr(IDN) = Ur(IM1)
          Fr(IM1) = Wr(IDN)*Wr(IVpara)**2 + ptotr - b1**2
          Fr(IM2) = Wr(IDN)*Wr(IVpara)*Wr(IVperp1) - b1*Wr(IBperp1)
          Fr(IM3) = Wr(IDN)*Wr(IVpara)*Wr(IVperp2) - b1*Wr(IBperp2)
          Fr(IEN) = ( Ur(IEN) + ptotr )*Wr(IVpara) &
                  - b1*( b1*Wr(IVpara) + Wr(IBperp1)*Wr(IVperp1) + Wr(IBperp2)*Wr(IVperp2) )
          Fr(IBperp1) = Wr(IBperp1)*Wr(IVpara) - b1*Wr(IVperp1)
          Fr(IBperp2) = Wr(IBperp2)*Wr(IVpara) - b1*Wr(IVperp2)

!         print*, Wr(IDN)*Wr(IVpara)*Wr(IVperp2) - b1*Wr(IBperp2)
!         print*, Wl(IDN)*Wl(IVpara)*Wl(IVperp2) - b1*Wl(IBperp2)


          cfl = dsqrt( (gam*Wl(IPR) + Wl(IBperp1)**2 + Wl(IBperp2)**2 + b1**2)/Wl(IDN))
          cfr = dsqrt( (gam*Wr(IPR) + Wr(IBperp1)**2 + Wr(IBperp2)**2 + b1**2)/Wr(IDN))

!          sl = min(Wl(IV1),Wr(IV1)) - max(cfl,cfr)
!          sr = max(Wl(IV1),Wr(IV1)) + max(cfl,cfr)
          sl = min(Wl(IVpara) - cfl,Wr(IVpara) - cfr)
          sr = max(Wl(IVpara) + cfl,Wr(IVpara) + cfr)
!          Fst(:)  = (sr*Fl(:) - sl*Fr(:) + sl*sr*( Ur(:) - Ul(:) ))/(sr - sl)
!          Ust(:) = ( sr*Ur(:) - sl*Ul(:) - Fr(:) + Fl(:) )/(sr - sl)

          if( sl > 0.0d0 ) then
               flx(IDN) = Fl(IDN)
               flx(IVpara) = Fl(IM1)
               flx(IVperp1) = Fl(IM2)
               flx(IVperp2) = Fl(IM3)
               flx(IEN) = Fl(IEN)
               flx(IBperp1) = Fl(IBperp1)
               flx(IBperp2) = Fl(IBperp2)
          else if (sr <= 0.0d0 ) then
               flx(IDN) = Fr(IDN)
               flx(IVpara) = Fr(IM1)
               flx(IVperp1) = Fr(IM2)
               flx(IVperp2) = Fr(IM3)
               flx(IEN) = Fr(IEN)
               flx(IBperp1) = Fr(IBperp1)
               flx(IBperp2) = Fr(IBperp2)
          else 
               flx(IDN)  = (sr*Fl(IDN) - sl*Fr(IDN) + sl*sr*( Ur(IDN) - Ul(IDN) ))/(sr - sl)
               flx(IVpara)  = (sr*Fl(IM1) - sl*Fr(IM1) + sl*sr*( Ur(IM1) - Ul(IM1) ))/(sr - sl)
               flx(IVperp1)  = (sr*Fl(IM2) - sl*Fr(IM2) + sl*sr*( Ur(IM2) - Ul(IM2) ))/(sr - sl)
               flx(IVperp2)  = (sr*Fl(IM3) - sl*Fr(IM3) + sl*sr*( Ur(IM3) - Ul(IM3) ))/(sr - sl)
               flx(IEN)  = (sr*Fl(IEN) - sl*Fr(IEN) + sl*sr*( Ur(IEN) - Ul(IEN) ))/(sr - sl)
               flx(IBperp1)  = (sr*Fl(IBperp1) - sl*Fr(IBperp1) + sl*sr*( Ur(IBperp1) - Ul(IBperp1) ))/(sr - sl)
               flx(IBperp2)  = (sr*Fl(IBperp2) - sl*Fr(IBperp2) + sl*sr*( Ur(IBperp2) - Ul(IBperp2) ))/(sr - sl)
          endif


!         do i=1,NFLX1D 
!         if( flx(i) .ne. flx(i) ) then 
!             print*,(sr*Fl(:) - sl*Fr(:) + sl*sr*( Ur(:) - Ul(:) ))/(sr - sl)
!         stop
!         endif
!         enddo
!
      return
      end subroutine HLL

!      subroutine HLLC(dl,vl,pl,dr,vr,pr)
!!=====================================================================
!!
!! HLLC Scheme
!!
!! Purpose
!! Calculation of Numerical Flux by HLLC method
!!
!! Reference
!!  Toro EF, Spruce M, Speares W. (1992,1994)
!!
!! Input
!! Output
!!=====================================================================
!      use fluxmod, only: dflux, mvflux, etflux 
!
!      implicit none
!      real(8),dimension(in),intent(in) :: dl, dr, vl, vr, pl, pr
!
!!----- U -----
!! qql :: left state
!! qqr :: right state
!      real(8) :: rol,vxl,vyl,vzl,ptl,eel
!      real(8) :: ror,vxr,vyr,vzr,ptr,eer
!      real(8) :: rxl,ryl,rzl
!      real(8) :: rxr,ryr,rzr
!      real(8) :: ptst
!
!!----- U* ----
!! qqlst ::  left state
!! qqrst :: right state
!      real(8) :: rolst,vxlst,vylst,vzlst,eelst
!      real(8) :: rorst,vxrst,vyrst,vzrst,eerst
!      real(8) :: rxlst,rylst,rzlst
!      real(8) :: rxrst,ryrst,rzrst
!
!!----- flux ---
!! fqql ::  left physical flux
!! fqqr :: right physical flux
!      real(8) :: frol,frxl,fryl,frzl,feel
!      real(8) :: fror,frxr,fryr,frzr,feer
!
!!----- wave speed ---
!! sl ::  left-going fastest signal velocity
!! sr :: right-going fastest signal velocity
!! sm :: contact discontinuity velocity
!! slst ::  left-going alfven velocity
!! srst :: right-going alfven velocity
!      real(8) :: sm,sl,sr
!
!! cfl :: left-state Fast wave velocity
!! cfr :: right-sate Fast wave velocity
!      real(8) :: cfl,cfr
!
!!--------------------
!! temporary variables
!      real(8) :: sdl,sdr,sdml,sdmr,isdml,isdmr,rosdl,rosdr
!      real(8) :: temp
!  
!! no if
!      real(8) :: sign1,maxs1,mins1
!      real(8) :: msl,msr
!
!!----- Step 0. ----------------------------------------------------------|
!
!!---- Left state
!        
!        rol = leftst(mudn)
!        eel = leftst(muet)
!        rxl = leftst(muvu)
!        ryl = leftst(muvv)
!        rzl = leftst(muvw)
!        vxl = leftst(muvu)/leftst(mudn)
!        vyl = leftst(muvv)/leftst(mudn)
!        vzl = leftst(muvw)/leftst(mudn)
!        ptl = leftst(mpre)
!
!!---- Right state
!        
!        ror = rigtst(mudn)
!        eer = rigtst(muet)
!        rxr = rigtst(muvu)
!        ryr = rigtst(muvv)
!        rzr = rigtst(muvw)
!        vxr = rigtst(muvu)/rigtst(mudn)
!        vyr = rigtst(muvv)/rigtst(mudn)
!        vzr = rigtst(muvw)/rigtst(mudn)
!        ptr = rigtst(mpre)
!!----- Step 1. ----------------------------------------------------------|
!! Compute wave left & right wave speed
!!
!         
!        cfl = leftst(mcsp)
!        cfr = rigtst(mcsp)
!
!        sl = min(vxl,vxr)-max(cfl,cfr) ! note sl is negative
!        sr = max(vxl,vxr)+max(cfl,cfr)
!!----- Step 2. ----------------------------------------------------------|
!! compute L/R fluxs
!!
!! Left value
!        frol = leftst(mfdn)
!        feel = leftst(mfet)
!        frxl = leftst(mfvu)
!        fryl = leftst(mfvv)
!        frzl = leftst(mfvw)
!
!! Right value
!! Left value
!        fror = rigtst(mfdn)
!        feer = rigtst(mfet)
!        frxr = rigtst(mfvu)
!        fryr = rigtst(mfvv) 
!        frzr = rigtst(mfvw)
!
!!----- Step 4. ----------------------------------------------------------|
!! compute middle and alfven wave
!!
!        sdl = sl - vxl
!        sdr = sr - vxr
!        rosdl = rol*sdl
!        rosdr = ror*sdr
!
!        temp = 1.0d0/(rosdr - rosdl)
!! Eq. 45
!        sm = (rosdr*vxr - rosdl*vxl - ptr + ptl)*temp
!           
!        sdml = sl - sm; isdml = 1.0d0/sdml
!        sdmr = sr - sm; isdmr = 1.0d0/sdmr
!        
!!----- Step 5. ----------------------------------------------------------|
!! compute intermediate states
!!
!! Eq. 49
!        ptst = (rosdr*ptl-rosdl*ptr+rosdl*rosdr*(vxr-vxl))*temp
!
!!----- Step 5A. ----------------------------------------------------------|
!! compute Ul*
!!
!
!        rolst = rol*sdl   *isdml
!        vxlst = sm
!        rxlst = rolst*vxlst
!           
!        vylst = vyl
!        rylst = rolst*vylst
!        vzlst = vzl
!        rzlst = rolst*vzlst
!
!        eelst =(sdl*eel - ptl*vxl + ptst*sm  )*isdml
!
!!----- Step 5B. ----------------------------------------------------------|
!! compute Ur*
!!
!
!        rorst   = rosdr   *isdmr
!        vxrst = sm
!        rxrst = rorst*vxrst
!        vyrst = vyr
!        ryrst = rorst*vyrst
!        vzrst = vzr
!        rzrst = rorst*vzrst
!           
!        eerst = (sdr*eer - ptr*vxr  + ptst*sm  )*isdmr
!              
!!----- Step 6. ----------------------------------------------------------|
!! compute flux
!        sign1 = sign(1.0d0,sm)    ! 1 for sm>0, -1 for sm<0
!        maxs1 =  max(0.0d0,sign1) ! 1 sm>0, 0 for sm<0
!        mins1 = -min(0.0d0,sign1) ! 0 sm>0,-1 for sm<0
!
!        msl   = min(sl  ,0.0d0)   ! 0 for sl > 0, sl for sl < 0
!        msr   = max(sr  ,0.0d0)   ! S_R > 0
!
!        nflux(mden) = (frol+msl*(rolst-rol))*maxs1 &
!     &               +(fror+msr*(rorst-ror))*mins1
!        nflux(meto) = (feel+msl*(eelst-eel))*maxs1 &
!     &               +(feer+msr*(eerst-eer))*mins1
!        nflux(mrvu) = (frxl+msl*(rxlst-rxl))*maxs1 &
!     &               +(frxr+msr*(rxrst-rxr))*mins1
!        nflux(mrvv) = (fryl+msl*(rylst-ryl))*maxs1 &
!     &               +(fryr+msr*(ryrst-ryr))*mins1
!        nflux(mrvw) = (frzl+msl*(rzlst-rzl))*maxs1 &
!     &               +(frzr+msr*(rzrst-rzr))*mins1
!
!      return
!      end subroutine HLLC

      subroutine UpdateConsv( x1a, x2a, x3a, flux1, flux2, flux3, U, B )
      use commons
      implicit none
      real(8), intent(in)  :: x1a(:), x2a(:), x3a(:)
      real(8), intent(in)  :: flux1(:,:,:,:), flux2(:,:,:,:), flux3(:,:,:,:)
      real(8), intent(out) :: U(:,:,:,:), B(:,:,:,:)
      integer::i,n,j,k

      do n=1,NVAR
      do k=ks,ke
      do j=js,je
      do i=is,ie
         U(i,j,k,n) = U(i,j,k,n) + dt*(- flux1(i+1,j,k,n) + flux1(i,j,k,n))/(x1a(i+1)-x1a(i)) &
                                 + dt*(- flux2(i,j+1,k,n) + flux2(i,j,k,n))/(x2a(j+1)-x2a(j)) 
      enddo
      enddo
      enddo
      enddo

      do n=1,3
      do k=ks,ke
      do j=js,je
      do i=is,ie
         B(i,j,k,n) = B(i,j,k,n) + dt*(- flux1(i+1,j,k,NVAR+n) + flux1(i,j,k,NVAR+n))/(x1a(i+1)-x1a(i))  &
                                 + dt*(- flux2(i,j+1,k,NVAR+n) + flux2(i,j,k,NVAR+n))/(x2a(j+1)-x2a(j)) 
      enddo
      enddo
      enddo
      enddo
      

      return
      end subroutine UpdateConsv

      subroutine Output( x1a, x1b, x2a, x2b, W, B )
      use commons
      implicit none
      real(8), intent(in) :: x1a(:), x1b(:), x2a(:), x2b(:), W(:,:,:,:), B(:,:,:,:)
      integer::i,j,k
      character(20),parameter::dirname="bindata_y/"
      character(40)::filename
      real(8),save::tout
      data tout / 0.0d0 /
      integer::nout
      data nout / 1 /
      integer,parameter:: unitout=17
      integer,parameter:: unitbin=13
      integer,parameter:: gs=1

      logical, save:: is_inited
      data is_inited /.false./

      if (.not. is_inited) then
         call makedirs(dirname)
         is_inited =.true.
      endif

!      print*, time, tout+dtout, time+1.0d-14 .lt. tout+dtout
      if(time + 1.0d-14.lt. tout+dtout) return

      write(filename,'(a3,i5.5,a4)')"bin",nout,".dat"
      filename = trim(dirname)//filename
!      open(unitbin,file=filename,status='replace',form='formatted') 
      open(unitbin,file=filename,form='formatted',action="write")
      write(unitbin,*) "# time = ",time
      do k=ks,ke
      do j=1,jn-1
      do i=1,in-1
          write(unitbin,*) x1b(i), x2b(j), W(i,j,k,IDN), W(i,j,k,IV1), W(i,j,k,IV2), W(i,j,k,IV3), W(i,j,k,IPR), B(i,j,k,1), &
          B(i,j,k,2), B(i,j,k,3)
!          write(*,*) x1b(i), d(i), v(i), p(i)
      enddo
      enddo
      enddo
      close(unitbin)
!      open(unitbin,file=filename,status='replace',form='binary') 
!      open(unitbin,file=filename,status='replace',form='unformatted') 
!      write(unitbin) x1out(:,:)
!      write(unitbin) hydout(:,:)
!      close(unitbin)
!
      write(6,*) "output:",nout,time

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
end program main
