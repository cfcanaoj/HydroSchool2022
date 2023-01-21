module commons
implicit none
integer::ntime                        ! counter of the timestep
integer,parameter::ntimemax=200000     ! the maximum timesteps
real(8)::time,dt                    ! time, timewidth
data time / 0.0d0 /
real(8),parameter:: timemax=100.0d0   
real(8),parameter:: dtout=0.1d0

integer,parameter::nx=128*2        ! the number of grids in the simulation box
integer,parameter::ny=128*2        ! the number of grids in the simulation box
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

real(8),parameter::Ccfl=0.2d0

integer, parameter :: IDN = 1
integer, parameter :: IM1 = 2
integer, parameter :: IM2 = 3
integer, parameter :: IM3 = 4
integer, parameter :: IPR = 5
integer, parameter :: IB1 = 6
integer, parameter :: IB2 = 7
integer, parameter :: IB3 = 8
integer, parameter :: IPS = 9
integer, parameter :: NVAR = 9
integer, parameter :: NFLX = 9

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
      real(8),dimension(in,jn,kn,NFLX) :: flux1
      real(8),dimension(in,jn,kn,NFLX) :: flux2
      real(8),dimension(in,jn,kn,NFLX) :: flux3

!      write(6,*) "setup grids and initial condition"
      call GenerateGrid(x1a, x1b, x2a, x2b, x3a, x3b)
      call GenerateProblem(x1b, x2b, x3b, W)
      call ConsvVariable(W, U)
      call BoundaryCondition( W)
      call Output( x1a, x1b, x2a, x2b, W )

      open(1,file="dvy_256.dat",action="write")
! main loop
      mloop: do !ntime=1,ntimemax
         call TimestepControl(x1a, x2a, x3a, W)
         if( time + dt > timemax ) dt = timemax - time
         call BoundaryCondition( W )

         write(1,*) time, dvy( x1a, x2a, W )

         call NumericalFlux( W, flux1, flux2, flux3 )
         call UpdateConsv( x1a, x2a, x3a, flux1, flux2, flux3, U )
         call PrimVariable( U, W )
         time=time+dt
         call Output( x1a, x1b, x2a, x2b, W)

         if(time >= timemax) exit mloop
      enddo mloop
      close(1)
      call Output( x1a, x1b, x2a, x2b, W)

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

      subroutine GenerateProblem(x1b, x2b, x3b, W )
      use commons
      use eosmod
      implicit none
      integer::i, j, k
      real(8), intent(in ) :: x1b(:), x2b(:), x3b(:)
      real(8), intent(out) :: W(:,:,:,:)
      real(8) :: rho1,rho2,pre,v1,v2,dv
      real(8)::pi, B0, fac, h
      pi=acos(-1.0d0)

      B0 = 0.0d0

      rho1 = 1.0d0
      rho2 = 2.0d0
      dv = 2.00d0
      v1 = dv*rho2/(rho1+rho2)
      v2 = -dv*rho1/(rho1+rho2)

      h = 0.2d0

      do k=ks,ke
      do j=js,je
      do i=is,ie
         W(i,j,k,IDN) = 1.0d0 + 1.0d0*( dtanh( (x2b(j)+0.25d0)/0.025d0 ) - tanh( (x2b(j)-0.25d0)/0.025d0) )
         W(i,j,k,IV1)  = dtanh( (x2b(j)+0.25d0)/0.025d0 ) - dtanh( (x2b(j) - 0.25d0)/0.025d0 ) - 1.0d0
         W(i,j,k,IV2)  = 0.01d0*dsin(2.0d0*pi*x1b(i))* &
             ( dexp( - (x2b(j) + 0.25d0)**2/0.1d0**2 ) +  &
               dexp( - (x2b(j) - 0.25d0)**2/0.1d0**2 ) )
         W(i,j,k,IPR) = 10.0d0

!         fac = 0.5d0*(tanh( (x2b(j) + 0.25d0)/0.025d0) - tanh( (x2b(j) - 0.25d0)/0.025d0))
!!
!         W(i,j,k,IDN) = (rho1 - rho2)*fac + rho2
!         W(i,j,k,IPR) = 10.0d0
!         W(i,j,k,IV1) = (v1 - v2)*fac + v2
!         W(i,j,k,IV2)  = 0.01d0*dsin(2.0d0*pi*x1b(i))* &
!             ( dexp( - (x2b(j) + 0.25d0)**2/0.1d0**2 ) +  &
!               dexp( - (x2b(j) - 0.25d0)**2/0.1d0**2 ) )
!         W(i,j,k,IV3) = 0.0d0
!         W(i,j,k,IB1) = -B0*dsin(2.0d0*pi*x2b(j))
!         W(i,j,k,IB2) =  B0*dsin(4.0d0*pi*x1b(i))
!         W(i,j,k,IB3) =  0.0d0
!         W(i,j,k,IPS) = 0.0d0
      enddo
      enddo
      enddo

!      do k=ks,ke
!      do j=js,je
!      do i=is,ie
!
!         if( x2b(j) < 0.0d0 ) then 
!             W(i,j,k,IDN) = 1.08d0
!             W(i,j,k,IV1) = 0.5d0
!             W(i,j,k,IV2) = 1.2d0 
!             W(i,j,k,IV3) = 0.01d0 ! 0.01d0
!             W(i,j,k,IPR) = 0.95d0 !0.95d0
!
!             B(i,j,k,1) = 2.0d0/dsqrt(4.0d0*3.141592653589793)
!             B(i,j,k,2) = 4.0d0/dsqrt(4.0d0*3.141592653589793)
!             B(i,j,k,3) = 3.6d0/dsqrt(4.0d0*3.141592653589793)
!
!        else 
!             W(i,j,k,IDN) = 1.0d0
!             W(i,j,k,IV1) = 0.0d0
!             W(i,j,k,IV2) = 0.0d0
!             W(i,j,k,IV3) = 0.0d0
!             W(i,j,k,IPR) = 1.0d0
!
!             B(i,j,k,1) = 2.0d0/dsqrt(4.0d0*3.141592653589793)
!             B(i,j,k,2) = 4.0d0/dsqrt(4.0d0*3.141592653589793)
!             B(i,j,k,3) = 4.0d0/dsqrt(4.0d0*3.141592653589793)
!         endif
!
!!         if( x1b(i) < 0.0d0 ) then 
!!             W(i,j,k,IDN) = 1.08d0
!!             W(i,j,k,IV1) = 1.2d0
!!             W(i,j,k,IV2) = 0.01d0
!!             W(i,j,k,IV3) = 0.5d0
!!             W(i,j,k,IPR) = 0.95d0
!!
!!             B(i,j,k,1) = 4.0d0/dsqrt(4.0d0*3.141592653589793)
!!             B(i,j,k,2) = 3.6d0/dsqrt(4.0d0*3.141592653589793)
!!             B(i,j,k,3) = 2.0d0/dsqrt(4.0d0*3.141592653589793)
!!
!!         else 
!!             W(i,j,k,IDN) = 1.0d0
!!             W(i,j,k,IV1) = 0.0d0
!!             W(i,j,k,IV2) = 0.0d0
!!             W(i,j,k,IV3) = 0.0d0
!!             W(i,j,k,IPR) = 1.0d0
!!
!!             B(i,j,k,1) = 4.0d0/dsqrt(4.0d0*3.141592653589793)
!!             B(i,j,k,2) = 4.0d0/dsqrt(4.0d0*3.141592653589793)
!!             B(i,j,k,3) = 2.0d0/dsqrt(4.0d0*3.141592653589793)
!!         endif
!!
!
!!         if      ( x2b(j) .gt. 0.25d0 )then
!!             W(i,j,k,IV1) =    u1 - (  u1-  u2)/2.0d0*exp(-( x2b(j)-0.25d0)/Lsm)
!!             W(i,j,k,IDN) =  rho1 - (rho1-rho2)/2.0d0*exp(-( x2b(j)-0.25d0)/Lsm)
!!         else if (x2b(j) .gt.  0.0d0 )then
!!             W(i,j,k,IV1) =    u2 + (  u1-  u2)/2.0d0*exp(-( 0.25d0-x2b(j))/Lsm)
!!             W(i,j,k,IDN) =  rho2 + (rho1-rho2)/2.0d0*exp(-( 0.25d0-x2b(j))/Lsm)
!!         else if (x2b(j) .gt. -0.25d0)then
!!             W(i,j,k,IV1) =    u2 + (  u1-  u2)/2.0d0*exp(-( x2b(j)+0.25d0)/Lsm)
!!             W(i,j,k,IDN) =  rho2 + (rho1-rho2)/2.0d0*exp(-( x2b(j)+0.25d0)/Lsm)
!!          else
!!             W(i,j,k,IV1) =    u1 - (  u1-  u2)/2.0d0*exp(-(-0.25d0-x2b(j))/Lsm)
!!             W(i,j,k,IDN) =  rho1 - (rho1-rho2)/2.0d0*exp(-(-0.25d0-x2b(j))/Lsm)
!!         endif
!!
!!          W(i,j,k,IPR) = 2.5d0
!!          W(i,j,k,IV2) = 0.01d0*sin(4.0d0*pi*x1b(i))
!!          W(i,j,k,IV3) = 0.0d0
!      enddo
!      enddo
!      enddo
!

      
!      call BoundaryCondition

      return
      end subroutine GenerateProblem

      subroutine BoundaryCondition(W)
      use commons
      implicit none
      real(8), intent(inout) :: W(:,:,:,:)
      integer::i,j,k

      do k=ks,ke
      do j=1,jn-1
      do i=1,mgn
!          W(is-i,j,k,IDN)  = W(is-1+i,j,k,IDN)
!         W(is-i,j,k,IV1)  = W(is-1+i,j,k,IV1)
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
          W(is-i,j,k,IB1)  = W(ie+1-i,j,k,IB1)
          W(is-i,j,k,IB2)  = W(ie+1-i,j,k,IB2)
          W(is-i,j,k,IB3)  = W(ie+1-i,j,k,IB3)
          W(is-i,j,k,IPS)  = W(ie+1-i,j,k,IPS)
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=1,jn-1
      do i=1,mgn
!          W(ie+i,j,k,IDN) = W(ie-i+1,j,k,IDN)
!         W(ie+i,j,k,IV1) = W(ie-i+1,j,k,IV1)
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
          W(ie+i,j,k,IB1) = W(is+i-1,j,k,IB1)
          W(ie+i,j,k,IB2) = W(is+i-1,j,k,IB2)
          W(ie+i,j,k,IB3) = W(is+i-1,j,k,IB3)
          W(ie+i,j,k,IPS) = W(is+i-1,j,k,IPS)
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
          W(i,js-j,k,IB1)  = W(i,je+1-j,k,IB1)
          W(i,js-j,k,IB2)  = W(i,je+1-j,k,IB2)
          W(i,js-j,k,IB3)  = W(i,je+1-j,k,IB3)
          W(i,js-j,k,IPS)  = W(i,je+1-j,k,IPS)
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
!         B(i,je+j,k,1)  = B(i,je-j+1,k,1)
!          B(i,je+j,k,2)  = B(i,je-j+1,k,2)
!          B(i,je+j,k,3)  = B(i,je-j+1,k,3)
!
          W(i,je+j,k,IDN)  = W(i,js+j-1,k,IDN)
          W(i,je+j,k,IV1)  = W(i,js+j-1,k,IV1)
          W(i,je+j,k,IV2)  = W(i,js+j-1,k,IV2)
          W(i,je+j,k,IV3)  = W(i,js+j-1,k,IV3)
          W(i,je+j,k,IPR)  = W(i,js+j-1,k,IPR)
          W(i,je+j,k,IB1)  = W(i,js+j-1,k,IB1)
          W(i,je+j,k,IB2)  = W(i,js+j-1,k,IB2)
          W(i,je+j,k,IB3)  = W(i,js+j-1,k,IB3)
          W(i,je+j,k,IPS)  = W(i,js+j-1,k,IPS)
      enddo
      enddo
      enddo



      return
      end subroutine BoundaryCondition
!
      subroutine ConsvVariable(W, U)
      use commons
      use eosmod
      implicit none
      real(8), intent(in) :: W(:,:,:,:)
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
                       + 0.5d0*( W(i,j,k,IB1)**2 + W(i,j,k,IB2)**2 + W(i,j,k,IB3)**2 ) &
                       + W(i,j,k,IPR)/(gam - 1.0d0)
          U(i,j,k,IB1) = W(i,j,k,IB1)
          U(i,j,k,IB2) = W(i,j,k,IB2)
          U(i,j,k,IB3) = W(i,j,k,IB3)
          U(i,j,k,IB3) = W(i,j,k,IB3)
          U(i,j,k,IPS) = W(i,j,k,IPS)
      enddo
      enddo
      enddo
      
      return
      end subroutine Consvvariable

      subroutine PrimVariable( U, W )
      use commons
      use eosmod
      implicit none
      real(8), intent(in) :: U(:,:,:,:)
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
                        - 0.5d0*(U(i,j,k,IB1)**2 + W(i,j,k,IB2)**2 + W(i,j,k,IB3)**2) )*(gam-1.0d0)
           W(i,j,k,IB1) = U(i,j,k,IB1)
           W(i,j,k,IB2) = U(i,j,k,IB2)
           W(i,j,k,IB3) = U(i,j,k,IB3)
           W(i,j,k,IPS) = U(i,j,k,IPS)
      enddo
      enddo
      enddo

      return
      end subroutine PrimVariable

      subroutine TimestepControl(x1a, x2a, x3a, W)
      use commons
      use eosmod
      implicit none
      real(8), intent(in) :: x1a(:), x2a(:), x3a(:), W(:,:,:,:)
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
         cf = dsqrt( (gam*W(i,j,k,IPR) + W(i,j,k,IB1)**2 + W(i,j,k,IB2)**2 + W(i,j,k,IB3)**2)/W(i,j,k,IDN))
         dtl1 =(x1a(i+1)-x1a(i))/(abs(W(i,j,k,IV1)) + cf)
         dtl2 =(x2a(j+1)-x2a(j))/(abs(W(i,j,k,IV2)) + cf)
!         dtl3 =(x3a(j+1)-x3a(j))/(abs(W(i,j,k,IV3)) + cf)
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
      subroutine NumericalFlux( W, flux1, flux2, flux3 )
      use commons !, only: is, ie, in
      implicit none
      integer::i,j,k
      real(8), intent(in) :: W(:,:,:,:)
      real(8), intent(out) :: flux1(:,:,:,:)
      real(8), intent(out) :: flux2(:,:,:,:)
      real(8), intent(out) :: flux3(:,:,:,:)
      real(8),dimension(in,jn,kn,NFLX):: Wl,Wr
      real(8),dimension(NFLX):: flx
      real(8) :: dWm(NFLX), dWp(NFLX), dWmon(NFLX)
      real(8) :: ddmon, dvmon, dpmon
      real(8) :: ch

      ch = Ccfl*min( x1a(is+1) - x1a(is), x2a(js+1) - x2a(js ) )/dt

      ! numerical flux in the x direction
      do k=ks,ke
      do j=js,je
      do i=is-1,ie+1
         dWp(1:NVAR) = W(i+1,j,k,1:NVAR) - W(i  ,j,k,1:NVAR)
         dWm(1:NVAR) = W(i  ,j,k,1:NVAR) - W(i-1,j,k,1:NVAR)

         call vanLeer(NFLX, dWp, dWm, dWmon)

         ! Wl(i,j,k) --> W_(i-1/2,j,k)
         ! Wr(i,j,k) --> W_(i-1/2,j,k)
         Wl(i+1,j,k,1:NVAR) = W(i,j,k,1:NVAR) + 0.5d0*dWmon(1:NVAR)*1.0d0
         Wr(i  ,j,k,1:NVAR) = W(i,j,k,1:NVAR) - 0.5d0*dWmon(1:NVAR)*1.0d0
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=js,je
      do i=is,ie+1
        call HLL(1,ch,Wl(i,j,k,:),Wr(i,j,k,:),flx)

         flux1(i,j,k,:)  = flx(:)
      enddo
      enddo
      enddo

      ! numerical flux in the y direction
      do k=ks,ke
      do j=js-1,je+1
      do i=is,ie
         dWp(1:NVAR) = W(i,j+1,k,1:NVAR) - W(i,j  ,k,1:NVAR)
         dWm(1:NVAR) = W(i,j  ,k,1:NVAR) - W(i,j-1,k,1:NVAR)

         call vanLeer(NFLX, dWp, dWm, dWmon)

         ! Wl(i,j,k) --> W_(i-1/2,j,k)
         ! Wr(i,j,k) --> W_(i-1/2,j,k)
         Wl(i,j+1,k,1:NVAR) = W(i,j,k,1:NVAR) + 0.5d0*dWmon(1:NVAR)*1.0d0
         Wr(i,j  ,k,1:NVAR) = W(i,j,k,1:NVAR) - 0.5d0*dWmon(1:NVAR)*1.0d0
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=js,je+1
      do i=is,ie
         call HLL(2,ch,Wl(i,j,k,:),Wr(i,j,k,:),flx)

         flux2(i,j,k,:)  = flx(:)
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
      subroutine HLL(idir,ch,Wl,Wr,flx)
      use commons !, only : is, ie, NVAR
      use eosmod
      implicit none
      integer, intent(in) :: idir
      real(8),intent(in)  :: ch
      real(8),intent(in)  ::Wl(:), Wr(:)
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
          b1 = 0.5d0*( Wl(IBpara) + Wr(IBpara) )
          
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
               flx(IBpara) = 0.0d0
               flx(IBperp1) = Fl(IBperp1)
               flx(IBperp2) = Fl(IBperp2)
          else if (sr <= 0.0d0 ) then
               flx(IDN) = Fr(IDN)
               flx(IVpara) = Fr(IM1)
               flx(IVperp1) = Fr(IM2)
               flx(IVperp2) = Fr(IM3)
               flx(IEN) = Fr(IEN)
               flx(IBpara) = 0.0d0
               flx(IBperp1) = Fr(IBperp1)
               flx(IBperp2) = Fr(IBperp2)
          else 
               flx(IDN)  = (sr*Fl(IDN) - sl*Fr(IDN) + sl*sr*( Ur(IDN) - Ul(IDN) ))/(sr - sl)
               flx(IVpara)  = (sr*Fl(IM1) - sl*Fr(IM1) + sl*sr*( Ur(IM1) - Ul(IM1) ))/(sr - sl)
               flx(IVperp1)  = (sr*Fl(IM2) - sl*Fr(IM2) + sl*sr*( Ur(IM2) - Ul(IM2) ))/(sr - sl)
               flx(IVperp2)  = (sr*Fl(IM3) - sl*Fr(IM3) + sl*sr*( Ur(IM3) - Ul(IM3) ))/(sr - sl)
               flx(IEN)  = (sr*Fl(IEN) - sl*Fr(IEN) + sl*sr*( Ur(IEN) - Ul(IEN) ))/(sr - sl)
               flx(IBpara) = 0.0d0
               flx(IBperp1)  = (sr*Fl(IBperp1) - sl*Fr(IBperp1) + sl*sr*( Ur(IBperp1) - Ul(IBperp1) ))/(sr - sl)
               flx(IBperp2)  = (sr*Fl(IBperp2) - sl*Fr(IBperp2) + sl*sr*( Ur(IBperp2) - Ul(IBperp2) ))/(sr - sl)
          endif

          flx(IBpara) = 0.5d0*(Wl(IPS) + Wr(IPS)) - 0.50d0*ch*(Wr(IBpara) - Wl(IBpara))
          flx(IPS) = ch*ch*( 0.5d0*( Wl(IBpara) + Wr(IBpara) ) - 0.5d0/ch*(Wr(IPS) - Wl(IPS)) )
!          flx(IBpara) = 0.0d0
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

!!=====================================================================

      subroutine UpdateConsv( x1a, x2a, x3a, flux1, flux2, flux3, U)
      use commons
      implicit none
      real(8), intent(in)  :: x1a(:), x2a(:), x3a(:)
      real(8), intent(in)  :: flux1(:,:,:,:), flux2(:,:,:,:), flux3(:,:,:,:)
      real(8), intent(out) :: U(:,:,:,:)
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

      do k=ks,ke
      do j=js,je
      do i=is,ie
         W(i,j,k,IPS) = W(i,j,k,IPS)*dexp( - 0.1d0*Ccfl )
      enddo
      enddo
      enddo


      return
      end subroutine UpdateConsv

      subroutine Output( x1a, x1b, x2a, x2b, W )
      use commons
      implicit none
      real(8), intent(in) :: x1a(:), x1b(:), x2a(:), x2b(:), W(:,:,:,:)
      integer::i,j,k
      character(20),parameter::dirname="snap/"
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

      write(filename,'(a4,i5.5,a8)')"snap",nout,"_256.dat"
      filename = trim(dirname)//filename
!      open(unitbin,file=filename,status='replace',form='formatted') 
      open(unitbin,file=filename,form='formatted',action="write")
      write(unitbin,"(a9,f4.2,a7,i3.3,a7,i3.3)") "# time = ",time, "  nx = ",nx, "  ny = ",ny
      do k=ks,ke
      do j=js,je
      do i=is,ie
          write(unitbin,*) x1b(i), x2b(j), W(i,j,k,IDN), W(i,j,k,IV1), W(i,j,k,IV2), W(i,j,k,IV3), W(i,j,k,IPR), W(i,j,k,IB1), &
          W(i,j,k,IB2), W(i,j,k,IB3),W(i,j,k,IPS)
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

      real(8) function divergenceB(x1a, x2a, W)
      use commons
      use eosmod
      implicit none
      real(8), intent(in) :: x1a(:), x2a(:), W(:,:,:,:)
      integer::i,j,k
      real(8) :: error

      do k=ks,ke
      do j=js,je
      do i=is,ie
         error = dabs( ( W(i+1,j,k,IB1) - W(i-1,j,k,IB1) )/(x1a(i+1)-x1a(i-1))  &
                    + ( W(i,j+1,k,IB2) - W(i,j-1,k,IB2) )/(x2a(j+1)-x2a(j-1)) )/ &
                    dsqrt( W(i,j,k,IB1)**2 + W(i,j,k,IB2)**2 )*min(x1a(i+1) - x1a(i),x2a(j+1)-x2a(j))

         divergenceB = divergenceB + error
      enddo
      enddo
      enddo
      divergenceB = divergenceB/dble((ie-is+1)*(je-js+1)*(ke-ks+1))
      
      return
      end function

      real(8) function dvy(x1a, x2a, W)
      use commons
      use eosmod
      implicit none
      real(8), intent(in) :: x1a(:), x2a(:), W(:,:,:,:)
      integer::i,j,k
      real(8) :: error

      dvy = 0.0d0
      do k=ks,ke
      do j=js,je
      do i=is,ie
           dvy = dvy + W(i,j,k,IV2)**2
      enddo
      enddo
      enddo
      dvy = dvy/dble(nx*ny)
      
      return
      end function
end program main
