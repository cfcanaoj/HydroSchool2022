module commons
implicit none
integer::ntime                        ! counter of the timestep
integer,parameter::ntimemax=200000     ! the maximum timesteps
real(8)::time,dt                    ! time, timewidth
data time / 0.0d0 /
real(8),parameter:: timemax=1.0d0/dsqrt(5.0d0)
real(8),parameter:: dtout=0.1d0/dsqrt(5.0d0)

integer,parameter::nx=1024/2/2/2/1/1/1 ! the number of grids in the simulation box
integer,parameter::ny=nx/2 ! the number of grids in the simulation box
integer,parameter::nz=1          ! the number of grids in the simulation box
integer,parameter::mgn=2         ! the number of ghost cells
integer,parameter::in=nx+2*mgn+1 ! the total number of grids including ghost cells
integer,parameter::jn=ny+2*mgn+1 ! the total number of grids including ghost cells
integer,parameter::kn=2 ! the total number of grids including ghost cells
integer,parameter::is=mgn+1         ! the index of the leftmost grid
integer,parameter::js=mgn+1         ! the index of the leftmost grid
integer,parameter::ks=1         ! the index of the leftmost grid
integer,parameter::ie=nx+mgn     ! the index of the rightmost grid
integer,parameter::je=ny+mgn     ! the index of the rightmost grid
integer,parameter::ke=1 ! the index of the rightmost grid
real(8),parameter::x1min=0.d0,x1max=1.0d0
real(8),parameter::x2min=0.d0,x2max=1.0d0*dble(ny)/dble(nx)
real(8),parameter::x3min=0.0d0,x3max=1.0d0

real(8),parameter::Ccfl=0.4d0/2/2/2/1/1/1

integer, parameter :: IDN = 1
integer, parameter :: IM1 = 2
integer, parameter :: IM2 = 3
integer, parameter :: IM3 = 4
integer, parameter :: IPR = 5
integer, parameter :: IB1 = 6
integer, parameter :: IB2 = 7
integer, parameter :: IB3 = 8
integer, parameter :: NVAR = 5
integer, parameter :: NFLX = 8

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
      real(8),dimension(in,jn,kn,NVAR) :: Uo
      real(8),dimension(in,jn,kn,NVAR) :: U
      real(8),dimension(in,jn,kn,NVAR) :: Q
      real(8),dimension(in,jn,kn,3) :: Bso ! magnetic field at cell surface
      real(8),dimension(in,jn,kn,3) :: Bs ! magnetic field at cell surface
      real(8),dimension(in,jn,kn,3) :: Bc ! magnetic field at cell center
      real(8),dimension(in,jn,kn,NVAR) :: F
      real(8),dimension(in,jn,kn,NVAR) :: G
      real(8),dimension(in,jn,kn,NVAR) :: H
      real(8),dimension(in,jn,kn,3) :: E ! electric field
      integer :: i, j,k

!      write(6,*) "setup grids and initial condition"
      call GenerateGrid(x1a, x1b, x2a, x2b, x3a, x3b)
      call GenerateProblem(x1b, x2b, x3b, Q, Bs, Bc)
      call ConsvVariable(Q, Bc, U)
      call BoundaryCondition( Q, Bs, Bc )
!      call Output( x1a, x1b, x2a, x2b, Q, Bc )


      print*, nx, Analysis( x1a, x1b, x2a, x2b, Q, Bc)

      open(1,file="Bpara.dat",action="write")
! main loop
      mloop: do !ntime=1,ntimemax
         call TimestepControl(x1a, x2a, x3a, Q, Bc)
         if( time + dt > timemax ) dt = timemax - time

         Uo(:,:,:,:) = U(:,:,:,:)
         Bso(:,:,:,:) = Bs(:,:,:,:)

         call NumericalFlux( x1a, x2a, x3a, Q, F, G, H, E )
         call UpdateConsv( 0.5d0*dt, x1a, x2a, x3a, F, G, H, E, Q, U, Bs, U, Bs )
         call PrimVariable( U, Bs, Q, Bc )
         call BoundaryCondition( Q, Bs, Bc )

         call NumericalFlux( x1a, x2a, x3a, Q, F, G, H, E )
         call UpdateConsv( dt, x1a, x2a, x3a, F, G, H, E, Q, Uo, Bso, U, Bs )
         call PrimVariable( U, Bs, Q, Bc )
         call BoundaryCondition( Q, Bs, Bc )



         time=time+dt
!         call Output( x1a, x1b, x2a, x2b, Q, Bc)
!         write(1,*) time, errorBpara(Q), divergenceB(x1a,x2a,Bc)

         print*,time

         if(time >= timemax) exit mloop
      enddo mloop
      close(1)
!      call Output( x1a, x1b, x2a, x2b, Q, Bc)

      print*, nx, Analysis( x1a, x1b, x2a, x2b, Q, Bc)

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

      subroutine GenerateProblem(x1b, x2b, x3b, Q, Bs, Bc )
      use commons
      use eosmod
      implicit none
      integer::i, j, k
      real(8), intent(in ) :: x1b(:), x2b(:), x3b(:)
      real(8), intent(out) :: Q(:,:,:,:)
      real(8), intent(out) :: Bs(:,:,:,:)
      real(8), intent(out) :: Bc(:,:,:,:)
      real(8) :: Az(in,jn,kn), k_dot_x
      real(8) :: cosa, sina, pi, x1m, Bpara, ca, kwave, omega, amp
      pi=acos(-1.0d0)

      cosa = 1.0d0/dsqrt(5.0d0)
      sina = 2.0d0/dsqrt(5.0d0)

      Bpara = 1.0d0
!      cosa = 1.0d0
!      sina = 0.0d0


      ca = Bpara
      kwave = 2.0d0*pi*sqrt(1.0d0 + (sina/cosa)**2)
      omega = ca*kwave
      amp = 0.1d0;


      do k=ks,ke
      do j=js,je
      do i=is,ie
          k_dot_x = kwave*( x1b(i)*cosa + x2b(j)*sina )
          Q(i,j,k,IDN) = 1.0d0
          Q(i,j,k,IV1) = - amp*dsin(k_dot_x)*sina
          Q(i,j,k,IV2) =   amp*dsin(k_dot_x)*cosa
          Q(i,j,k,IV3) =   amp*dcos(k_dot_x)
          Q(i,j,k,IPR) = 10.0d0
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=js-1,je+1
      do i=is-1,ie+1
          Az(i,j,k) = Bpara*( - x1a(i)*sina + x2a(j)*cosa ) &
                    + amp/kwave*dcos( kwave*( x1a(i)*cosa + x2a(j)*sina ) )
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=js,je
      do i=is,ie+1
          Bs(i,j,k,1) = (Az(i,j+1,k) - Az(i,j,k))/(x2a(j+1) - x2a(j));
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=js,je+1
      do i=is,ie
          Bs(i,j,k,2) = - (Az(i+1,j,k) - Az(i,j,k))/(x1a(i+1) - x1a(i));
      enddo
      enddo
      enddo

      do k=ks,ke+1
      do j=js,je
      do i=is,ie
          Bs(i,j,k,3) = 0.1*dcos(kwave*( x1b(i)*cosa + x2b(j)*sina ))
      enddo
      enddo
      enddo

      call CellCenterMagneticField(Bs, Bc)

!      do k=ks,ke
!      do j=js,je
!      do i=is,ie
!
!         if( x2b(j) < 0.0d0 ) then 
!             Q(i,j,k,IDN) = 1.08d0
!             Q(i,j,k,IV1) = 0.5d0
!             Q(i,j,k,IV2) = 1.2d0 
!             Q(i,j,k,IV3) = 0.01d0 ! 0.01d0
!             Q(i,j,k,IPR) = 0.95d0 !0.95d0
!
!             B(i,j,k,1) = 2.0d0/dsqrt(4.0d0*3.141592653589793)
!             B(i,j,k,2) = 4.0d0/dsqrt(4.0d0*3.141592653589793)
!             B(i,j,k,3) = 3.6d0/dsqrt(4.0d0*3.141592653589793)
!
!        else 
!             Q(i,j,k,IDN) = 1.0d0
!             Q(i,j,k,IV1) = 0.0d0
!             Q(i,j,k,IV2) = 0.0d0
!             Q(i,j,k,IV3) = 0.0d0
!             Q(i,j,k,IPR) = 1.0d0
!
!             B(i,j,k,1) = 2.0d0/dsqrt(4.0d0*3.141592653589793)
!             B(i,j,k,2) = 4.0d0/dsqrt(4.0d0*3.141592653589793)
!             B(i,j,k,3) = 4.0d0/dsqrt(4.0d0*3.141592653589793)
!         endif
!
!!         if( x1b(i) < 0.0d0 ) then 
!!             Q(i,j,k,IDN) = 1.08d0
!!             Q(i,j,k,IV1) = 1.2d0
!!             Q(i,j,k,IV2) = 0.01d0
!!             Q(i,j,k,IV3) = 0.5d0
!!             Q(i,j,k,IPR) = 0.95d0
!!
!!             B(i,j,k,1) = 4.0d0/dsqrt(4.0d0*3.141592653589793)
!!             B(i,j,k,2) = 3.6d0/dsqrt(4.0d0*3.141592653589793)
!!             B(i,j,k,3) = 2.0d0/dsqrt(4.0d0*3.141592653589793)
!!
!!         else 
!!             Q(i,j,k,IDN) = 1.0d0
!!             Q(i,j,k,IV1) = 0.0d0
!!             Q(i,j,k,IV2) = 0.0d0
!!             Q(i,j,k,IV3) = 0.0d0
!!             Q(i,j,k,IPR) = 1.0d0
!!
!!             B(i,j,k,1) = 4.0d0/dsqrt(4.0d0*3.141592653589793)
!!             B(i,j,k,2) = 4.0d0/dsqrt(4.0d0*3.141592653589793)
!!             B(i,j,k,3) = 2.0d0/dsqrt(4.0d0*3.141592653589793)
!!         endif
!!
!
!!         if      ( x2b(j) .gt. 0.25d0 )then
!!             Q(i,j,k,IV1) =    u1 - (  u1-  u2)/2.0d0*exp(-( x2b(j)-0.25d0)/Lsm)
!!             Q(i,j,k,IDN) =  rho1 - (rho1-rho2)/2.0d0*exp(-( x2b(j)-0.25d0)/Lsm)
!!         else if (x2b(j) .gt.  0.0d0 )then
!!             Q(i,j,k,IV1) =    u2 + (  u1-  u2)/2.0d0*exp(-( 0.25d0-x2b(j))/Lsm)
!!             Q(i,j,k,IDN) =  rho2 + (rho1-rho2)/2.0d0*exp(-( 0.25d0-x2b(j))/Lsm)
!!         else if (x2b(j) .gt. -0.25d0)then
!!             Q(i,j,k,IV1) =    u2 + (  u1-  u2)/2.0d0*exp(-( x2b(j)+0.25d0)/Lsm)
!!             Q(i,j,k,IDN) =  rho2 + (rho1-rho2)/2.0d0*exp(-( x2b(j)+0.25d0)/Lsm)
!!          else
!!             Q(i,j,k,IV1) =    u1 - (  u1-  u2)/2.0d0*exp(-(-0.25d0-x2b(j))/Lsm)
!!             Q(i,j,k,IDN) =  rho1 - (rho1-rho2)/2.0d0*exp(-(-0.25d0-x2b(j))/Lsm)
!!         endif
!!
!!          Q(i,j,k,IPR) = 2.5d0
!!          Q(i,j,k,IV2) = 0.01d0*sin(4.0d0*pi*x1b(i))
!!          Q(i,j,k,IV3) = 0.0d0
!      enddo
!      enddo
!      enddo
!

      
!      call BoundaryCondition

      return
      end subroutine GenerateProblem

      subroutine BoundaryCondition( Q, Bs, Bc )
      use commons
      implicit none
      real(8), intent(inout) :: Q(:,:,:,:)
      real(8), intent(inout) :: Bs(:,:,:,:)
      real(8), intent(inout) :: Bc(:,:,:,:)
      integer::i,j,k

      ! x inner boundary
      do k=ks,ke
      do j=js-mgn,je+mgn
      do i=1,mgn
          Q(is-i,j,k,IDN)  = Q(ie+1-i,j,k,IDN)
          Q(is-i,j,k,IV1)  = Q(ie+1-i,j,k,IV1)
          Q(is-i,j,k,IV2)  = Q(ie+1-i,j,k,IV2)
          Q(is-i,j,k,IV3)  = Q(ie+1-i,j,k,IV3)
          Q(is-i,j,k,IPR)  = Q(ie+1-i,j,k,IPR)
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
          Q(ie+i,j,k,IDN) = Q(is+i-1,j,k,IDN)
          Q(ie+i,j,k,IV1) = Q(is+i-1,j,k,IV1)
          Q(ie+i,j,k,IV2) = Q(is+i-1,j,k,IV2)
          Q(ie+i,j,k,IV3) = Q(is+i-1,j,k,IV3)
          Q(ie+i,j,k,IPR) = Q(is+i-1,j,k,IPR)
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
          Q(i,js-j,k,IDN)  = Q(i,je+1-j,k,IDN)
          Q(i,js-j,k,IV1)  = Q(i,je+1-j,k,IV1)
          Q(i,js-j,k,IV2)  = Q(i,je+1-j,k,IV2)
          Q(i,js-j,k,IV3)  = Q(i,je+1-j,k,IV3)
          Q(i,js-j,k,IPR)  = Q(i,je+1-j,k,IPR)
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

          Q(i,je+j,k,IDN)  = Q(i,js+j-1,k,IDN)
          Q(i,je+j,k,IV1)  = Q(i,js+j-1,k,IV1)
          Q(i,je+j,k,IV2)  = Q(i,js+j-1,k,IV2)
          Q(i,je+j,k,IV3)  = Q(i,js+j-1,k,IV3)
          Q(i,je+j,k,IPR)  = Q(i,js+j-1,k,IPR)
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


      ! cell center
      do k=ks,ke
      do j=js-mgn,je+mgn
      do i=is-mgn,is-1
          Bc(i,j,k,1) = 0.5d0*(Bs(i,j,k,1) + Bs(i+1,j,k,1))
          Bc(i,j,k,2) = 0.5d0*(Bs(i,j,k,2) + Bs(i,j+1,k,2))
          Bc(i,j,k,3) = 0.5d0*(Bs(i,j,k,3) + Bs(i,j,k+1,3))
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=js-mgn,je+mgn
      do i=ie+1,ie+mgn
          Bc(i,j,k,1) = 0.5d0*(Bs(i,j,k,1) + Bs(i+1,j,k,1))
          Bc(i,j,k,2) = 0.5d0*(Bs(i,j,k,2) + Bs(i,j+1,k,2))
          Bc(i,j,k,3) = 0.5d0*(Bs(i,j,k,3) + Bs(i,j,k+1,3))
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=js-mgn,js-1
      do i=is-mgn,ie+mgn
          Bc(i,j,k,1) = 0.5d0*(Bs(i,j,k,1) + Bs(i+1,j,k,1))
          Bc(i,j,k,2) = 0.5d0*(Bs(i,j,k,2) + Bs(i,j+1,k,2))
          Bc(i,j,k,3) = 0.5d0*(Bs(i,j,k,3) + Bs(i,j,k+1,3))
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=je+1,je+mgn
      do i=is-mgn,ie+mgn
          Bc(i,j,k,1) = 0.5d0*(Bs(i,j,k,1) + Bs(i+1,j,k,1))
          Bc(i,j,k,2) = 0.5d0*(Bs(i,j,k,2) + Bs(i,j+1,k,2))
          Bc(i,j,k,3) = 0.5d0*(Bs(i,j,k,3) + Bs(i,j,k+1,3))
      enddo
      enddo
      enddo

!      do k=ks,ke
!      do j=js-mgn,je+mgn+1
!      do i=is-mgn,ie+mgn
!          print*,x1b(i),x2a(j),Bs(i,j,k,2)
!      enddo
!      enddo
!      enddo
!      stop
!
      return
      end subroutine BoundaryCondition
!
      subroutine ConsvVariable(Q, Bc, U)
      use commons
      use eosmod
      implicit none
      real(8), intent(in) :: Q(:,:,:,:)
      real(8), intent(in) :: Bc(:,:,:,:)
      real(8), intent(out) :: U(:,:,:,:)
      integer::i,j,k

      do k=ks,ke
      do j=js,je
      do i=is,ie
          U(i,j,k,IDN) = Q(i,j,k,IDN)
          U(i,j,k,IM1) = Q(i,j,k,IDN)*Q(i,j,k,IV1)
          U(i,j,k,IM2) = Q(i,j,k,IDN)*Q(i,j,k,IV2)
          U(i,j,k,IM3) = Q(i,j,k,IDN)*Q(i,j,k,IV3)
          U(i,j,k,IEN) = 0.5d0*Q(i,j,k,IDN)*( Q(i,j,k,IV1)**2 + Q(i,j,k,IV2)**2 + Q(i,j,k,IV3)**2 ) &
                       + 0.5d0*( Bc(i,j,k,1)**2 + Bc(i,j,k,2)**2 + Bc(i,j,k,3)**2 ) &
                       + Q(i,j,k,IPR)/(gam - 1.0d0)
      enddo
      enddo
      enddo
      
      return
      end subroutine Consvvariable

      subroutine PrimVariable( U, Bs, Q, Bc )
      use commons
      use eosmod
      implicit none
      real(8), intent(in) :: U(:,:,:,:),Bs(:,:,:,:)
      real(8), intent(out) :: Q(:,:,:,:),Bc(:,:,:,:)
      integer::i,j,k
      real(8) :: inv_d;

      do k=ks,ke
      do j=js,je
      do i=is,ie
           Bc(i,j,k,1) = 0.5d0*( Bs(i+1,j,k,1) + Bs(i,j,k,1) )
           Bc(i,j,k,2) = 0.5d0*( Bs(i,j+1,k,2) + Bs(i,j,k,2) )
           Bc(i,j,k,3) = 0.5d0*( Bs(i,j,k+1,3) + Bs(i,j,k,3) )
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=js,je
      do i=is,ie
           Q(i,j,k,IDN) = U(i,j,k,IDN)
           inv_d = 1.0d0/U(i,j,k,IDN)
           Q(i,j,k,IV1) = U(i,j,k,IM1)*inv_d
           Q(i,j,k,IV2) = U(i,j,k,IM2)*inv_d
           Q(i,j,k,IV3) = U(i,j,k,IM3)*inv_d
           Q(i,j,k,IPR) = ( U(i,j,k,IEN) &
                        - 0.5d0*(U(i,j,k,IM1)**2 + U(i,j,k,IM2)**2 + U(i,j,k,IM3)**2)*inv_d  &
                        - 0.5d0*(Bc(i,j,k,1)**2 + Bc(i,j,k,2)**2 + Bc(i,j,k,3)**2) )*(gam-1.0d0)
      enddo
      enddo
      enddo

      return
      end subroutine PrimVariable

      subroutine CellCenterMagneticField( Bs, Bc )
      use commons
      use eosmod
      implicit none
      real(8), intent(in) :: Bs(:,:,:,:)
      real(8), intent(out) :: Bc(:,:,:,:)
      integer::i,j,k
      real(8) :: inv_d;

      do k=ks,ke
      do j=js,je
      do i=is,ie
             Bc(i,j,k,1) = 0.5d0*( Bs(i+1,j,k,1) + Bs(i,j,k,1) )
             Bc(i,j,k,2) = 0.5d0*( Bs(i,j+1,k,2) + Bs(i,j,k,2) )
             Bc(i,j,k,3) = 0.5d0*( Bs(i,j,k+1,3) + Bs(i,j,k,3) )
      enddo
      enddo
      enddo

      return
      end subroutine CellCenterMagneticField

      subroutine TimestepControl(x1a, x2a, x3a, Q, Bc)
      use commons
      use eosmod
      implicit none
      real(8), intent(in) :: x1a(:), x2a(:), x3a(:), Q(:,:,:,:), Bc(:,:,:,:)
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
         cf = dsqrt( (gam*Q(i,j,k,IPR) + Bc(i,j,k,1)**2 + Bc(i,j,k,2)**2 + Bc(i,j,k,3)**2)/Q(i,j,k,IDN))
         dtl1 =(x1a(i+1)-x1a(i))/(abs(Q(i,j,k,IV1)) + cf)
         dtl2 =(x2a(j+1)-x2a(j))/(abs(Q(i,j,k,IV2)) + cf)
!         dtl3 =(x3a(j+1)-x3a(j))/(abs(Q(i,j,k,IV3)) + cf)
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
!            sgn = 1.0d0
!            if (dvp(i) < 0.0d0) then
!                sgn = -1.0d0
!            endif
!            dv(i) = sgn*min( 2.0d0*dabs(dvm(i)), 2.0d0*dabs(dvp(i)), 0.5d0*dabs(dvp(i) + dvm(i)) )
         else
            dv(i) = 0.0d0
         endif
!         dv(i) = 0.5d0*(dvp(i) + dvm(i))
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
      subroutine NumericalFlux( x1a, x2a, x3a, Q, F, G, H, E)
      use commons !, only: is, ie, in
      implicit none
      integer::i,j,k
      real(8), intent(in) :: x1a(:), x2a(:), x3a(:)
      real(8), intent(in) :: Q(:,:,:,:)
      real(8), intent(out) :: F(:,:,:,:)
      real(8), intent(out) :: G(:,:,:,:)
      real(8), intent(out) :: H(:,:,:,:)
      real(8), intent(out) :: E(:,:,:,:)
      real(8),dimension(in,jn,kn,NFLX):: Ql,Qr
      real(8),dimension(NFLX):: flx
      real(8) :: dQm(NVAR), dQp(NVAR), dQmon(NVAR)
      real(8) :: ddmon, dvmon, dpmon

      real(8),dimension(in,jn,kn) :: e2_x1a, e3_x1a
      real(8),dimension(in,jn,kn) :: e1_x2a, e3_x2a
      real(8),dimension(in,jn,kn) :: weight1, weight2, weight3
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

      do k=ks,ke
      do j=js-1,je+1
      do i=is-1,ie+1
!        call HLL(1,Ql(i,j,k,:),Qr(i,j,k,:),flx)
        call HLLD(1,Ql(i,j,k,:),Qr(i,j,k,:),x1a(i+1)-x1a(i),flx,wghtCT)

         F(i,j,k,1:NVAR)  = flx(1:NVAR)
         e3_x1a(i,j,k) =  -flx(IB2)
         e2_x1a(i,j,k) =   flx(IB3)

         weight1(i,j,k) = wghtCT
      enddo
      enddo
      enddo


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


      do k=ks,ke
      do j=js-1,je+1
      do i=is-1,ie+1
!         call HLL(2,Ql(i,j,k,:),Qr(i,j,k,:),flx)
         call HLLD(2,Ql(i,j,k,:),Qr(i,j,k,:),x2a(j+1) - x2a(j), flx,wghtCT)

         G(i,j,k,1:NVAR) = flx(1:NVAR)
         e1_x2a(i,j,k) =  - flx(IB3)
         e3_x2a(i,j,k) =    flx(IB1)

         weight2(i,j,k) = wghtCT

      enddo
      enddo
      enddo

      call ElectricField( Q, Bc, e2_x1a, e3_x1a, e3_x2a, e1_x2a, weight1, weight2, E )

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
      use commons !, only : is, ie, NVAR
      use eosmod
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
          b1 = 0.5d0*( Ql(IBpara) + Qr(IBpara) )
          
          pbl = 0.5d0*(b1**2 + Ql(IBperp1)**2 + Ql(IBperp2)**2)
          pbr = 0.5d0*(b1**2 + Qr(IBperp1)**2 + Qr(IBperp2)**2)
          ptotl = Ql(IPR) + pbl
          ptotr = Qr(IPR) + pbr

          ! conserved variables in the left and right states
          Ul(IDN) = Ql(IDN)
          Ul(IM1) = Ql(IDN)*Ql(IVpara)
          Ul(IM2) = Ql(IDN)*Ql(IVperp1)
          Ul(IM3) = Ql(IDN)*Ql(IVperp2)
          Ul(IEN) = 0.5d0*Ql(IDN)*( Ql(IVpara)**2 + Ql(IVperp1)**2 + Ql(IVperp2)**2) & 
                  + pbl + Ql(IPR)/(gam - 1.0d0)
          Ul(IBperp1) = Ql(IBperp1)
          Ul(IBperp2) = Ql(IBperp2)

          Ur(IDN) = Qr(IDN)
          Ur(IM1) = Qr(IDN)*Qr(IVpara)
          Ur(IM2) = Qr(IDN)*Qr(IVperp1)
          Ur(IM3) = Qr(IDN)*Qr(IVperp2)
          Ur(IEN) = 0.5d0*Qr(IDN)*( Qr(IVpara)**2 + Qr(IVperp1)**2 + Qr(IVperp2)**2) & 
                  + pbr + Qr(IPR)/(gam - 1.0d0)
          Ur(IBperp1) = Qr(IBperp1)
          Ur(IBperp2) = Qr(IBperp2)

          ! flux in the left and right states
          Fl(IDN) = Ul(IM1)
          Fl(IM1) = Ql(IDN)*Ql(IVpara)**2 + ptotl - b1**2
          Fl(IM2) = Ql(IDN)*Ql(IVpara)*Ql(IVperp1) - b1*Ql(IBperp1)
          Fl(IM3) = Ql(IDN)*Ql(IVpara)*Ql(IVperp2) - b1*Ql(IBperp2)
          Fl(IEN) = ( Ul(IEN) + ptotl )*Ql(IVpara) &
                  - b1*( b1*Ql(IVpara) + Ql(IBperp1)*Ql(IVperp1) + Ql(IBperp2)*Ql(IVperp2) )
          Fl(IBperp1) = Ql(IBperp1)*Ql(IVpara) - b1*Ql(IVperp1)
          Fl(IBperp2) = Ql(IBperp2)*Ql(IVpara) - b1*Ql(IVperp2)



          Fr(IDN) = Ur(IM1)
          Fr(IM1) = Qr(IDN)*Qr(IVpara)**2 + ptotr - b1**2
          Fr(IM2) = Qr(IDN)*Qr(IVpara)*Qr(IVperp1) - b1*Qr(IBperp1)
          Fr(IM3) = Qr(IDN)*Qr(IVpara)*Qr(IVperp2) - b1*Qr(IBperp2)
          Fr(IEN) = ( Ur(IEN) + ptotr )*Qr(IVpara) &
                  - b1*( b1*Qr(IVpara) + Qr(IBperp1)*Qr(IVperp1) + Qr(IBperp2)*Qr(IVperp2) )
          Fr(IBperp1) = Qr(IBperp1)*Qr(IVpara) - b1*Qr(IVperp1)
          Fr(IBperp2) = Qr(IBperp2)*Qr(IVpara) - b1*Qr(IVperp2)

!         print*, Qr(IDN)*Qr(IVpara)*Qr(IVperp2) - b1*Qr(IBperp2)
!         print*, Ql(IDN)*Ql(IVpara)*Ql(IVperp2) - b1*Ql(IBperp2)


          cfl = dsqrt( (gam*Ql(IPR) + Ql(IBperp1)**2 + Ql(IBperp2)**2 + b1**2)/Ql(IDN))
          cfr = dsqrt( (gam*Qr(IPR) + Qr(IBperp1)**2 + Qr(IBperp2)**2 + b1**2)/Qr(IDN))

!          sl = min(Ql(IV1),Qr(IV1)) - max(cfl,cfr)
!          sr = max(Ql(IV1),Qr(IV1)) + max(cfl,cfr)
          sl = min(Ql(IVpara) - cfl,Qr(IVpara) - cfr)
          sr = max(Ql(IVpara) + cfl,Qr(IVpara) + cfr)
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
!

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
      subroutine HLLD(idir,Ql,Qr,dx,flx,wghtCT)
      use commons !, only : is, ie, NVAR
      use eosmod
      implicit none
      integer, intent(in) :: idir
      real(8),intent(in)  :: Ql(:), Qr(:)
      real(8),intent(in)  :: dx
      real(8),intent(out) :: flx(:)
      real(8),intent(out) :: wghtCT
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
      real(8) :: v_over_c
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
      subroutine ElectricField( Q, Bc, e2_x1a, e3_x1a, e3_x2a, e1_x2a, weight1, weight2, E )
      use commons !, only: is, ie, in
      implicit none
      integer::i,j,k
      real(8), intent(in)  :: Q(:,:,:,:), Bc(:,:,:,:)
      real(8), intent(in)  :: e2_x1a(:,:,:)
      real(8), intent(in)  :: e3_x1a(:,:,:)
      real(8), intent(in)  :: e1_x2a(:,:,:)
      real(8), intent(in)  :: e3_x2a(:,:,:)
      real(8), intent(in)  :: weight1(:,:,:)
      real(8), intent(in)  :: weight2(:,:,:)
      real(8), intent(out) :: E(:,:,:,:)
      real(8) :: Etmp(in,jn,kn) 
      real(8) :: de3_l1, de3_r1, de3_l2, de3_r2

      do k=ks,ke
      do j=js-1, je+1
      do i=is-1, ie+1
           Etmp(i,j,k) = Q(i,j,k,IV2)*Bc(i,j,k,1) - Q(i,j,k,IV1)*Bc(i,j,k,2)
      enddo
      enddo
      enddo

      do j=js, je
      do i=is, ie+1
           E(i,j,ke+1,2) = e2_x1a(i,j,ks)
           E(i,j,ks  ,2) = e2_x1a(i,j,ks)
      enddo
      enddo

      do j=js, je+1
      do i=is, ie
           E(i,j,ke+1,1) = e1_x2a(i,j,ks)
           E(i,j,ks  ,1) = e1_x2a(i,j,ks)
      enddo
      enddo

      do k=ks,ke
      do j=js, je+1
      do i=is, ie+1
          de3_l2 = (1.0-weight1(i,j-1,k))*(e3_x2a(i  ,j,k) - Etmp(i  ,j-1,k)) + &
                   (    weight1(i,j-1,k))*(e3_x2a(i-1,j,k) - Etmp(i-1,j-1,k))

          de3_r2 = (1.0-weight1(i,j,k))*(e3_x2a(i,j,k) - Etmp(i,j,k)) + &
                   (    weight1(i,j,k))*(e3_x2a(i-1,j,k) - Etmp(i-1,j,k))

          de3_l1 = (1.0-weight2(i-1,j,k))*(e3_x1a(i,j,k) - Etmp(i-1,j,k)) + &
                   (    weight2(i-1,j,k))*(e3_x1a(i,j-1,k) - Etmp(i-1,j-1,k))

          de3_r1 = (1.0-weight2(i,j,k))*(e3_x1a(i,j,k) - Etmp(i,j,k)) + &
                   (    weight2(i,j,k))*(e3_x1a(i,j-1,k) - Etmp(i,j-1,k  ))

           E(i,j,k,3) = 0.25d0*( de3_l1 + de3_r1 + de3_l2 + de3_r2 + &
                              e3_x2a(i-1,j,k) + e3_x2a(i,j,k) + e3_x1a(i,j-1,k) + e3_x1a(i,j,k))
!                          print*,i,j,k,e3_x2a(i-1,j,k) , e3_x2a(i,j,k) , e3_x1a(i,j-1,k) , e3_x1a(i,j,k)

!          E(i,j,k,3) = 0.25d0*( e3_x2a(i-1,j,k) + e3_x2a(i,j,k) + e3_x1a(i,j-1,k) + e3_x1a(i,j,k))
      enddo
      enddo
      enddo

!      do k=ks,ke
!      do j=js, je-1
!      do i=is, ie
!           print*, "write",E(i,j+1,k,3) - E(i,j,k,3)
!      enddo
!      enddo
!      enddo
!


      return
      end subroutine ElectricField 
!!=====================================================================

      subroutine UpdateConsv( dt1, x1a, x2a, x3a, F, G, H, E, Q, Uo, Bso, U, Bs)
      use commons
      implicit none
      real(8), intent(in) :: dt1
      real(8), intent(in)  :: x1a(:), x2a(:), x3a(:)
      real(8), intent(in)  :: F(:,:,:,:), G(:,:,:,:), H(:,:,:,:)
      real(8), intent(in)  :: Uo(:,:,:,:), Q(:,:,:,:), Bso(:,:,:,:)
      real(8), intent(out) :: U(:,:,:,:), Bs(:,:,:,:), E(:,:,:,:)

      real(8) :: divB 
      integer::i,n,j,k

      do n=1,NVAR
      do k=ks,ke
      do j=js,je
      do i=is,ie
         U(i,j,k,n) = Uo(i,j,k,n) + dt1*(- F(i+1,j,k,n) + F(i,j,k,n))/(x1a(i+1)-x1a(i)) &
                                  + dt1*(- G(i,j+1,k,n) + G(i,j,k,n))/(x2a(j+1)-x2a(j)) 

!        if(n==1) print*,x1b(i),- F(i+1,j,k,n) + F(i,j,k,n)
      enddo
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=js,je
      do i=is,ie+1
           Bs(i,j,k,1) = Bso(i,j,k,1) &
                       - dt1*(E(i,j+1,k,3) - E(i,j,k,3))/(x2b(j+1) - x2b(j))
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=js,je+1
      do i=is,ie
           Bs(i,j,k,2) = Bso(i,j,k,2) &
                       + dt1*(E(i+1,j,k,3) - E(i,j,k,3))/(x1b(i+1) - x1b(i))
      enddo
      enddo
      enddo

      do k=ks,ke+1
      do j=js,je
      do i=is,ie
           Bs(i,j,k,3) = Bso(i,j,k,3) &
                       - dt1*(E(i+1,j,k,2) - E(i,j,k,2))/(x1b(i+1) - x1b(i)) &
                       + dt1*(E(i,j+1,k,1)  - E(i,j,k,1))/(x2b(j+1) - x2b(j))
      enddo
      enddo
      enddo



      return
      end subroutine UpdateConsv

      subroutine Output( x1a, x1b, x2a, x2b, Q, Bc )
      use commons
      implicit none
      real(8), intent(in) :: x1a(:), x1b(:), x2a(:), x2b(:), Q(:,:,:,:), Bc(:,:,:,:)
      integer::i,j,k
      character(20),parameter::dirname="snap_alf/"
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
      real(8) :: divB

      if (.not. is_inited) then
         call makedirs(dirname)
         is_inited =.true.
      endif

!      print*, time, tout+dtout, time+1.0d-14 .lt. tout+dtout
      if(time + 1.0d-14.lt. tout+dtout) return

      write(filename,'(a4,i5.5,a7)')"snap",nout,"_ct.dat"
      filename = trim(dirname)//trim(filename)
!      open(unitbin,file=filename,status='replace',form='formatted') 
      open(unitbin,file=filename,form='formatted',action="write")
      write(unitbin,"(a9,f4.2,a7,i3.3,a7,i3.3)") "# time = ",time, "  nx = ",nx, "  ny = ",ny
      do k=ks,ke
      do j=js,je
      do i=is-1,ie+1
         divB = ( Bc(i+1,j,k,1) - Bc(i-1,j,k,1) )/(x1b(i+1) - x1b(i-1)) &
              + ( Bc(i,j+1,k,2) - Bc(i,j-1,k,2) )/(x2b(j+1) - x2b(j-1)) 

          write(unitbin,*) x1b(i), x2b(j), Q(i,j,k,IDN), Q(i,j,k,IV1), Q(i,j,k,IV2), Q(i,j,k,IV3), Q(i,j,k,IPR), &
              Bc(i,j,k,1), Bc(i,j,k,2), Bc(i,j,k,3), divB
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

      real(8) function errorBpara(Q)
      use commons
      use eosmod
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

      real(8) function divergenceB(x1a, x2a, Bc)
      use commons
      use eosmod
      implicit none
      real(8), intent(in) :: x1a(:), x2a(:), Bc(:,:,:,:)
      integer::i,j,k
      real(8) :: error

      do k=ks,ke
      do j=js,je
      do i=is,ie
         error = dabs( ( Bc(i+1,j,k,1) - Bc(i-1,j,k,1) )/(x1a(i+1)-x1a(i-1))  &
                    + ( Bc(i,j+1,k,2) - Bc(i,j-1,k,2) )/(x2a(j+1)-x2a(j-1)) )/ &
                    dsqrt( Bc(i,j,k,IB1)**2 + Bc(i,j,k,2)**2 )*min(x1a(i+1) - x1a(i),x2a(j+1)-x2a(j))

         divergenceB = divergenceB + error
      enddo
      enddo
      enddo
      divergenceB = divergenceB/dble((ie-is+1)*(je-js+1)*(ke-ks+1))
      
      return
      end function

      real(8) function dvy(x1a, x2a, Q)
      use commons
      use eosmod
      implicit none
      real(8), intent(in) :: x1a(:), x2a(:), Q(:,:,:,:)
      integer::i,j,k
      real(8) :: error

      dvy = 0.0d0
      do k=ks,ke
      do j=js,je
      do i=is,ie
           dvy = dvy + Q(i,j,k,IV2)**2
      enddo
      enddo
      enddo
      dvy = dvy/dble(nx*ny)
      
      return
      end function

      real(8) function Analysis( x1a, x1b, x2a, x2b, Q, Bc)
      use commons
      use eosmod
      implicit none
      real(8), intent(in) :: x1a(:), x2a(:), x1b(:), x2b(:), Bc(:,:,:,:), Q(:,:,:,:)
      integer::i,j,k
      real(8) :: error, sina, cosa, k_dot_x, Bperp, kwave, vperp

      kwave = 2.0d0*dacos(-1.0d0)*dsqrt(1.0d0 + 2.0d0**2)
      cosa = 1.0d0/dsqrt(5.0d0)
      sina = 2.0d0/dsqrt(5.0d0)

      error= 0.0d0
      do k=ks,ke
      do j=js,je
      do i=is,ie
         k_dot_x = kwave*( x1b(i)*cosa + x2b(j)*sina ) 
         Bperp = - sina*Bc(i,j,k,1) + cosa*Bc(i,j,k,2)
         vperp = - sina*Q(i,j,k,IV1) + cosa*Q(i,j,k,IV2)

         error = error + ( Bperp - 0.1d0*dsin(k_dot_x) )**2 &
                       + ( vperp - 0.1d0*dsin(k_dot_x) )**2
      enddo
      enddo
      enddo

      Analysis = dsqrt(error/( dble(ie - is+1)*dble(je - js + 1) ))
      
      return
      end function

end program main
