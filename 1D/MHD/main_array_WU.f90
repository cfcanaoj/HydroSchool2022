module commons
implicit none
integer::ntime                        ! counter of the timestep
integer,parameter::ntimemax=200000     ! the maximum timesteps
real(8)::time,dt                    ! time, timewidth
data time / 0.0d0 /
real(8),parameter:: timemax=0.2d0   
real(8),parameter:: dtout=5.0d-3

integer,parameter::ngrid=128       ! the number of grids in the simulation box
integer,parameter::mgn=2            ! the number of ghost cells
integer,parameter::in=ngrid+2*mgn+1 ! the total number of grids including ghost cells
integer,parameter::is=mgn+1         ! the index of the leftmost grid
integer,parameter::ie=ngrid+mgn     ! the index of the rightmost grid
integer,parameter::nflux = 3
real(8),parameter:: x1min=-0.5d0,x1max=0.5d0

integer, parameter :: IDN = 1
integer, parameter :: IM1 = 2
integer, parameter :: IM2 = 3
integer, parameter :: IM3 = 4
integer, parameter :: IPR = 5
integer, parameter :: NVAR = 5

integer, parameter :: IB2 = 6
integer, parameter :: IB3 = 7
integer, parameter :: NFLX = 7

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
      real(8),dimension(in,NVAR) :: U
      real(8),dimension(in,NVAR) :: W
      real(8),dimension(in,3) :: B
      real(8),dimension(in,NFLX) :: flux

!      write(6,*) "setup grids and initial condition"
      call GenerateGrid(x1a, x1b)
      call GenerateProblem(x1b, W, B)
      call ConsvVariable(W, B, U)
      call Output( x1a, x1b, W, B )

! main loop
      mloop: do ntime=1,ntimemax
         call TimestepControl(x1a, W, B)
         if( time + dt > timemax ) dt = timemax - time
         call BoundaryCondition( W, B )
         call NumericalFlux( W, B, flux )
         call UpdateConsv( flux, U, B )
         call PrimVariable( U, B, W )
         time=time+dt
         call Output( x1a, x1b, W, B)

         if(time >= timemax) exit mloop
      enddo mloop
      call Output( x1a, x1b, W, B)

!      write(6,*) "program has been finished"
contains

      subroutine GenerateGrid(x1a, x1b)
      use commons
      use eosmod
      implicit none
      real(8), intent(out) :: x1a(:), x1b(:)
      real(8) :: dx,dy
      integer::i
      dx=(x1max-x1min)/ngrid
      do i=1,in
         x1a(i) = dx*(i-(mgn+1))+x1min
      enddo
      do i=1,in-1
         x1b(i) = 0.5d0*(x1a(i+1)+x1a(i))
      enddo

      return
      end subroutine GenerateGrid

      subroutine GenerateProblem(x1b, W, B)
      use commons
      use eosmod
      implicit none
      integer::i
      real(8), intent(in ) :: x1b(:)
      real(8), intent(out) :: W(:,:), B(:,:)
      real(8) :: rho1,rho2,Lsm,u1,u2

      real(8)::pi
      pi=acos(-1.0d0)

      do i=is,ie
         if( x1b(i) < 0.0d0 ) then 
             W(i,IDN) = 1.08d0
             W(i,IV1) = 1.2d0
             W(i,IV2) = 0.01d0
             W(i,IV3) = 0.5d0
             W(i,IPR) = 0.95d0

             B(i,1) = 4.0d0/dsqrt(4.0d0*3.141592653589793)
             B(i,2) = 3.6d0/dsqrt(4.0d0*3.141592653589793)
             B(i,3) = 2.0d0/dsqrt(4.0d0*3.141592653589793)

         else 
             W(i,IDN) = 1.0d0
             W(i,IV1) = 0.0d0
             W(i,IV2) = 0.0d0
             W(i,IV3) = 0.0d0
             W(i,IPR) = 1.0d0

             B(i,1) = 4.0d0/dsqrt(4.0d0*3.141592653589793)
             B(i,2) = 4.0d0/dsqrt(4.0d0*3.141592653589793)
             B(i,3) = 2.0d0/dsqrt(4.0d0*3.141592653589793)
         endif
      enddo

      
!      call BoundaryCondition

      return
      end subroutine GenerateProblem

      subroutine BoundaryCondition(W,B)
      use commons
      implicit none
      real(8), intent(inout) :: W(:,:), B(:,:)
      integer::i

      do i=1,mgn
          W(is-i,IDN)  = W(is-1+i,IDN)
          W(is-i,IV1)  = W(is-1+i,IV1)
          W(is-i,IV2)  = W(is-1+i,IV2)
          W(is-i,IV3)  = W(is-1+i,IV3)
          W(is-i,IPR)  = W(is-1+i,IPR)
          B(is-i,1)  = B(is-1+i,1)
          B(is-i,2)  = B(is-1+i,2)
          B(is-i,3)  = B(is-1+i,3)
      enddo

      do i=1,mgn
          W(ie+i,IDN) = W(ie-i+1,IDN)
          W(ie+i,IV1) = W(ie-i+1,IV1)
          W(ie+i,IV2) = W(ie-i+1,IV2)
          W(ie+i,IV3) = W(ie-i+1,IV3)
          W(ie+i,IPR) = W(ie-i+1,IPR)
          B(ie+i,1) = B(ie-i+1,1)
          B(ie+i,2) = B(ie-i+1,2)
          B(ie+i,3) = B(ie-i+1,3)
      enddo

      return
      end subroutine BoundaryCondition
!
      subroutine ConsvVariable(W, B, U)
      use commons
      use eosmod
      implicit none
      real(8), intent(in) :: W(:,:), B(:,:)
      real(8), intent(out) :: U(:,:)
      integer::i

      do i=is,ie
          U(i,IDN) = W(i,IDN)
          U(i,IM1) = W(i,IDN)*W(i,IV1)
          U(i,IM2) = W(i,IDN)*W(i,IV2)
          U(i,IM3) = W(i,IDN)*W(i,IV3)
          U(i,IEN) = 0.5d0*W(i,IDN)*( W(i,IV1)**2 + W(i,IV2)**2 + W(i,IV3)**2 ) &
                   + 0.5d0*( B(i,1)**2 + B(i,2)**2 + B(i,3)**2 ) &
                   + W(i,IPR)/(gam - 1.0d0)
      enddo
      
      return
      end subroutine Consvvariable

      subroutine PrimVariable( U, B, W )
      use commons
      use eosmod
      implicit none
      real(8), intent(in) :: U(:,:), B(:,:)
      real(8), intent(out) :: W(:,:)
      integer::i
      real(8) :: inv_d;

      do i=is,ie
           W(i,IDN) = U(i,IDN)
           inv_d = 1.0d0/U(i,IDN)
           W(i,IV1) = U(i,IM1)*inv_d
           W(i,IV2) = U(i,IM2)*inv_d
           W(i,IV3) = U(i,IM3)*inv_d
           W(i,IPR) = ( U(i,IEN) &
                    - 0.5d0*(U(i,IM1)**2 + U(i,IM2)**2 + U(i,IM3)**2)*inv_d  &
                    - 0.5d0*(B(i,1)**2 + B(i,2)**2 + B(i,3)**2) )*(gam-1.0d0)
      enddo

      return
      end subroutine PrimVariable

      subroutine TimestepControl(x1a, W,B)
      use commons
      use eosmod
      implicit none
      real(8), intent(in) :: x1a(:), W(:,:), B(:,:)
      real(8)::dtl1
      real(8)::dtl2
      real(8)::dtl3
      real(8)::dtlocal
      real(8)::dtmin,cf
      integer::i

      dtmin=1.0d90

      do i=is,ie
         cf = dsqrt( (gam*W(i,IPR) + B(i,1)**2 + B(i,2)**2 + B(i,3)**2)/W(i,IDN))
         dtlocal =(x1a(i+1)-x1a(i))/(abs(W(i,IV1)) + cf)
         if(dtlocal .lt. dtmin) dtmin = dtlocal
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
      subroutine NumericalFlux( W, B, flux )
      use commons !, only: is, ie, in
      implicit none
      integer::i
      real(8), intent(in) :: W(:,:), B(:,:)
      real(8), intent(out) :: flux(:,:)
      real(8),dimension(in,NFLX):: Wl,Wr
      real(8),dimension(NFLX):: flx
      real(8) :: dWm(NFLX), dWp(NFLX), dWmon(NFLX)
      real(8) :: ddmon, dvmon, dpmon

      do i=is-1,ie+1
         dWp(1:NVAR) = W(i+1,1:NVAR) - W(i  ,1:NVAR)
         dWm(1:NVAR) = W(i  ,1:NVAR) - W(i-1,1:NVAR)

         dWp(IB2:IB3) = B(i+1,2:3) - B(i,  2:3)
         dWm(IB2:IB3) = B(i  ,2:3) - B(i-1,2:3)

         call vanLeer(NFLX, dWp, dWm, dWmon)

         Wl(i+1,1:NVAR) = W(i,1:NVAR) + 0.5d0*dWmon(1:NVAR)*1.0d0
         Wr(i  ,1:NVAR) = W(i,1:NVAR) - 0.5d0*dWmon(1:NVAR)*1.0d0

         Wl(i+1,IB2:IB3) = B(i,2:3) + 0.5d0*dWmon(IB2:IB3)*1.0d0
         Wr(i  ,IB2:IB3) = B(i,2:3) - 0.5d0*dWmon(IB2:IB3)*1.0d0
      enddo

      do i=is,ie+1
         call HLL(Wl(i,:),Wr(i,:),0.5d0*(B(i-1,1) + B(i,1)),flx)
         flux(i,:)  = flx(:)
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
!            1D array (IDN, IV1, IV2, IV3, IPR, IB2, IB3)
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
!            index: (IDN, IV1, IV2, IV3, IPR, IB2, IB3)
!---------------------------------------------------------------------
      subroutine HLL(Wl,Wr,b1,flx)
      use commons !, only : is, ie, NVAR
      use eosmod
      implicit none
      real(8),intent(in)::Wl(:), Wr(:)
      real(8),intent(in) :: b1
      real(8),intent(out) :: flx(:)
      real(8):: Ul(NFLX), Ur(NFLX)
      real(8):: Fl(NFLX), Fr(NFLX)
      real(8):: Ust(NFLX)
      real(8):: Fst(NFLX)
      real(8):: cfl,cfr
      real(8):: sl, sr
      real(8):: pbl, pbr, ptotl, ptotr
      integer :: i, n

          pbl = 0.5d0*(b1**2 + Wl(IB2)**2 + Wl(IB3)**2)
          pbr = 0.5d0*(b1**2 + Wr(IB2)**2 + Wr(IB3)**2)
          ptotl = Wl(IPR) + pbl
          ptotr = Wr(IPR) + pbr

          ! conserved variables in the left and right states
          Ul(IDN) = Wl(IDN)
          Ul(IM1) = Wl(IDN)*Wl(IV1)
          Ul(IM2) = Wl(IDN)*Wl(IV2)
          Ul(IM3) = Wl(IDN)*Wl(IV3)
          Ul(IEN) = 0.5d0*Wl(IDN)*( Wl(IV1)**2 + Wl(IV2)**2 + Wl(IV3)**2) & 
                  + pbl + Wl(IPR)/(gam - 1.0d0)
          Ul(IB2) = Wl(IB2)
          Ul(IB3) = Wl(IB3)

          Ur(IDN) = Wr(IDN)
          Ur(IM1) = Wr(IDN)*Wr(IV1)
          Ur(IM2) = Wr(IDN)*Wr(IV2)
          Ur(IM3) = Wr(IDN)*Wr(IV3)
          Ur(IEN) = 0.5d0*Wr(IDN)*( Wr(IV1)**2 + Wr(IV2)**2 + Wr(IV3)**2) & 
                  + pbr + Wr(IPR)/(gam - 1.0d0)
          Ur(IB2) = Wr(IB2)
          Ur(IB3) = Wr(IB3)

          ! flux in the left and right states
          Fl(IDN) = Ul(IM1)
          Fl(IM1) = Wl(IDN)*Wl(IV1)**2 + ptotl - b1**2
          Fl(IM2) = Wl(IDN)*Wl(IV1)*Wl(IV2) - b1*Wl(IB2)
          Fl(IM3) = Wl(IDN)*Wl(IV1)*Wl(IV3) - b1*Wl(IB3)
          Fl(IEN) = ( Ul(IEN) + ptotl )*Wl(IV1) &
                  - b1*( b1*Wl(IV1) + Wl(IB2)*Wl(IV2) + Wl(IB3)*Wl(IV3) )
          Fl(IB2) = Wl(IB2)*Wl(IV1) - b1*Wl(IV2)
          Fl(IB3) = Wl(IB3)*Wl(IV1) - b1*Wl(IV3)

          Fr(IDN) = Ur(IM1)
          Fr(IM1) = Wr(IDN)*Wr(IV1)**2 + ptotr - b1**2
          Fr(IM2) = Wr(IDN)*Wr(IV1)*Wr(IV2) - b1*Wr(IB2)
          Fr(IM3) = Wr(IDN)*Wr(IV1)*Wr(IV3) - b1*Wr(IB3)
          Fr(IEN) = ( Ur(IEN) + ptotr )*Wr(IV1) &
                  - b1*( b1*Wr(IV1) + Wr(IB2)*Wr(IV2) + Wr(IB3)*Wr(IV3) )
          Fr(IB2) = Wr(IB2)*Wr(IV1) - b1*Wr(IV2)
          Fr(IB3) = Wr(IB3)*Wr(IV1) - b1*Wr(IV3)

          cfl = dsqrt( (gam*Wl(IPR) + Wl(IB2)**2 + Wl(IB3)**2 + b1**2)/Wl(IDN))
          cfr = dsqrt( (gam*Wr(IPR) + Wr(IB2)**2 + Wr(IB3)**2 + b1**2)/Wr(IDN))

!          sl = min(Wl(IV1),Wr(IV1)) - max(cfl,cfr)
!          sr = max(Wl(IV1),Wr(IV1)) + max(cfl,cfr)
          sl = min(Wl(IV1) - cfl,Wr(IV1) - cfr)
          sr = max(Wl(IV1) + cfl,Wr(IV1) + cfr)
!          Fst(:)  = (sr*Fl(:) - sl*Fr(:) + sl*sr*( Ur(:) - Ul(:) ))/(sr - sl)
!          Ust(:) = ( sr*Ur(:) - sl*Ul(:) - Fr(:) + Fl(:) )/(sr - sl)


          if( sl > 0.0d0 ) then
               flx(:)  = Fl(:) 
          else if (sr <= 0.0d0 ) then
               flx(:)  = Fr(:) 
          else 
               flx(:)  = (sr*Fl(:) - sl*Fr(:) + sl*sr*( Ur(:) - Ul(:) ))/(sr - sl)
          endif

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

      subroutine UpdateConsv( flux, U, B )
      use commons
      implicit none
      real(8), intent(in)  :: flux(:,:)
      real(8), intent(out) :: U(:,:), B(:,:)
      integer::i,n

      do n=1,NVAR
      do i=is,ie
         U(i,n) = U(i,n) + dt*(- flux(i+1,n) + flux(i,n))/(x1a(i+1)-x1a(i)) 
      enddo
      enddo

      do i=is,ie
         B(i,2) = B(i,2) + dt*(- flux(i+1,IB2) + flux(i,IB2))/(x1a(i+1)-x1a(i)) 
         B(i,3) = B(i,3) + dt*(- flux(i+1,IB3) + flux(i,IB3))/(x1a(i+1)-x1a(i)) 
      enddo
      

      return
      end subroutine UpdateConsv

      subroutine Output( x1a, x1b, W, B )
      use commons
      implicit none
      real(8), intent(in) :: x1a(:), x1b(:), W(:,:), B(:,:)
      integer::i
      character(20),parameter::dirname="bindata/"
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
         call makedirs("bindata")
         is_inited =.true.
      endif

!      print*, time, tout+dtout, time+1.0d-14 .lt. tout+dtout
      if(time + 1.0d-14.lt. tout+dtout) return

      write(filename,'(a3,i5.5,a4)')"bin",nout,".dat"
      filename = trim(dirname)//filename
!      open(unitbin,file=filename,status='replace',form='formatted') 
      open(unitbin,file=filename,form='formatted',action="write")
      write(unitbin,*) "# time = ",time
      do i=1,in-1
          write(unitbin,*) x1b(i), W(i,IDN), W(i,IV1), W(i,IV2), W(i,IV3), W(i,IPR), B(i,1), &
          B(i,2), B(i,3)
!          write(*,*) x1b(i), d(i), v(i), p(i)
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
