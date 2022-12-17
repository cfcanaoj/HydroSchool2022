      module commons
      implicit none
      integer::ntime                        ! counter of the timestep
      integer,parameter::ntimemax=20000     ! the maximum timesteps
      real(8)::time,dt                    ! time, timewidth
      data time / 0.0d0 /
      real(8),parameter:: timemax=0.2d0   
      real(8),parameter:: dtout=5.0d-3

      integer,parameter::ngrid=150        ! the number of grids in the simulation box
      integer,parameter::mgn=2            ! the number of ghost cells
      integer,parameter::in=ngrid+2*mgn+1 ! the total number of grids including ghost cells
      integer,parameter::is=mgn+1         ! the index of the leftmost grid
      integer,parameter::ie=ngrid+mgn     ! the index of the rightmost grid
      integer,parameter::nflux = 3

      real(8),parameter:: x1min=-0.5d0,x1max=0.5d0

      end module commons
     
      module eosmod
      implicit none
! adiabatic
!      real(8),parameter::gam=5.0d0/3.0d0 !! adiabatic index
      real(8),parameter::gam=1.4d0 !! adiabatic index
! isothermal
!      real(8)::csiso  !! isothemal sound speed
      end module eosmod

      program main
      use commons
      implicit none

      real(8),dimension(in)::x1a,x1b
! arrays of the conserved variables
      real(8),dimension(in) :: d,et
      real(8),dimension(in,3) :: mv
! arrays of the primitive variables
      real(8),dimension(in) :: p
      real(8),dimension(in,3) :: v

      real(8),dimension(in,3) :: B

      real(8),dimension(in) :: dflux, etflux
      real(8),dimension(in,3) :: mvflux, bflux


      write(6,*) "setup grids and initial condition"
      call GenerateGrid(x1a, x1b)
      call GenerateProblem(x1b, d, v, p, B)
      call ConsvVariable(d, v, p, B, et, mv)
      call Output( x1a, x1b, d, v, p, B )
! main loop
                                  write(6,*)"step","time","dt"
      mloop: do ntime=1,ntimemax
         call TimestepControl(x1a, d, v, p, B)
         if( time + dt > timemax ) dt = timemax - time
!         if(mod(ntime,100) .eq. 0 ) write(6,*)ntime,time,dt
         call BoundaryCondition(d,v,p,B)
         call NumericalFlux(d, v, p, B, dflux, mvflux, etflux, Bflux)
         call UpdateConsv( dflux, mvflux, etflux, d, mv, et )
         call PrimVariable( d, mv, et, v, p )
         time=time+dt
         call Output( x1a, x1b, d, v, p )

         if(time >= timemax) exit mloop
      enddo mloop
      call Output( x1a, x1b, d, v, p )

      write(6,*) "program has been finished"
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

      subroutine GenerateProblem(x1b, d, v, p, B)
      use commons
      use eosmod
      implicit none
      integer::i
      real(8), intent(in ) :: x1b(:)
      real(8), intent(out) :: d(:), B(:,:), v(:,:), p(:)
      real(8) :: rho1,rho2,Lsm,u1,u2
      data rho1  /  1.0d0 /
      data rho2  /  2.0d0 /
      data u1    /  0.5d0 /
      data u2    / -0.5d0 /
      data Lsm   /  0.025d0 /

      real(8)::pi
      pi=acos(-1.0d0)

      do i=is,ie
         if( x1b(i) < 0.0d0 ) then 
             d(i) = 1.0d0
             v(i,1) = 0.0d0
             v(i,2) = 0.0d0
             v(i,3) = 0.0d0
             p(i) = 1.0d0

             B(i,1) = 0.75d0
             B(i,2) = 1.0d0
             B(i,3) = 0.0d0
         else 
             d(i) = 0.125d0
             v(i,1) = 0.0d0
             v(i,2) = 0.0d0
             v(i,3) = 0.0d0
             p(i) = 0.1d0

             B(i,1) = 0.75d0
             B(i,2) = -1.0d0
             B(i,3) = 0.0d0
         endif
      enddo

      
!      call BoundaryCondition

      return
      end subroutine GenerateProblem

      subroutine BoundaryCondition(d,v,p,B)
      use commons
      implicit none
      real(8), intent(out) :: d(:), v(:,:), p(:), B(:,:)
      integer::i

      do i=1,mgn
!          d(i) =  d(ie-mgn+i)
!          v(i) = v(ie-mgn+i)
!          p(i) = p(ie-mgn+i)
          d(is-i)  = d(is-1+i)
          v(is-i,1)  = v(is-1+i,1)
          v(is-i,2)  = v(is-1+i,2)
          v(is-i,3)  = v(is-1+i,3)
          p(is-i)  = p(is-1+i)
          B(is-i,1)  = B(is-1+i,1)
          B(is-i,2)  = B(is-1+i,2)
          B(is-i,3)  = B(is-1+i,3)
      enddo

      do i=1,mgn
!           d(ie+i) =  d(is+i-1)
!          v(ie+i) = v(is+i-1)
!          p(ie+i) = p(is+i-1)

          d(ie+i) =  d(ie-i+1)
          v(ie+i,1) = v(ie-i+1,1)
          v(ie+i,2) = v(ie-i+1,2)
          v(ie+i,3) = v(ie-i+1,3)
          p(ie+i) = p(ie-i+1)
          B(ie+i,1) = B(ie-i+1,1)
          B(ie+i,2) = B(ie-i+1,2)
          B(ie+i,3) = B(ie-i+1,3)
      enddo

      return
      end subroutine BoundaryCondition
!
      subroutine ConsvVariable(d, v, p, B, et, mv)
      use commons
      use eosmod
      implicit none
      real(8), intent(in) :: d(:), v(:,:), p(:), B(:,:)
      real(8), intent(out) :: et(:), mv(:,:)
      integer::i

      do i=is,ie
          et(i)  = 0.5d0*d(i)*( v(i,1)**2 + v(i,2)**2 + v(i,3)**2 )  &
                 + 0.5d0*( B(i,1)**2 + B(i,2)**2 + B(i,3)**2 ) &
                 + p(i)/(gam - 1.0d0)
          mv(i,1) = d(i)*v(i,1)
          mv(i,2) = d(i)*v(i,2)
          mv(i,3) = d(i)*v(i,3)
      enddo
      
      return
      end subroutine Consvvariable

      subroutine PrimVariable( d, mv, et, v, p )
      use commons
      use eosmod
      implicit none
      real(8), intent(in) :: d(:), mv(:), et(:)
      real(8), intent(out) :: v(:), p(:)
      integer::i

      do i=is,ie
           v(i) = mv(i)/d(i)

!! adiabatic
           p(i) =  ( et(i) - 0.5d0*d(i)*v(i)**2 )*(gam-1.0d0)
!          cs(i) =  dsqrt(gam*p(i)/d(i))
! isotermal
!           p(i) =  d(i)*csiso**2
!          cs(i) =  csiso
      enddo

      return
      end subroutine PrimVariable

      subroutine TimestepControl(x1a, d, v, p, B)
      use commons
      use eosmod
      implicit none
      real(8), intent(in) :: x1a(:), d(:), v(:,:), p(:), B(:,:)
      real(8)::dtl1
      real(8)::dtl2
      real(8)::dtl3
      real(8)::dtlocal
      real(8)::dtmin, cf
      integer::i

      dtmin=1.0d90

      do i=is,ie
         cf = dsqrt(( gam*p(i) + B(i,1)**2 + B(i,2)**2 + B(i,3)**2 )/d(i))
         dtlocal =(x1a(i+1)-x1a(i))/(abs(v(i,1)) + cf )
         if(dtlocal .lt. dtmin) dtmin = dtlocal
      enddo

      dt = 0.05d0 * dtmin
!      write(6,*)"dt",dt
      return
      end subroutine TimestepControl

!!      subroutine StateVevtor
!!      use commons
!!      use fluxmod
!!      use eosmod
!!      implicit none
!!      integer::i
!!
!!      do i=1,in-1
!!         svc(nden,i) =  d(i)
!!         svc(nve1,i) = v(i)
!!! adiabatic
!!         svc(nene,i) = ei(i)/d(i)
!!         svc(npre,i) = ei(i)*(gam-1.0d0)
!!         svc(ncsp,i) = sqrt(gam*(gam-1.0d0)*ei(i)/d(i))
!!! isotermal
!!!         svc(nene,i) = csiso**2
!!!         svc(npre,i) = d(i)*csiso**2
!!!         svc(ncsp,i) = csiso
!!         p(i) = svc(npre,i)  ! for output boundary 
!!      enddo
!!
!!      return
!!      end subroutine StateVevtor
!
!      subroutine minmod(a,b,d)
!      use fluxmod, only : ntimed
!      implicit none
!      real(8),dimension(ntimed),intent(in)::a,b
!      real(8),dimension(ntimed),intent(out)::d
!      integer:: n
!
!      do n=1,ntimed
!         d(n) = sign(1.0d0,a(n))*max(0.0d0,min(abs(a(n))                &
!     &                                        ,sign(1.0d0,a(n))*b(n)))
!      enddo
!
!      return
!      end subroutine minmod
!
!
      subroutine vanLeer(dvp,dvm,dv)
      implicit none
      real(8),intent(in)::dvp,dvm
      real(8),intent(out)::dv
      integer:: n

         if(dvp*dvm .gt. 0.0d0)then
            dv =2.0d0*dvp*dvm/(dvp+dvm)
         else
            dv = 0.0d0
         endif

      return
      end subroutine vanLeer
!
!
!
!      subroutine MClimiter(a,b,c,d)
!      use fluxmod, only : ntimed
!      implicit none
!      real(8),dimension(ntimed),intent(in)::a,b,c
!      real(8),dimension(ntimed),intent(out)::d
!      integer:: n
!
!      do n=1,ntimed
!         d(n) = sign(1.0d0,a(n))*max(0.0d0,min(abs(a(n))         &
!     &                                  ,sign(1.0d0,a(n))*b(n)   &
!     &                                  ,sign(1.0d0,a(n))*c(n))) 
!      enddo
!
!      return
!      end subroutine MClimiter
!
      subroutine NumericalFlux(d, v, p, B, dflux, mvflux, etflux, Bflux)
      use commons, only: is, ie, in
      implicit none
      integer::i
      real(8), intent(in) :: d(:), v(:,:), p(:), B(:,:)
      real(8), intent(out) :: dflux(:), mvflux(:,:), etflux(:), Bflux(:,:)
      real(8),dimension(in):: dl,dr
      real(8),dimension(in,8):: Wl,Wr
      real(8),dimension(8):: flux
      real(8) :: dWm(8), dWp(8), dWmon(8)

      do i=is-1,ie+1
         dWp(1) = d(i+1) - d(i)
         dWm(1) = d(i) - d(i-1)

         dWp(2) = v(i+1,1) - v(i,1)
         dWm(2) = v(i,1) - v(i-1,1)

         dWp(3) = v(i+1,2) - v(i,2)
         dWm(3) = v(i,2) - v(i-1,2)

         dWp(4) = v(i+1,3) - v(i,3)
         dWm(4) = v(i,3) - v(i-1,3)

         dWp(5) = p(i+1,2) - p(i,2)
         dWm(5) = p(i,2) - p(i-1,2)

         dWp(6) = B(i+1,1) - B(i,1)
         dWm(6) = B(i,1) - B(i-1,1)

         dWp(7) = B(i+1,2) - B(i,2)
         dWm(7) = B(i,2) - B(i-1,2)

         dWp(8) = B(i+1,3) - B(i,3)
         dWm(8) = B(i,3) - B(i-1,3)

         call vanLeer(dWp,dWm,dWmon)

         dl(i+1) = d(i) + 0.5d0*ddmon*0.0d0
         dr(i  ) = d(i) - 0.5d0*ddmon*0.0d0

         vl(i+1,1) = v(i,1) + 0.5d0*dvmon1*0.0d0
         vr(i  ,1) = v(i,1) - 0.5d0*dvmon1*0.0d0

         vl(i+1,2) = v(i,2) + 0.5d0*dvmon2*0.0d0
         vr(i  ,2) = v(i,2) - 0.5d0*dvmon2*0.0d0

         vl(i+1,3) = v(i,3) + 0.5d0*dvmon3*0.0d0
         vr(i  ,3) = v(i,3) - 0.5d0*dvmon3*0.0d0

         pl(i+1) = p(i) + 0.5d0*dpmon*0.0d0
         pr(i  ) = p(i) - 0.5d0*dpmon*0.0d0

         Bl(i+1,1) = B(i,1) + 0.5d0*dbmon1*0.0d0
         Br(i  ,1) = B(i,1) - 0.5d0*dbmon1*0.0d0

         Bl(i+1,2) = B(i,2) + 0.5d0*dbmon2*0.0d0
         Br(i  ,2) = B(i,2) - 0.5d0*dbmon2*0.0d0

         Bl(i+1,3) = B(i,3) + 0.5d0*dbmon3*0.0d0
         Br(i  ,3) = B(i,3) - 0.5d0*dbmon3*0.0d0
      enddo

!         call HLLE(leftst,rigtst,nflux)
      do i=is,ie+1
         call HLLE(dl(i),vl(i,:),pl(i),dr(i),vr(i),pr(i),flux)
         dflux(i)  = flux(1)
         mvflux(i) = flux(2)
         etflux(i) = flux(3)
      enddo

      return
      end subroutine Numericalflux

!      subroutine HLLE(leftst,rigtst,nflux)
      subroutine HLLE(dl,vl,pl,dr,vr,pr,flux)
      use commons, only : is, ie
      use eosmod
      implicit none
      real(8),intent(in)::dl,vl,pl,dr,vr,pr
      real(8),intent(out) :: flux(:)
      real(8):: csl,csr
      real(8):: sl, sr
      real(8):: dfluxl, dfluxr
      real(8):: mvl, mvr, mvfluxl, mvfluxr
      real(8):: etl, etr, etfluxl, etfluxr
      integer :: i

          ! conserved variables in the left and right states
          mvl = dl*vl
          mvr = dr*vr
          etl = pl/(gam - 1.0d0) + 0.5d0*dl*vl**2
          etr = pr/(gam - 1.0d0) + 0.5d0*dr*vr**2

          ! flux in the left and right states
          dfluxl = mvl
          dfluxr = mvr

          mvfluxl = pl + dl*vl**2
          mvfluxr = pr + dr*vr**2

          etfluxl = ( gam*pl/(gam - 1.0d0) + 0.5d0*dl*vl**2 )*vl
          etfluxr = ( gam*pr/(gam - 1.0d0) + 0.5d0*dr*vr**2 )*vr

          csl = dsqrt(gam*pl/dl)
          csr = dsqrt(gam*pr/dr)

          sl = min(vl,vr) - max(csl,csr)
!          sl = min(0.0d0,sl)
          sr = max(vl,vr) + max(csl,csr)
!          sr = max(0.0d0,sr)

          flux(1)  = (sr*dfluxl - sl*dfluxr + sl*sr*( dr - dl ))/(sr - sl)
          flux(2) = (sr*mvfluxl - sl*mvfluxr + sl*sr*( mvr - mvl ))/(sr - sl)
          flux(3) = (sr*etfluxl - sl*etfluxr + sl*sr*( etr - etl ))/(sr - sl)

      return
      end subroutine HLLE
!
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

      subroutine UpdateConsv( dflux, mvflux, etflux, d, mv, et )
      use commons
      implicit none
      real(8), intent(in)  :: dflux(:), mvflux(:), etflux(:)
      real(8), intent(out) :: d(:), mv(:), et(:)
      integer::i

      do i=is,ie
         d(i) = d(i) & 
     & +dt*(                                       &
     &  (- dflux(i+1)                    &
     &   + dflux(i  ))/(x1a(i+1)-x1a(i)) &
     &  )

         mv(i) = mv(i)                   &
     & +dt*(                                       &
     &  (- mvflux(i+1)                    &
     &   + mvflux(i  ))/(x1a(i+1)-x1a(i)) &
     &  )

          et(i) = et(i)                    &
     & +dt*(                                       &
     &  (- etflux(i+1)                    &
     &   + etflux(i  ))/(x1a(i+1)-x1a(i)) &
     &      )
      enddo

      return
      end subroutine UpdateConsv

      subroutine Output( x1a, x1b, d, v, p, B )
      use commons
      implicit none
      real(8), intent(in) :: x1a(:), x1b(:), d(:), v(:,:), p(:), B(:,:)
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
      integer,parameter:: nvar=3
      real(8)::x1out(is-gs:ie+gs,2)
      real(8)::hydout(is-gs:ie+gs,nvar)

      logical, save:: is_inited
      data is_inited /.false./

      if (.not. is_inited) then
         call makedirs("bindata")
         is_inited =.true.
      endif

      print*, time, tout+dtout, time+1.0d-14 .lt. tout+dtout
      if(time + 1.0d-14.lt. tout+dtout) return

      write(filename,'(a3,i5.5,a4)')"bin",nout,".dat"
      filename = trim(dirname)//filename
!      open(unitbin,file=filename,status='replace',form='formatted') 
      open(unitbin,file=filename,form='formatted',action="write")
      write(unitbin,*) "# time = ",time
      do i=1,in
          write(unitbin,*) x1b(i), d(i), v(i,1), v(i,2), v(i,3), p(i), &
                            B(i,1), B(i,2), B(i,3)
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
      write(*, *) trim(command)
      call system(command)
      end subroutine makedirs
end program main
