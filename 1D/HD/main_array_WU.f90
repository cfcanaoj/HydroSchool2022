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

      integer, parameter :: IDN = 1
      integer, parameter :: IM1 = 2
      integer, parameter :: IPR = 3
      integer, parameter :: NVAR = 3

      integer, parameter :: IV1 = 2
      integer, parameter :: IEN = 3

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
      real(8),dimension(in,NVAR) :: U
      real(8),dimension(in,NVAR) :: W
      real(8),dimension(in,NVAR) :: flux

      write(6,*) "setup grids and initial condition"
      call GenerateGrid(x1a, x1b)
      call GenerateProblem(x1b, W)
      call ConsvVariable(W, U)
      call Output( x1a, x1b, W )

! main loop
      mloop: do ntime=1,ntimemax
         call TimestepControl(x1a, W)
         if( time + dt > timemax ) dt = timemax - time
         call BoundaryCondition(W)
         call NumericalFlux(W, flux)
         call UpdateConsv( flux, U )
         call PrimVariable( U, W )
         time=time+dt
         call Output( x1a, x1b, W)

         if(time >= timemax) exit mloop
      enddo mloop
      call Output( x1a, x1b, W)

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

      subroutine GenerateProblem(x1b, W)
      use commons
      use eosmod
      implicit none
      integer::i
      real(8), intent(in ) :: x1b(:)
      real(8), intent(out) :: W(:,:)
      real(8) :: rho1,rho2,Lsm,u1,u2

      real(8)::pi
      pi=acos(-1.0d0)

      do i=is,ie
         if( x1b(i) < 0.0d0 ) then 
             W(i,IDN) = 1.0d0
             W(i,IV1) = 0.0d0
             W(i,IPR) = 1.0d0
         else 
             W(i,IDN) = 0.125d0
             W(i,IV1) = 0.0d0
             W(i,IPR) = 0.1d0
         endif
      enddo

      
!      call BoundaryCondition

      return
      end subroutine GenerateProblem

      subroutine BoundaryCondition(W)
      use commons
      implicit none
      real(8), intent(inout) :: W(:,:)
      integer::i

      do i=1,mgn
!          d(i) =  d(ie-mgn+i)
!          v(i) = v(ie-mgn+i)
!          p(i) = p(ie-mgn+i)
          W(is-i,IDN)  = W(is-1+i,IDN)
          W(is-i,IV1)  = W(is-1+i,IV1)
          W(is-i,IPR)  = W(is-1+i,IPR)
      enddo

      do i=1,mgn
!           d(ie+i) =  d(is+i-1)
!          v(ie+i) = v(is+i-1)
!          p(ie+i) = p(is+i-1)

          W(ie+i,IDN) = W(ie-i+1,IDN)
          W(ie+i,IV1) = W(ie-i+1,IV1)
          W(ie+i,IPR) = W(ie-i+1,IPR)
      enddo

      return
      end subroutine BoundaryCondition
!
      subroutine ConsvVariable(W, U)
      use commons
      use eosmod
      implicit none
      real(8), intent(in) :: W(:,:)
      real(8), intent(out) :: U(:,:)
      integer::i

      do i=is,ie
          U(i,IDN) = W(i,IDN)
          U(i,IM1) = W(i,IDN)*W(i,IV1)
          U(i,IEN)  = 0.5d0*W(i,IDN)*W(i,IV1)**2 + W(i,IPR)/(gam - 1.0d0)
      enddo
      
      return
      end subroutine Consvvariable

      subroutine PrimVariable( U, W )
      use commons
      use eosmod
      implicit none
      real(8), intent(in) :: U(:,:)
      real(8), intent(out) :: W(:,:)
      integer::i

      do i=is,ie
           W(i,IDN) = U(i,IDN)
           W(i,IV1) = U(i,IM1)/U(i,IDN)
           W(i,IPR) = ( U(i,IEN) - 0.5d0*U(i,IM1)**2/U(i,IDN) )*(gam-1.0d0)
      enddo

      return
      end subroutine PrimVariable

      subroutine TimestepControl(x1a, W)
      use commons
      use eosmod
      implicit none
      real(8), intent(in) :: x1a(:), W(:,:)
      real(8)::dtl1
      real(8)::dtl2
      real(8)::dtl3
      real(8)::dtlocal
      real(8)::dtmin
      integer::i

      dtmin=1.0d90

      do i=is,ie
         dtlocal =(x1a(i+1)-x1a(i))/(abs(W(i,IV1)) + dsqrt(gam*W(i,IPR)/W(i,IDN)))
         if(dtlocal .lt. dtmin) dtmin = dtlocal
      enddo

      dt = 0.3d0 * dtmin
!      write(6,*)"dt",dt
      return
      end subroutine TimestepControl

      subroutine vanLeer(n,dvp,dvm,dv)
      implicit none
      real(8),intent(in)::dvp(:),dvm(:)
      real(8),intent(out)::dv(:)
      integer,intent(in) :: n
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
      subroutine NumericalFlux( W, flux )
      use commons !, only: is, ie, in
      implicit none
      integer::i
      real(8), intent(in) :: W(:,:)
      real(8), intent(out) :: flux(:,:)
      real(8),dimension(in,NVAR):: Wl,Wr
      real(8),dimension(NVAR):: flx
      real(8) :: dWm(NVAR), dWp(NVAR), dWmon(NVAR)
      real(8) :: ddmon, dvmon, dpmon

      do i=is-1,ie+1
         dWp(:) = W(i+1,:) - W(i  ,:)
         dWm(:) = W(i  ,:) - W(i-1,:)

         call vanLeer(NVAR, dWp,dWm,dWmon)

         Wl(i+1,:) = W(i,:) + 0.5d0*dWmon(:)*1.0d0
         Wr(i  ,:) = W(i,:) - 0.5d0*dWmon(:)*1.0d0
      enddo

!         call HLLE(leftst,rigtst,nflux)
      do i=is,ie+1
         call HLLE(Wl(i,:),Wr(i,:),flx)
         flux(i,:)  = flx(:)
      enddo

      return
      end subroutine Numericalflux

!      subroutine HLLE(leftst,rigtst,nflux)
      subroutine HLLE(Wl,Wr,flx)
      use commons !, only : is, ie, NVAR
      use eosmod
      implicit none
      real(8),intent(in)::Wl(:), Wr(:)
      real(8),intent(out) :: flx(:)
      real(8):: Ul(NVAR), Ur(NVAR)
      real(8):: Fl(NVAR), Fr(NVAR)
      real(8):: csl,csr
      real(8):: sl, sr
      real(8):: dfluxl, dfluxr
      real(8):: mvl, mvr, mvfluxl, mvfluxr
      real(8):: etl, etr, etfluxl, etfluxr
      integer :: i, n

          ! conserved variables in the left and right states
          Ul(IDN) = Wl(IDN)
          Ur(IDN) = Wr(IDN)

          Ul(IM1) = Wl(IDN)*Wl(IV1)
          Ur(IM1) = Wr(IDN)*Wr(IV1)

          Ul(IEN) = 0.5d0*Wl(IDN)*Wl(IV1)**2 + Wl(IPR)/(gam - 1.0d0)
          Ur(IEN) = 0.5d0*Wr(IDN)*Wr(IV1)**2 + Wr(IPR)/(gam - 1.0d0)

          ! flux in the left and right states
          Fl(IDN) = Ul(IM1)
          Fr(IDN) = Ur(IM1)

          Fl(IM1) = Wl(IPR) + Wl(IDN)*Wl(IV1)**2 
          Fr(IM1) = Wr(IPR) + Wr(IDN)*Wr(IV1)**2 


          Fl(IEN) = ( gam*Wl(IPR)/(gam - 1.0d0) + 0.5d0*Wl(IDN)*Wl(IV1)**2)*Wl(IV1)
          Fr(IEN) = ( gam*Wr(IPR)/(gam - 1.0d0) + 0.5d0*Wr(IDN)*Wr(IV1)**2)*Wr(IV1)

          csl = dsqrt(gam*Wl(IPR)/Wl(IDN))
          csr = dsqrt(gam*Wr(IPR)/Wr(IDN))

          sl = min(Wl(IV1),Wr(IV1)) - max(csl,csr)
          sr = max(Wl(IV1),Wr(IV1)) + max(csl,csr)

          do n=1,NVAR
             flx(n)  = (sr*Fl(n) - sl*Fr(n) + sl*sr*( Ur(n) - Ul(n) ))/(sr - sl)
          enddo

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

      subroutine UpdateConsv( flux, U )
      use commons
      implicit none
      real(8), intent(in)  :: flux(:,:)
      real(8), intent(out) :: U(:,:)
      integer::i,n

      do n=1,NVAR
      do i=is,ie
         U(i,n) = U(i,n) + dt*(- flux(i+1,n) + flux(i,n))/(x1a(i+1)-x1a(i)) 
      enddo
      enddo

      return
      end subroutine UpdateConsv

      subroutine Output( x1a, x1b, W )
      use commons
      implicit none
      real(8), intent(in) :: x1a(:), x1b(:), W(:,:)
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

      print*, time, tout+dtout, time+1.0d-14 .lt. tout+dtout
      if(time + 1.0d-14.lt. tout+dtout) return

      write(filename,'(a3,i5.5,a4)')"bin",nout,".dat"
      filename = trim(dirname)//filename
!      open(unitbin,file=filename,status='replace',form='formatted') 
      open(unitbin,file=filename,form='formatted',action="write")
      write(unitbin,*) "# time = ",time
      do i=1,in-1
          write(unitbin,*) x1b(i), W(i,IDN), W(i,IV1), W(i,IPR)
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
