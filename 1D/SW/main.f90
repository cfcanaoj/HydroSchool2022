      module commons
      implicit none
      integer::nhy                        ! counter of the timestep
      integer,parameter::nhymax=20000     ! the maximum timesteps
      real(8)::time,dt                    ! time, timewidth
      data time / 0.0d0 /
      real(8),parameter:: timemax=0.2d0   
      real(8),parameter:: dtout=5.0d-3

      integer,parameter::ngrid=150        ! the number of grids in the simulation box
      integer,parameter::mgn=2            ! the number of ghost cells
      integer,parameter::in=ngrid+2*mgn+1 ! the total number of grids including ghost cells
      integer,parameter::is=mgn+1         ! the index of the leftmost grid
      integer,parameter::ie=ngrid+mgn     ! the index of the rightmost grid

      real(8),parameter:: x1min=-0.5d0,x1max=0.5d0
      real(8),dimension(in)::x1a,x1b

! arrays of the conserved variables
      real(8),dimension(in)::d,et,mv
! arrays of the primitive variables
      real(8),dimension(in)::p,ei,v,cs

      end module commons
     
      module eosmod
      implicit none
! adiabatic
      real(8),parameter::gam=5.0d0/3.0d0 !! adiabatic index
! isothermal
!      real(8)::csiso  !! isothemal sound speed
end module eosmod

      module fluxmod
      use commons, only : in
      implicit none
      integer,parameter::nden=1,nve1=2,nene=3,npre=4,ncsp=5
      integer,parameter::nhyd=5
      real(8),dimension(nhyd,in):: svc

      integer,parameter::mudn=1,muvu=2,muvv=3,muvw=4,muet=5  &
     &                  ,mfdn=6,mfvu=7,mfvv=8,mfvw=9,mfet=10 &
     &                  ,mcsp=11,mvel=12,mpre=13
      integer,parameter:: mflx=5,madd=3

      integer,parameter:: mden=1,mrv=2,meto=3 &
     &                          ,mrvu=muvu,mrvv=muvv,mrvw=muvw
      real(8),dimension(in) :: dflux, mvflux, etflux

      end module fluxmod

      program main
      use commons
      implicit none
      write(6,*) "setup grids and initial condition"
      call GenerateGrid
      call GenerateProblem
      call ConsvVariable
      write(6,*) "entering main loop"
      call Output
! main loop
                                  write(6,*)"step","time","dt"
      mloop: do nhy=1,nhymax
         call TimestepControl
         if(mod(nhy,100) .eq. 0 ) write(6,*)nhy,time,dt
         call BoundaryCondition
!         call StateVevtor
         call NumericalFlux
         call UpdateConsv
         call PrimVariable
         time=time+dt
         call Output
         if(time > timemax) exit mloop
         print*,time
      enddo mloop

      write(6,*) "program has been finished"
      end program main

      subroutine GenerateGrid
      use commons
      use eosmod
      implicit none
      real(8)::dx,dy
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

      subroutine GenerateProblem
      use commons
      use eosmod
      implicit none
      integer::i
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
             v(i) = 0.0d0
             p(i) = 1.0d0
         else 
             d(i) = 0.125d0
             v(i) = 0.0d0
             p(i) = 0.1d0
         endif
      enddo

      do i=is,ie
          ei(i) = p(i)/(gam-1.0d0)
          cs(i) = sqrt(gam*p(i)/d(i))
      enddo
      
      call BoundaryCondition

      return
      end subroutine GenerateProblem

      subroutine BoundaryCondition
      use commons
      implicit none
      integer::i

      do i=1,mgn
!          d(i) =  d(ie-mgn+i)
!          ei(i) = ei(ie-mgn+i)
!          v(i) = v(ie-mgn+i)
!          p(i) = p(ie-mgn+i)
          d(is-i)  = d(is-1+i)
          ei(is-i) = ei(is-1+i)
          v(is-i)  = v(is-1+i)
          p(is-i)  = p(is-1+i)
      enddo

      do i=1,mgn
!           d(ie+i) =  d(is+i-1)
!          ei(ie+i) = ei(is+i-1)
!          v(ie+i) = v(is+i-1)
!          p(ie+i) = p(is+i-1)

          d(ie+i) =  d(ie-i+1)
          ei(ie+i) = ei(ie-i+1)
          v(ie+i) = v(ie-i+1)
          p(ie+i) = p(ie-i+1)
      enddo

      return
      end subroutine BoundaryCondition

      subroutine ConsvVariable
      use commons
      implicit none
      integer::i

      do i=is,ie
          et(i)  = 0.5d0*d(i)*v(i)**2 + ei(i)
          mv(i) = d(i)*v(i)
      enddo
      
      return
      end subroutine Consvvariable

      subroutine PrimVariable
      use commons
      use eosmod
      implicit none
      integer::i
      do i=is,ie
          v(i) = mv(i)/d(i)

          ei(i) = et(i) - 0.5d0*d(i)*v(i)**2  

! adiabatic
           p(i) =  ei(i)*(gam-1.0d0)
          cs(i) =  dsqrt(gam*p(i)/d(i))
! isotermal
!           p(i) =  d(i)*csiso**2
!          cs(i) =  csiso
      enddo

      return
      end subroutine PrimVariable

      subroutine TimestepControl
      use commons
      implicit none
      real(8)::dtl1
      real(8)::dtl2
      real(8)::dtl3
      real(8)::dtlocal
      real(8)::dtmin
      integer::i
      dtmin=1.0d90
      do i=is,ie
         dtlocal =(x1a(i+1)-x1a(i))/(abs(v(i)) +cs(i))
         if(dtlocal .lt. dtmin) dtmin = dtlocal
      enddo

      dt = 0.05d0 * dtmin
!      write(6,*)"dt",dt
      return
      end subroutine TimestepControl

!      subroutine StateVevtor
!      use commons
!      use fluxmod
!      use eosmod
!      implicit none
!      integer::i
!
!      do i=1,in-1
!         svc(nden,i) =  d(i)
!         svc(nve1,i) = v(i)
!! adiabatic
!         svc(nene,i) = ei(i)/d(i)
!         svc(npre,i) = ei(i)*(gam-1.0d0)
!         svc(ncsp,i) = sqrt(gam*(gam-1.0d0)*ei(i)/d(i))
!! isotermal
!!         svc(nene,i) = csiso**2
!!         svc(npre,i) = d(i)*csiso**2
!!         svc(ncsp,i) = csiso
!         p(i) = svc(npre,i)  ! for output boundary 
!      enddo
!
!      return
!      end subroutine StateVevtor

      subroutine minmod(a,b,d)
      use fluxmod, only : nhyd
      implicit none
      real(8),dimension(nhyd),intent(in)::a,b
      real(8),dimension(nhyd),intent(out)::d
      integer:: n

      do n=1,nhyd
         d(n) = sign(1.0d0,a(n))*max(0.0d0,min(abs(a(n))                &
     &                                        ,sign(1.0d0,a(n))*b(n)))
      enddo

      return
      end subroutine minmod


      subroutine vanLeer(dvp,dvm,dv)
      use fluxmod, only : nhyd
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



      subroutine MClimiter(a,b,c,d)
      use fluxmod, only : nhyd
      implicit none
      real(8),dimension(nhyd),intent(in)::a,b,c
      real(8),dimension(nhyd),intent(out)::d
      integer:: n

      do n=1,nhyd
         d(n) = sign(1.0d0,a(n))*max(0.0d0,min(abs(a(n))         &
     &                                  ,sign(1.0d0,a(n))*b(n)   &
     &                                  ,sign(1.0d0,a(n))*c(n))) 
      enddo

      return
      end subroutine MClimiter

      subroutine NumericalFlux
      use commons, only: is, ie, in, d, v, p
      use fluxmod
      implicit none
      integer::i
      real(8),dimension(in):: dl,dr
      real(8),dimension(in):: vl,vr
      real(8),dimension(in):: pl,pr
      real(8),dimension(nhyd):: dsvp,dsvm,dsvc,dsv
      real(8),dimension(nhyd,in):: leftpr,rigtpr
      real(8),dimension(2*mflx+madd,in):: leftco,rigtco
      real(8),dimension(2*mflx+madd):: leftst,rigtst
      real(8),dimension(mflx):: nflux
      real(8) :: ddm, ddp, dvm, dvp, dpm, dpp
      real(8) :: ddmon, dvmon, dpmon

      do i=is-1,ie+1
         ddp = d(i+1) - d(i)
         ddm = d(i) - d(i-1)
         call vanLeer(ddp,ddm,ddmon)

         dvp = v(i+1) - v(i)
         dvm = v(i) - v(i-1)
         call vanLeer(dvp,dvm,dvmon)

         dpp = p(i+1) - p(i)
         dpm = p(i) - p(i-1)
         call vanLeer(dpp,dpm,dpmon)


         dl(i+1) = d(i) + 0.5d0*ddmon
         dr(i  ) = d(i) - 0.5d0*ddmon

         vl(i+1) = v(i) + 0.5d0*dvmon
         vr(i  ) = v(i) - 0.5d0*dvmon

         pl(i+1) = p(i) + 0.5d0*dpmon
         pr(i  ) = p(i) - 0.5d0*dpmon
      enddo

!         call HLLE(leftst,rigtst,nflux)
      call HLLE(dl,vl,pl,dr,vr,pr)

      return
      end subroutine Numericalflux

!      subroutine HLLE(leftst,rigtst,nflux)
      subroutine HLLE(dl,vl,pl,dr,vr,pr)
      use commons, only : is, ie
      use eosmod
      use fluxmod
      implicit none
      real(8),dimension(in),intent(in)::dl,vl,pl,dr,vr,pr
      real(8):: csl,csr
      real(8):: sl, sr
      real(8):: dfluxl, dfluxr
      real(8):: mvl, mvr, mvfluxl, mvfluxr
      real(8):: etl, etr, etfluxl, etfluxr
      integer :: i

      do i=is,ie+1
          ! conserved variables in the left and right states
          mvl = dl(i)*vl(i)
          mvr = dr(i)*vr(i)
          etl = pl(i)/(gam - 1.0d0) + 0.5d0*dl(i)*vl(i)**2
          etr = pr(i)/(gam - 1.0d0) + 0.5d0*dr(i)*vr(i)**2

          ! flux in the left and right states
          dfluxl = mvl
          dfluxr = mvr

          mvfluxl = pl(i) + dl(i)*vl(i)**2
          mvfluxr = pr(i) + dr(i)*vr(i)**2

          etfluxl = ( gam*pl(i)/(gam - 1.0d0) + 0.5d0*dl(i)*vl(i)**2 )*vl(i)
          etfluxr = ( gam*pr(i)/(gam - 1.0d0) + 0.5d0*dr(i)*vr(i)**2 )*vr(i)

          csl = dsqrt(gam*pl(i)/dl(i))
          csr = dsqrt(gam*pr(i)/dr(i))

          sl = min(vl(i),vr(i)) - max(csl,csr)
!          sl = min(0.0d0,sl)
          sr = max(vl(i),vr(i)) + max(csl,csr)
!          sr = max(0.0d0,sr)

          dflux(i)  = (sr*dfluxl - sl*dfluxr + sl*sr*( dr(i) - dl(i) ))/(sr - sl)
          mvflux(i) = (sr*mvfluxl - sl*mvfluxr + sl*sr*( mvr - mvl ))/(sr - sl)
          etflux(i) = (sr*etfluxl - sl*etfluxr + sl*sr*( etr - etl ))/(sr - sl)
      enddo

      return
      end subroutine HLLE

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

      subroutine UpdateConsv
      use commons
      use fluxmod
      implicit none
      integer::i

      do i=is,ie
         
         d(i) = d(i)                       &
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

      subroutine Output
      use commons
      implicit none
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

      if(time .lt. tout+dtout) return

      write(filename,'(a3,i5.5,a4)')"unf",nout,".dat"
      filename = trim(dirname)//filename

      open(unitout,file=filename,status='replace',form='formatted')
      write(unitout,*) "# ",time,dt
      write(unitout,*) "# ",ngrid,gs
      write(unitout,*) "# ",ngrid,gs
      close(unitout)

      x1out(is-gs:ie+gs,1) = x1b(is-gs:ie+gs)
      x1out(is-gs:ie+gs,2) = x1a(is-gs:ie+gs)

      hydout(is-gs:ie+gs,1) =  d(is-gs:ie+gs)
      hydout(is-gs:ie+gs,2) = v(is-gs:ie+gs)
      hydout(is-gs:ie+gs,3) =  p(is-gs:ie+gs)

      write(filename,'(a3,i5.5,a4)')"bin",nout,".dat"
      filename = trim(dirname)//filename
!      open(unitbin,file=filename,status='replace',form='formatted') 
      open(unitbin,file=filename,form='formatted',action="write")
      do i=1,in
          write(unitbin,*) x1b(i), d(i), v(i), p(i)
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
      tout=time

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
