module commons
implicit none
integer::ntime                        ! counter of the timestep
integer,parameter::ntimemax=20000     ! the maximum timesteps
real(8)::time,dt                    ! time, timewidth
data time / 0.0d0 /
real(8),parameter:: timemax=1.0d0/dsqrt(5.0d0/3.0d0)
real(8),parameter:: dtout=5.0d-3

integer,parameter::nx=128  ! the number of grids in the simulation box
integer,parameter::mgn=2   ! the number of ghost cells
integer,parameter::nxtot=nx+2*mgn+1 ! the total number of grids including ghost cells
integer,parameter::is=mgn+1         ! the index of the leftmost grid
integer,parameter::ie=nx+mgn     ! the index of the rightmost grid
integer,parameter::nflux = 3
real(8),parameter:: x1min=-0.5d0,x1max=0.5d0

integer, parameter :: IDN = 1
integer, parameter :: IMX = 2
integer, parameter :: IPR = 3
integer, parameter :: NVAR = 3

integer, parameter :: IVX = 2
integer, parameter :: IEN = 3

real(8),parameter::gam=5.0d0/3.0d0 !! adiabatic index

end module commons
     
program main
use commons
implicit none

integer :: i
real(8),dimension(nxtot)::xf,xv
real(8),dimension(nxtot,NVAR) :: U ! conservative variables
real(8),dimension(nxtot,NVAR) :: Q ! primitive variables
real(8),dimension(nxtot,NVAR) :: F ! numerical flux
real(8) :: error

      write(6,*) "setup grids and initial condition"
      call GenerateGrid(xf, xv)
      call GenerateProblem(xv, Q)
      call ConsvVariable(Q, U)
      call Output( xf, xv, Q )

!      open(1,file="sound128.dat",action="write")
! main loop
      do !ntime=1,ntimemax
         call TimestepControl(xf, Q)
         if( time + dt > timemax ) dt = timemax - time
         call BoundaryCondition(Q)
         call NumericalFlux(Q, F)
         call UpdateConsv( F, U )
         call PrimVariable( U, Q )
         time=time+dt
         call Output( xf, xv, Q)

!         write(1,*) time, 0.5d0*(Q(2,IDN) + Q(3,IDN))
         if(time >= timemax) exit 
      enddo 

      call Output( xf, xv, Q)
!      do i=is,ie
!         print*,xv(i),Q(i,IDN),1.0d0 + 1.0d-4*dsin(2.0d0*dcos(-1.0d0)*xv(i))
!      enddo
!      close(1)

      write(6,*) "program has been finished"
contains

subroutine GenerateGrid(xf, xv)
use commons
implicit none
real(8), intent(out) :: xf(:), xv(:)
real(8) :: dx,dy
integer::i

    dx=(x1max-x1min)/nx
    do i=1,nxtot
         xf(i) = dx*(i-(mgn+1))+x1min
    enddo
    do i=1,nxtot-1
         xv(i) = 0.5d0*(xf(i+1)+xf(i))
    enddo

return
end subroutine GenerateGrid

subroutine GenerateProblem(xv, Q)
use commons
implicit none
integer::i
real(8), intent(in ) :: xv(:)
real(8), intent(out) :: Q(:,:)
real(8) :: rho1,rho2,Lsm,u1,u2
real(8)::pi,amp

      pi=acos(-1.0d0)
      amp = 1.0d-4

      do i=is,ie
             Q(i,IDN) = 1.0d0 + amp*dsin(2.0d0*pi*xv(i))
             Q(i,IVX) = amp*dsqrt(gam)*dsin(2.0d0*pi*xv(i))
             Q(i,IPR) = 1.0d0 + amp*gam*dsin(2.0d0*pi*xv(i))
      enddo

return
end subroutine GenerateProblem

subroutine BoundaryCondition(Q)
use commons
implicit none
real(8), intent(inout) :: Q(:,:)
integer::i

    do i=1,mgn 
         Q(is-i,IDN)  = Q(ie+1-i,IDN)
         Q(is-i,IVX)  = Q(ie+1-i,IVX)
         Q(is-i,IPR)  = Q(ie+1-i,IPR)
    enddo

    do i=1,mgn
         Q(ie+i,IDN) = Q(is+i-1,IDN)
         Q(ie+i,IVX) = Q(is+i-1,IVX)
         Q(ie+i,IPR) = Q(is+i-1,IPR)
    enddo

return
end subroutine BoundaryCondition
!
subroutine ConsvVariable(Q, U)
use commons
implicit none
real(8), intent(in) :: Q(:,:)
real(8), intent(out) :: U(:,:)
integer::i

    do i=is,ie
        U(i,IDN) = Q(i,IDN)
        U(i,IMX) = Q(i,IDN)*Q(i,IVX)
        U(i,IEN)  = 0.5d0*Q(i,IDN)*Q(i,IVX)**2 + Q(i,IPR)/(gam - 1.0d0)
    enddo
      
return
end subroutine Consvvariable

subroutine PrimVariable( U, Q )
use commons
implicit none
real(8), intent(in) :: U(:,:)
real(8), intent(out) :: Q(:,:)
integer::i

    do i=is,ie
        Q(i,IDN) = U(i,IDN)
        Q(i,IVX) = U(i,IMX)/U(i,IDN)
        Q(i,IPR) = ( U(i,IEN) - 0.5d0*U(i,IMX)**2/U(i,IDN) )*(gam-1.0d0)
    enddo

return
end subroutine PrimVariable

subroutine TimestepControl(xf, Q)
use commons
implicit none
real(8), intent(in) :: xf(:), Q(:,:)
real(8)::dtl1
real(8)::dtlocal
real(8)::dtmin
integer::i

    dtmin=1.0d90

    do i=is,ie
        dtlocal =(xf(i+1)-xf(i))/(abs(Q(i,IVX)) + dsqrt(gam*Q(i,IPR)/Q(i,IDN)))
         if(dtlocal .lt. dtmin) dtmin = dtlocal
    enddo

    dt = 0.0125d0 * dtmin

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
subroutine NumericalFlux( Q, F )
use commons !, only: is, ie, in
implicit none
real(8), intent(in) :: Q(:,:)
real(8), intent(out) :: F(:,:)
integer::i
real(8),dimension(nxtot,NVAR):: Ql,Qr
real(8),dimension(NVAR):: flx
real(8) :: dQm(NVAR), dQp(NVAR), dQmon(NVAR)
real(8) :: ddmon, dvmon, dpmon

    do i=is-1,ie+1
        dQp(:) = Q(i+1,:) - Q(i  ,:)
        dQm(:) = Q(i  ,:) - Q(i-1,:)

        call vanLeer(NVAR, dQp,dQm,dQmon)

        Ql(i,:) = Q(i,:) + 0.5d0*dQmon(:)*1.0d0!*1.0d0*(1.0d0 - dsqrt(gam*Q(i,IPR)/Q(i,IDN))*dt)
        Qr(i-1,:) = Q(i,:) - 0.5d0*dQmon(:)*1.0d0!*1.0d0*(1.0d0 - dsqrt(gam*Q(i,IPR)/Q(i,IDN))*dt)
    enddo

!         call HLLE(leftst,rigtst,nflux)
    do i=is-1,ie
        call HLL(Ql(i,:),Qr(i,:),flx)
!        call HLLC(Ql(i,:),Qr(i,:),flx)
!        call Lax((xf(i) - xf(i-1))/dt,Ql(i,:),Qr(i,:),flx)
        F(i,:)  = flx(:)
    enddo

return
end subroutine Numericalflux

subroutine HLL(Ql,Qr,flx)
use commons !, only : is, ie, NVAR
implicit none
real(8),intent(in)::Ql(:), Qr(:)
real(8),intent(out) :: flx(:)
integer :: i, n
real(8):: Ul(NVAR), Ur(NVAR)
real(8):: Fl(NVAR), Fr(NVAR)
real(8):: csl,csr
real(8):: sl, sr

    ! conserved variables in the left and right states
    Ul(IDN) = Ql(IDN)
    Ur(IDN) = Qr(IDN)

    Ul(IMX) = Ql(IDN)*Ql(IVX)
    Ur(IMX) = Qr(IDN)*Qr(IVX)

    Ul(IEN) = 0.5d0*Ql(IDN)*Ql(IVX)**2 + Ql(IPR)/(gam - 1.0d0)
    Ur(IEN) = 0.5d0*Qr(IDN)*Qr(IVX)**2 + Qr(IPR)/(gam - 1.0d0)

    ! flux in the left and right states
    Fl(IDN) = Ul(IMX)
    Fr(IDN) = Ur(IMX)

    Fl(IMX) = Ql(IPR) + Ql(IDN)*Ql(IVX)**2 
    Fr(IMX) = Qr(IPR) + Qr(IDN)*Qr(IVX)**2 


    Fl(IEN) = ( gam*Ql(IPR)/(gam - 1.0d0) + 0.5d0*Ql(IDN)*Ql(IVX)**2)*Ql(IVX)
    Fr(IEN) = ( gam*Qr(IPR)/(gam - 1.0d0) + 0.5d0*Qr(IDN)*Qr(IVX)**2)*Qr(IVX)

    csl = dsqrt(gam*Ql(IPR)/Ql(IDN))
    csr = dsqrt(gam*Qr(IPR)/Qr(IDN))

    sl = min(Ql(IVX),Qr(IVX)) - max(csl,csr)
    sr = max(Ql(IVX),Qr(IVX)) + max(csl,csr)

    if( sl > 0.0d0 ) then 
        do n=1,NVAR
            flx(n) = Fl(n)
        enddo
    else if (sr < 0.0d0 ) then
        do n=1,NVAR
             flx(n) = Fr(n)
        enddo
    else 
        do n=1,NVAR
            flx(n)  = (sr*Fl(n) - sl*Fr(n) + sl*sr*( Ur(n) - Ul(n) ))/(sr - sl)
        enddo
    endif

return
end subroutine HLL
!
subroutine HLLC(Ql,Qr,flx)
use commons !, only : is, ie, NVAR
implicit none
real(8),intent(in)::Ql(:), Qr(:)
real(8),intent(out) :: flx(:)
integer :: i, n
real(8):: Ul(NVAR), Ur(NVAR), Ulst(NVAR), Urst(NVAR)
real(8):: Fl(NVAR), Fr(NVAR)
real(8):: csl,csr
real(8):: sl, sr, cl, cr, sm, pst

    csl = dsqrt(gam*Ql(IPR)/Ql(IDN))
    csr = dsqrt(gam*Qr(IPR)/Qr(IDN))
    sl = min(Ql(IVX),Qr(IVX)) - max(csl,csr)
    sr = max(Ql(IVX),Qr(IVX)) + max(csl,csr)

    cl = Ql(IDN)*(sl - Ql(IVX))
    cr = Qr(IDN)*(sr - Qr(IVX))

    sm = ( -cl*Ql(IVX) + cr*Qr(IVX) - Qr(IPR) + Ql(IPR) )/(-cl + cr)

    pst = Ql(IPR) + cl*(sm - Ql(IVX))
    if( abs( pst /( Qr(IPR) + cr*(sm - Qr(IVX)) ) - 1.0d0 ) > 1.0d-10 ) then
    print*,"error", pst
    stop
    endif

    ! conserved variables in the left and right states
    Ul(IDN) = Ql(IDN)
    Ur(IDN) = Qr(IDN)
    Ulst(IDN) = cl/(sl - sm)
    Urst(IDN) = cr/(sr - sm)

!    print*,Ul(IDN),Ulst(IDN)

    Ul(IMX) = Ql(IDN)*Ql(IVX)
    Ur(IMX) = Qr(IDN)*Qr(IVX)
    Ulst(IMX) = Ulst(IDN)*sm
    Urst(IMX) = Urst(IDN)*sm

    Ul(IEN) = 0.5d0*Ql(IDN)*Ql(IVX)**2 + Ql(IPR)/(gam - 1.0d0)
    Ur(IEN) = 0.5d0*Qr(IDN)*Qr(IVX)**2 + Qr(IPR)/(gam - 1.0d0)
    Ulst(IEN) = ( (sl - Ql(IVX))*Ul(IEN) - Ql(IPR)*Ql(IVX) + pst*sm )/(sl - sm)
    Urst(IEN) = ( (sr - Qr(IVX))*Ur(IEN) - Qr(IPR)*Qr(IVX) + pst*sm )/(sr - sm)

    ! flux in the left and right states
    Fl(IDN) = Ul(IMX)
    Fr(IDN) = Ur(IMX)

    Fl(IMX) = Ql(IPR) + Ql(IDN)*Ql(IVX)**2 
    Fr(IMX) = Qr(IPR) + Qr(IDN)*Qr(IVX)**2 

    Fl(IEN) = ( gam*Ql(IPR)/(gam - 1.0d0) + 0.5d0*Ql(IDN)*Ql(IVX)**2)*Ql(IVX)
    Fr(IEN) = ( gam*Qr(IPR)/(gam - 1.0d0) + 0.5d0*Qr(IDN)*Qr(IVX)**2)*Qr(IVX)


    if( sl > 0.0d0 ) then 
        do n=1,NVAR
            flx(n) = Fl(n)
        enddo
    else if (sr < 0.0d0 ) then
        do n=1,NVAR
             flx(n) = Fr(n)
        enddo
    else if (sm > 0.0d0 ) then
        do n=1,NVAR
            flx(n)  = Fl(n) + sl*(Ulst(n) - Ul(n))
        enddo
    else 
        do n=1,NVAR
            flx(n)  = Fr(n) + sr*(Urst(n) - Ur(n))
        enddo
    endif

return
end subroutine HLLC
!
subroutine Lax(dxdt,Ql,Qr,flx)
use commons !, only : is, ie, NVAR
implicit none
real(8),intent(in)::Ql(:), Qr(:)
real(8),intent(in)::dxdt
real(8),intent(out) :: flx(:)
integer :: i, n
real(8):: Ul(NVAR), Ur(NVAR)
real(8):: Fl(NVAR), Fr(NVAR)
real(8):: csl,csr
real(8):: sl, sr

    ! conserved variables in the left and right states
    Ul(IDN) = Ql(IDN)
    Ur(IDN) = Qr(IDN)

    Ul(IMX) = Ql(IDN)*Ql(IVX)
    Ur(IMX) = Qr(IDN)*Qr(IVX)

    Ul(IEN) = 0.5d0*Ql(IDN)*Ql(IVX)**2 + Ql(IPR)/(gam - 1.0d0)
    Ur(IEN) = 0.5d0*Qr(IDN)*Qr(IVX)**2 + Qr(IPR)/(gam - 1.0d0)

    ! flux in the left and right states
    Fl(IDN) = Ul(IMX)
    Fr(IDN) = Ur(IMX)

    Fl(IMX) = Ql(IPR) + Ql(IDN)*Ql(IVX)**2 
    Fr(IMX) = Qr(IPR) + Qr(IDN)*Qr(IVX)**2 

    Fl(IEN) = ( gam*Ql(IPR)/(gam - 1.0d0) + 0.5d0*Ql(IDN)*Ql(IVX)**2)*Ql(IVX)
    Fr(IEN) = ( gam*Qr(IPR)/(gam - 1.0d0) + 0.5d0*Qr(IDN)*Qr(IVX)**2)*Qr(IVX)

    do n=1,NVAR 
        flx(n)  = 0.5d0*(Fl(n) + Fr(n)) - 0.5d0*dxdt*(Ur(n) - Ul(n))
    enddo


return
end subroutine Lax
!
!
subroutine UpdateConsv( F, U )
use commons
implicit none
real(8), intent(in)  :: F(:,:)
real(8), intent(out) :: U(:,:)
integer::i,n

    do n=1,NVAR
        do i=is,ie
!            U(i,n) = U(i,n) + dt*(- F(i+1,n) + F(i,n))/(xf(i+1)-xf(i)) 
            U(i,n) = U(i,n) + dt*(- F(i,n) + F(i-1,n))/(xf(i)-xf(i-1)) 
        enddo
    enddo

return
end subroutine UpdateConsv

subroutine Output( xf, xv, Q )
use commons
implicit none
real(8), intent(in) :: xf(:), xv(:), Q(:,:)
integer::i
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
        call makedirs("snap")
        is_inited =.true.
    endif

    print*, time, tout+dtout, time+1.0d-14 .lt. tout+dtout
    if(time + 1.0d-14.lt. tout+dtout) return

!    write(filename,'(a3,i5.5,a4)') "lax",nout,".dat"
!    write(filename,'(a3,i5.5,a4)') "hll",nout,".dat"
!    write(filename,'(a4,i5.5,a4)') "hllc",nout,".dat"

    write(filename,'(a6,i5.5,a4)') "hll2nd",nout,".dat"
!    write(filename,'(a6,i5.5,a4)') "lax2nd",nout,".dat"
!    write(filename,'(a7,i5.5,a4)') "hllc2nd",nout,".dat"
    filename = trim(dirname)//filename
!      open(unitbin,file=filename,status='replace',form='formatted') 
    open(unitbin,file=filename,form='formatted',action="write")
    write(unitbin,"(a2,f6.4)") "# ",time
    write(unitbin,*) "# x, density, velocity, pressure"
    do i=is,ie
         write(unitbin,*) xv(i), Q(i,IDN), Q(i,IVX), Q(i,IPR)
    enddo
    close(unitbin)

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
