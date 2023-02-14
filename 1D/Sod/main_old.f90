program main
implicit none

! time evolution
integer :: ntime = 0    ! counter of the timestep
real(8) :: time = 0.0d0  ! time 
real(8) :: dt   = 0.0d0  ! time width
real(8),parameter:: timemax=0.2d0 ! simulation end time

! coordinate 
integer,parameter :: nx=64            ! the number of grids in the simulation box
integer,parameter :: ngh=2            ! the number of ghost cells in each side
integer,parameter :: nxtot=nx+2*ngh+1 ! the total number of grids including ghost cells
integer,parameter :: is=ngh+1         ! the index of the leftmost grid
integer,parameter :: ie=nx+ngh        ! the index of the rightmost grid
integer,parameter :: nflux = 3        ! the number of variables used in the flux calculation
real(8),parameter :: x1min=-0.5d0,x1max=0.5d0  ! ths size of the simulation box 

! indices of the primitive variables
integer, parameter :: IDN = 1
integer, parameter :: IMX = 2
integer, parameter :: IPR = 3

! indices of the conservative variables
integer, parameter :: IVX = 2
integer, parameter :: IEN = 3

! the number of the indices of the primitive and conservative variables
integer, parameter :: NVAR = 3

real(8),parameter::gam=1.4d0 !! adiabatic index

! definition of arrays 
real(8),dimension(nxtot) :: xf     ! xf(i) <=> x_{i-1/2}
real(8),dimension(nxtot) :: xv     ! xv(i) <=> x_i
real(8),dimension(nxtot,NVAR) :: U ! U(i) <=> U_i conservative variables
real(8),dimension(nxtot,NVAR) :: Q ! Q(i) <=> Q_i primitive variables
real(8),dimension(nxtot,NVAR) :: F ! F(i) <=> F_{i-1/2} numerical flux

! output 
character(20),parameter::dirname="sod_lax" ! directory name

! snapshot
integer, parameter :: unitsnap = 17
real(8),parameter:: dtout=5.0d-3  ! timewidth for snapshots

! realtime
integer, parameter :: unitevo = 11
integer, parameter :: nevo = 3    ! the number of variables derived in the realtime analysis
real(8) ::      phys_evo(nevo)    ! variables derived in the realtime analysis

    call makedirs(dirname) ! create the directory for outout

    write(6,*) "setup grids and initial condition"
    call GenerateGrid(xf, xv)    ! set coordinate 
    call GenerateProblem(xv, Q)  ! set initial condition of Q
    call Prim2Consv(Q, U)     ! Q --> U
    call Output( .true., dirname, xf, xv, Q ) !output the initial function

    write(6,*) "Start the simulation"
    open(unitevo,file=trim(dirname)//'/'//'realtime_ana.dat', action="write")
    do 
        dt = TimestepControl(xf, Q) ! determine the time width
        if( time + dt > timemax ) dt = timemax - time  ! to terminate the simulation at t = timemax
        call BoundaryCondition( Q )   ! set boundary condition 
        call NumericalFlux( Q, F )    ! calculate the numerical flux
        call UpdateConsv( dt, F, U )    ! update conservative variables
        call Consv2Prim( U, Q )   ! U --> Q

        time=time+dt
        ntime = ntime + 1

        write(6,"(a8,i5.5,a10,e15.8,a7,e15.8)") "ntime = ",ntime, ",  time = ",time, ",  dt = ",dt

        if( mod(ntime,10) .eq. 0 ) then
            call RealtimeAnalysis(xf,xv,Q,phys_evo)
            write(unitevo,*) time, phys_evo(1:nevo)
        endif

        call Output( .false., dirname, xf, xv, Q)

         if(time >= timemax) exit 
    enddo 
    close(unitevo)

    call AnalysisAfterSimu(time,xf,xv,Q)      ! analysis after simulation
    call Output( .true., dirname, xf, xv, Q)

    write(6,*) "the simulation has been finished"
contains

!-------------------------------------------------------------------
!       generate grid 
!-------------------------------------------------------------------
subroutine GenerateGrid(xf, xv)
implicit none
real(8), intent(out) :: xf(:), xv(:)
real(8) :: dx,dy
integer::i

    dx=(x1max-x1min)/nx
    do i=1,nxtot
         xf(i) = dx*(i-(ngh+1))+x1min
    enddo
    do i=1,nxtot-1
         xv(i) = 0.5d0*(xf(i+1)+xf(i))
    enddo

return
end subroutine GenerateGrid
!-------------------------------------------------------------------
!       generate initial condition 
!-------------------------------------------------------------------
subroutine GenerateProblem(xv, Q)
implicit none
integer::i
real(8), intent(in ) :: xv(:)
real(8), intent(out) :: Q(:,:)

      do i=is,ie
         if( xv(i) < 0.0d0 ) then 
             Q(i,IDN) = 1.0d0
             Q(i,IVX) = 0.0d0
             Q(i,IPR) = 1.0d0
        else 
             Q(i,IDN) = 0.125d0
             Q(i,IVX) = 0.0d0
             Q(i,IPR) = 0.1d0
         endif
      enddo

return
end subroutine GenerateProblem
!-------------------------------------------------------------------
!       Boundary Condition 
!-------------------------------------------------------------------
subroutine BoundaryCondition(Q)
implicit none
real(8), intent(inout) :: Q(:,:)
integer::i

    do i=1,ngh 
         Q(is-i,IDN)  = Q(is-1+i,IDN)
         Q(is-i,IVX)  = Q(is-1+i,IVX)
         Q(is-i,IPR)  = Q(is-1+i,IPR)
    enddo

    do i=1,ngh
         Q(ie+i,IDN) = Q(ie-i+1,IDN)
         Q(ie+i,IVX) = Q(ie-i+1,IVX)
         Q(ie+i,IPR) = Q(ie-i+1,IPR)
    enddo

return
end subroutine BoundaryCondition
!-------------------------------------------------------------------
!       Primitive variables ===> Conservative variables
!       Input  : Q
!       Output : U
!-------------------------------------------------------------------
subroutine Consv2Prim(Q, U)
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
end subroutine Consv2Prim

!-------------------------------------------------------------------
!       Conservative variables ===> Primitive variables
!       Input  : U
!       Output : Q
!-------------------------------------------------------------------
subroutine Prim2Consv( U, Q )
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
end subroutine Prim2Consv

!-------------------------------------------------------------------
!       Determine the time width
!       dt = min_i( dx_{i}/(|v_i| + c_i) )
!-------------------------------------------------------------------
real(8) function TimestepControl(xf, Q) 
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

    TimestepControl = 0.3d0 * dtmin

return
end function TimestepControl
!-------------------------------------------------------------------
!       Conservative variables ===> Primitive variables
!       Input  : U
!       Output : Q
!-------------------------------------------------------------------
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
!-------------------------------------------------------------------
!       Calculate the Numerical Flux at the cell surface from the primitive variables
!       Input  : Q
!       Output : F
!-------------------------------------------------------------------
subroutine NumericalFlux( Q, F )
implicit none
real(8), intent(in) :: Q(:,:)
real(8), intent(out) :: F(:,:)
integer::i
real(8),dimension(nxtot,NVAR):: Ql,Qr
real(8),dimension(NVAR):: flx
real(8) :: dQm(NVAR), dQp(NVAR), dQmon(NVAR)
real(8) :: ddmon, dvmon, dpmon

    do i=is-1,ie+1
        Ql(i+1,:) = Q(i,:)
        Qr(i  ,:) = Q(i,:)
    enddo

    do i=is,ie+1
          call Lax((xf(i) - xf(i-1))/dt,Ql(i,:),Qr(i,:),flx)
!        call HLL(Ql(i,:),Qr(i,:),flx)
        F(i,:)  = flx(:)
    enddo

return
end subroutine Numericalflux
!-------------------------------------------------------------------
!       The Lax Flux derived from the left and right quantities
!       Input  : dxdt, Ql, Qr
!       Output : flx
!-------------------------------------------------------------------
subroutine Lax(dxdt,Ql,Qr,flx)
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
    Ul(IMX) = Ql(IDN)*Ql(IVX)
    Ul(IEN) = 0.5d0*Ql(IDN)*Ql(IVX)**2 + Ql(IPR)/(gam - 1.0d0)

    Ur(IDN) = Qr(IDN)
    Ur(IMX) = Qr(IDN)*Qr(IVX)
    Ur(IEN) = 0.5d0*Qr(IDN)*Qr(IVX)**2 + Qr(IPR)/(gam - 1.0d0)

    ! flux in the left and right states
    Fl(IDN) = Ul(IMX)
    Fl(IMX) = Ql(IPR) + Ql(IDN)*Ql(IVX)**2 
    Fl(IEN) = ( gam*Ql(IPR)/(gam - 1.0d0) + 0.5d0*Ql(IDN)*Ql(IVX)**2)*Ql(IVX)

    Fr(IDN) = Ur(IMX)
    Fr(IMX) = Qr(IPR) + Qr(IDN)*Qr(IVX)**2 
    Fr(IEN) = ( gam*Qr(IPR)/(gam - 1.0d0) + 0.5d0*Qr(IDN)*Qr(IVX)**2)*Qr(IVX)

    do n=1,NVAR 
        flx(n)  = 0.5d0*(Fl(n) + Fr(n)) - 0.5d0*dxdt*(Ur(n) - Ul(n))
    enddo

return
end subroutine Lax
!-------------------------------------------------------------------
!       The HLL Flux derived from the left and right quantities
!       Input  : dxdt, Ql, Qr
!       Output : flx
!-------------------------------------------------------------------
subroutine HLL(Ql,Qr,flx)
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
!-------------------------------------------------------------------
!       The HLLC Flux derived from the left and right quantities
!       Input  : dxdt, Ql, Qr
!       Output : flx
!-------------------------------------------------------------------
subroutine HLLC(Ql,Qr,flx)
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
!-------------------------------------------------------------------
!       Update the conservative variables
!       Input  : dt1, F
!       Output : U
!-------------------------------------------------------------------
subroutine UpdateConsv( dt, F, U )
implicit none
real(8), intent(in)  :: dt
real(8), intent(in)  :: F(:,:)
real(8), intent(out) :: U(:,:)
integer::i,n

    do n=1,NVAR
        do i=is,ie
            U(i,n) = U(i,n) + dt*(- F(i+1,n) + F(i,n))/(xf(i+1)-xf(i)) 
        enddo
    enddo

return
end subroutine UpdateConsv
!-------------------------------------------------------------------
!       Output snapshot files 
!       Input  : flag, xf, xv, Q
!
!       flag = .true.  --> output snapshot when calling this subroutine
!       flag = .false. --> output snapshot every dtout
!-------------------------------------------------------------------
subroutine Output( flag, dirname, xf, xv, Q )
implicit none
logical,       intent(in) :: flag
character(20), intent(in) :: dirname 
real(8),       intent(in) :: xf(:), xv(:), Q(:,:)
integer :: i
character(40) :: filename
real(8), save :: tout = - dtout
integer :: nout = 0
integer,parameter:: gs=1

logical, save:: is_inited
data is_inited /.false./

!    if (.not. is_inited) then
!        call makedirs("snap")
!        is_inited =.true.
!    endif

    if( .not.flag) then
          if( time + 1.0d-14.lt. tout+dtout) return
    endif

    write(filename,'(i5.5)') nout
    filename = trim(dirname)//"/snap"//trim(filename)//".dat"
    open(unitsnap,file=filename,form='formatted',action="write")
    write(unitsnap,"(a2,f6.4)") "# ",time
    do i=is,ie
         write(unitsnap,*) xv(i), Q(i,IDN), Q(i,IVX), Q(i,IPR)
    enddo
    close(unitsnap)

    nout=nout+1
    tout=tout + dtout

return
end subroutine Output
!-------------------------------------------------------------------
!       Realtime Analysis
!       Input  : xf, xv
!       Output : phys_evo(nevo)
!-------------------------------------------------------------------
subroutine RealtimeAnalysis(xf,xv,Q,phys_evo)
real(8), intent(in)  :: xf(:), xv(:), Q(:,:)
real(8), intent(out) :: phys_evo(:)
integer :: i,j,k
real(8) :: tmp

      tmp = 0.0d0
      do i=is,ie
           tmp = tmp + 1.0d0
      enddo
      phys_evo(1) = tmp/dble(nx)
      
return
end subroutine
!-------------------------------------------------------------------
!       Analysis after simulation
!       Input  : xf, xv
!       Output : phys_evo(nevo)
!-------------------------------------------------------------------
subroutine AnalysisAfterSimu(time,xf,xv,Q)
real(8), intent(in)  :: xf(:), xv(:), Q(:,:)
real(8), intent(in)  :: time
integer :: i
real(8) :: tmp

      tmp = 0.0d0
      do i=is,ie
           tmp = tmp + 1.0d0
      enddo
      
return
end subroutine
!-------------------------------------------------------------------
!       create directory
!       Input  : the directory to be created
!-------------------------------------------------------------------
subroutine makedirs(outdir)
implicit none
character(len=*), intent(in) :: outdir
character(len=256) command
    write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
    write(*, *) trim(command)
    call system(command)
end subroutine makedirs

end program main
