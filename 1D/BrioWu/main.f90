program main
implicit none

! time evolution
integer :: ntime = 0    ! counter of the timestep
real(8) :: time = 0.0d0  ! time 
real(8) :: dt   = 0.0d0  ! time width
real(8),parameter:: timemax=0.1d0 ! simulation end time

integer,parameter::nx=256*1       ! the number of grids in the simulation box
integer,parameter::ngh=2            ! the number of ghost cells
integer,parameter::nxtot=nx+2*ngh+1 ! the total number of grids including ghost cells
integer,parameter::is=ngh+1         ! the index of the leftmost grid
integer,parameter::ie=nx+ngh     ! the index of the rightmost grid
real(8),parameter:: xmin=-0.5d0,xmax=0.5d0

real(8),parameter:: Bx=0.75

! indices of the primitive variables
integer, parameter :: IDN = 1
integer, parameter :: IMX = 2
integer, parameter :: IMY = 3
integer, parameter :: IMZ = 4
integer, parameter :: IPR = 5
integer, parameter :: IBY = 6
integer, parameter :: IBZ = 7
integer, parameter :: NVAR = 7
integer, parameter :: NFLX = 7

! indices of the conservative variables
integer, parameter :: IVX = 2
integer, parameter :: IVY = 3
integer, parameter :: IVZ = 4
integer, parameter :: IEN = 5

! adiabatic index
real(8),parameter::gam=5.0d0/3.0d0 


! definition of arrays 
real(8),dimension(nxtot)::xf,xv
real(8),dimension(nxtot,NVAR) :: U
real(8),dimension(nxtot,NVAR) :: Q
real(8),dimension(nxtot,NFLX) :: flux

! output 
character(20),parameter::dirname="lax" ! directory name

! snapshot
integer, parameter :: unitsnap = 17

! realtime analysis 
integer, parameter :: unitevo = 11
integer, parameter :: nevo = 1    ! the number of variables derived in the realtime analysis
real(8) ::  phys_evo(nevo)    ! variables derived in the realtime analysis

     ! make the directory for output
     call makedirs(trim(dirname))

      write(6,*) "setup grids and initial condition"
      call GenerateGrid(xf, xv)
      call GenerateProblem(xv, Q)
      call Prim2Consv(Q, U)
      call Output( .TRUE., dirname, xf, xv, Q )

      write(6,*) "Start the simulation"
      !open file to output the result of the realtime analysis
      open(unitevo,file=trim(dirname)//'/'//'ana.dat', action="write")
! main loop
      do 
         dt = TimestepControl(xf, Q)
         if( time + dt > timemax ) dt = timemax - time
         call BoundaryCondition( Q)
         call NumericalFlux( Q, flux )
         call UpdateConsv( dt, xf, flux, U)
         call Consv2Prim( U, Q )
         time=time+dt
         print*,"time = ",time, "dt = ",dt
         ntime = ntime + 1
         call Output( .FALSE., dirname, xf, xv, Q )

         if( mod(ntime,10) .eq. 0 ) then
            call RealtimeAnalysis(xf,xv,Q,phys_evo)
            write(unitevo,*) time, phys_evo(1:nevo)
         endif

         if(time >= timemax) exit 
      enddo 
      call Output( .TRUE., dirname, xf, xv, Q )
      call AnalysisAfterSimu(time,xf,xv,Q)

!      write(6,*) "program has been finished"
contains

!-------------------------------------------------------------------
!       Generate coordiantes
!       xf --> cell boundary xf(i) <==> x_{i-1/2}
!       xv --> cell center   xv(i) <==> x_{i}
!-------------------------------------------------------------------
subroutine GenerateGrid(xf, xv)
implicit none
real(8), intent(out) :: xf(:), xv(:)
real(8) :: dx,dy
integer::i

    dx=(xmax-xmin)/nx
    do i=1,nxtot 
        xf(i) = dx*(i-(ngh+1))+xmin
    enddo

    do i=1,nxtot-1
        xv(i) = 0.5d0*(xf(i+1)+xf(i))
    enddo

return
end subroutine GenerateGrid

!-------------------------------------------------------------------
!       Generate initial condition of the primitive variables
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
             Q(i,IVY) = 0.0d0
             Q(i,IVZ) = 0.0d0
             Q(i,IPR) = 1.0d0
             Q(i,IBY) = 1.0d0
             Q(i,IBZ) = 0.0d0
        else 
             Q(i,IDN) = 0.125d0
             Q(i,IVX) = 0.0d0
             Q(i,IVY) = 0.0d0
             Q(i,IVZ) = 0.0d0
             Q(i,IPR) = 0.1d0
             Q(i,IBY) = -1.0d0
             Q(i,IBZ) = 0.0d0
         endif
    enddo

return
end subroutine GenerateProblem

!-------------------------------------------------------------------
!       Boundary Condition of the primitive variables
!-------------------------------------------------------------------
subroutine BoundaryCondition(Q)
implicit none
real(8), intent(inout) :: Q(:,:)
integer::i,ihy

      do ihy=1,NVAR
      do i=1,ngh
          Q(is-i,ihy)  = Q(is-1+i,ihy)
      enddo
      enddo

      do ihy=1,NVAR
      do i=1,ngh
          Q(ie+i,ihy) = Q(ie-i+1,ihy)
      enddo
      enddo

return
end subroutine BoundaryCondition
!-------------------------------------------------------------------
!       Primitive variables ===> Conservative variables
!       Input  : Q
!       Output : U
!-------------------------------------------------------------------
subroutine Prim2Consv(Q, U)
implicit none
real(8), intent(in) :: Q(:,:)
real(8), intent(out) :: U(:,:)
integer::i

      do i=is,ie
          U(i,IDN) = Q(i,IDN)
          U(i,IMX) = Q(i,IDN)*Q(i,IVX)
          U(i,IMY) = Q(i,IDN)*Q(i,IVY)
          U(i,IMZ) = Q(i,IDN)*Q(i,IVZ)
          U(i,IEN) = 0.5d0*Q(i,IDN)*( Q(i,IVX)**2 + Q(i,IVY)**2 + Q(i,IVZ)**2 ) &
                   + 0.5d0*( Bx**2 + Q(i,IBY)**2 + Q(i,IBZ)**2 ) &
                   + Q(i,IPR)/(gam - 1.0d0)
          U(i,IBY) = Q(i,IBY)
          U(i,IBZ) = Q(i,IBZ)
      enddo
      
return
end subroutine Prim2Consv
!-------------------------------------------------------------------
!       Conservative variables ===> Primitive variables
!       Input  : U
!       Output : Q
!-------------------------------------------------------------------
subroutine Consv2Prim( U, Q )
implicit none
real(8), intent(in) :: U(:,:)
real(8), intent(out) :: Q(:,:)
integer::i
real(8) :: inv_d;

    do i=is,ie
        Q(i,IDN) = U(i,IDN)
        inv_d = 1.0d0/U(i,IDN)
        Q(i,IVX) = U(i,IMX)*inv_d
        Q(i,IVY) = U(i,IMY)*inv_d
        Q(i,IVZ) = U(i,IMZ)*inv_d
        Q(i,IPR) = ( U(i,IEN) &
                    - 0.5d0*(U(i,IMX)**2 + U(i,IMY)**2 + U(i,IMZ)**2)*inv_d  &
                    - 0.5d0*(Bx**2 + U(i,IBY)**2 + U(i,IBZ)**2) )*(gam-1.0d0)
        Q(i,IBY) = U(i,IBY)
        Q(i,IBZ) = U(i,IBZ)
    enddo

return
end subroutine Consv2Prim
!-------------------------------------------------------------------
!       determine dt 
!-------------------------------------------------------------------
Real(8) Function TimestepControl(xf, Q) 
implicit none
real(8), intent(in) :: xf(:), Q(:,:)
real(8)::dtl1
real(8)::dtl2
real(8)::dtl3
real(8)::dtlocal
real(8)::dtmin,cf
integer::i

    dtmin=1.0d90

    do i=is,ie
         cf = dsqrt( (gam*Q(i,IPR) + Bx**2 + Q(i,IBY)**2 + Q(i,IBZ)**2)/Q(i,IDN))
         dtlocal =(xf(i+1)-xf(i))/(abs(Q(i,IVX)) + cf)
         if(dtlocal .lt. dtmin) dtmin = dtlocal
    enddo

      dt = 0.1d0 * dtmin

return
end function TimestepControl
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
subroutine NumericalFlux( Q, flux )
implicit none
real(8), intent(in) :: Q(:,:)
real(8), intent(out) :: flux(:,:)
real(8),dimension(nxtot,NFLX):: Ql,Qr
real(8),dimension(NFLX):: flx
integer::i,ihy
real(8) :: dQm, dQp, dQ

      do ihy=1,NVAR
      do i=is-1,ie+1

         dQ = 0.0d0

         Ql(i+1,ihy) = Q(i,ihy) + 0.5d0*dQ
         Qr(i  ,ihy) = Q(i,ihy) - 0.5d0*dQ
      enddo
      enddo

      do i=is,ie+1
         call Lax((xv(i) - xv(i-1))/dt,Ql(i,:),Qr(i,:),flx(:))
         flux(i,:)  = flx(:)
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
real(8),intent(in)  :: Ql(:), Qr(:)
real(8),intent(in)  :: dxdt
real(8),intent(out) :: flx(:)
integer :: i, n
real(8):: Ul(NVAR), Ur(NVAR)
real(8):: Fl(NVAR), Fr(NVAR)
real(8):: csl,csr,pbl,pbr,ptotl,ptotr
real(8):: sl, sr


    pbl = 0.5d0*(Bx**2 + Ql(IBY)**2 + Ql(IBZ)**2)
    pbr = 0.5d0*(Bx**2 + Qr(IBY)**2 + Qr(IBZ)**2)
    ptotl = Ql(IPR) + pbl
    ptotr = Qr(IPR) + pbr

    ! conserved variables in the left and right states
    Ul(IDN) = Ql(IDN)
    Ul(IMX) = Ql(IDN)*Ql(IVX)
    Ul(IMY) = Ql(IDN)*Ql(IVY)
    Ul(IMZ) = Ql(IDN)*Ql(IVZ)
    Ul(IEN) = 0.5d0*Ql(IDN)*( Ql(IVX)**2 + Ql(IVY)**2 + Ql(IVZ)**2) & 
                  + pbl + Ql(IPR)/(gam - 1.0d0)
    Ul(IBY) = Ql(IBY)
    Ul(IBZ) = Ql(IBZ)

    Ur(IDN) = Qr(IDN)
    Ur(IMX) = Qr(IDN)*Qr(IVX)
    Ur(IMY) = Qr(IDN)*Qr(IVY)
    Ur(IMZ) = Qr(IDN)*Qr(IVZ)
    Ur(IEN) = 0.5d0*Qr(IDN)*( Qr(IVX)**2 + Qr(IVY)**2 + Qr(IVZ)**2) & 
                  + pbr + Qr(IPR)/(gam - 1.0d0)
    Ur(IBY) = Qr(IBY)
    Ur(IBZ) = Qr(IBZ)

    ! flux in the left and right states
    Fl(IDN) = Ul(IMX)
    Fl(IMX) = Ul(IMX)*Ql(IVX) + ptotl - Bx**2
    Fl(IMY) = Ul(IMY)*Ql(IVX) - Bx*Ql(IBY)
    Fl(IMZ) = Ul(IMZ)*Ql(IVX) - Bx*Ql(IBZ)
    Fl(IEN) = ( Ul(IEN) + ptotl )*Ql(IVX) &
            - Bx*( Bx*Ql(IVX) + Ql(IBY)*Ql(IVY) + Ql(IBZ)*Ql(IVZ) )
    Fl(IBY) = Ql(IBY)*Ql(IVX) - Bx*Ql(IVY)
    Fl(IBZ) = Ql(IBZ)*Ql(IVX) - Bx*Ql(IVZ)

    Fr(IDN) = Ur(IMX)
    Fr(IMX) = Ur(IMX)*Qr(IVX) + ptotr - Bx**2
    Fr(IMY) = Ur(IMY)*Qr(IVX) - Bx*Qr(IBY)
    Fr(IMZ) = Ur(IMZ)*Qr(IVX) - Bx*Qr(IBZ)
    Fr(IEN) = ( Ur(IEN) + ptotr )*Qr(IVX) &
                  - Bx*( Bx*Qr(IVX) + Qr(IBY)*Qr(IVY) + Qr(IBZ)*Qr(IVZ) )
    Fr(IBY) = Qr(IBY)*Qr(IVX) - Bx*Qr(IVY)
    Fr(IBZ) = Qr(IBZ)*Qr(IVX) - Bx*Qr(IVZ)

    do n=1,NVAR 
        flx(n)  = 0.5d0*(Fl(n) + Fr(n)) - 0.5d0*dxdt*(Ur(n) - Ul(n))
    enddo


return
end subroutine Lax

!---------------------------------------------------------------------
!     HLL Riemann Solver
!---------------------------------------------------------------------
subroutine HLL(Ql,Qr,flx)
implicit none
real(8),intent(in)::Ql(:), Qr(:)
real(8),intent(out) :: flx(:)
real(8):: Ul(NFLX), Ur(NFLX)
real(8):: Fl(NFLX), Fr(NFLX)
real(8):: Fst(NFLX)
real(8):: cfl,cfr
real(8):: sl, sr
real(8):: pbl, pbr, ptotl, ptotr
integer :: i, n




return
end subroutine HLL

!-------------------------------------------------------------------
!       Update consevative variables U using numerical flux F
!-------------------------------------------------------------------
subroutine UpdateConsv( dt, xf, flux, U )
implicit none
real(8), intent(in)  :: flux(:,:), dt, xf(:)
real(8), intent(out) :: U(:,:)
integer::i,n

      do n=1,NVAR
      do i=is,ie
         U(i,n) = U(i,n) + dt*(- flux(i+1,n) + flux(i,n))/(xf(i+1)-xf(i)) 
      enddo
      enddo


return
end subroutine UpdateConsv
!-------------------------------------------------------------------
!       Output snapshot files 
!       Input  : flag, dirname, xf, xv, Q
!
!       flag = .true.  --> output snapshot when calling this subroutine
!       flag = .false. --> output snapshot every dtsnap
!-------------------------------------------------------------------
subroutine Output( flag, dirname, xf, xv, Q )
implicit none
logical,       intent(in) :: flag
character(20), intent(in) :: dirname 
real(8),       intent(in) :: xf(:), xv(:), Q(:,:)
real(8), parameter:: dtsnap=5.0d-3
integer::i
character(40) :: filename
real(8), save :: tsnap = - dtsnap
integer, save :: nsnap = 0


    if( .not.flag) then
          if( time + 1.0d-14.lt. tsnap+dtsnap) return
    endif

    write(filename,'(i5.5)') nsnap
    filename = trim(dirname)//"/snap"//trim(filename)//".dat"
    open(unitsnap,file=filename,form='formatted',action="write")
    write(unitsnap,"(a2,f6.4)") "# ",time
    do i=is,ie
          write(unitsnap,*) xv(i), Q(i,IDN), Q(i,IVX), Q(i,IVY), Q(i,IVZ), Q(i,IPR), Bx, &
          Q(i,IBY), Q(i,IBZ)
    enddo
    close(unitsnap)

    print*, "output time=", time, trim(filename)

    nsnap=nsnap+1
    tsnap=tsnap + dtsnap

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

      phys_evo(1) = 0.0d0
      
return
end subroutine
!-------------------------------------------------------------------
!       Analysis after simulation
!       Input  : xf, xv, Q
!-------------------------------------------------------------------
subroutine AnalysisAfterSimu(time,xf,xv,Q)
real(8), intent(in)  :: xf(:), xv(:), Q(:,:)
real(8), intent(in)  :: time
integer :: i
real(8) :: error

      error = 0.0d0
      do i=is,ie
           error = error + 1.0d0
      enddo
      print*, nx, error
      
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
!      write(*, *) trim(command)
      call system(command)
end subroutine makedirs

end program main
