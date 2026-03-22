program main2d
!$ use omp_lib
use params,      only: nxtot, nytot, NVAR, NFLX, dirname, unitevo, timemax, nevo, gam
use prim_mod,    only: Prim2Consv, Consv2Prim
use flux_mod,    only: NumericalFlux, TimestepControl, UpdateConsv, BoundaryCondition, SrcTerms
implicit none

! time evolution
integer :: ntime = 0    ! counter of the timestep
real(8) :: time  = 0.0d0  ! time
real(8) :: dt    = 0.0d0  ! time width

! definition of arrays
real(8), dimension(nxtot)             :: xf, xv
real(8), dimension(nytot)             :: yf, yv
real(8), dimension(NVAR,nxtot,nytot) :: Uo
real(8), dimension(NVAR,nxtot,nytot) :: U
real(8), dimension(NVAR,nxtot,nytot) :: Q
real(8), dimension(3,nxtot,nytot)    :: Bso
real(8), dimension(3,nxtot,nytot)    :: Bs
real(8), dimension(3,nxtot,nytot)    :: Bc
real(8), dimension(NVAR,nxtot,nytot) :: F
real(8), dimension(NVAR,nxtot,nytot) :: G
real(8), dimension(3,nxtot,nytot)    :: E

real(8) :: phys_evo(nevo)
integer :: i, j

  ! make the directory for output
  call makedirs(trim(dirname))

  write(6,*) "setup grids and initial condition"
  call GenerateGrid(xf, xv, yf, yv)
  call GenerateProblem(xv, yv, Q, Bs, Bc)
  call Prim2Consv(Q, Bc, U)
!$omp parallel default(shared)
  call BoundaryCondition(Q, Bs, Bc)
!$omp end parallel
  call Output(time, .TRUE., dirname, xv, yv, Q, Bc)

  write(6,*) "Start the simulation"
  open(unitevo, file=trim(dirname)//'/'//'ana.dat', action="write")

  ! main loop
  do
    dt = TimestepControl(xf, yf, Q, Bc)
    if (time + dt > timemax) dt = timemax - time

!$omp parallel default(shared)
!$omp do collapse(2) schedule(static) private(i,j)
    do j = 1, nytot
    do i = 1, nxtot
      Uo(:,i,j)  = U(:,i,j)
      Bso(:,i,j) = Bs(:,i,j)
    enddo
    enddo
!$omp end do

    call NumericalFlux(0.5d0*dt, xf, yf, Q, Bc, Bs, F, G, E)
    call UpdateConsv(0.5d0*dt, xf, xv, yf, yv, F, G, E, Uo, Bso, U, Bs)
    call Consv2Prim(U, Bs, Q, Bc)
    call BoundaryCondition(Q, Bs, Bc)

    call NumericalFlux(dt, xf, yf, Q, Bc, Bs, F, G, E)
    call UpdateConsv(dt, xf, xv, yf, yv, F, G, E, Uo, Bso, U, Bs)
    call Consv2Prim(U, Bs, Q, Bc)
    call BoundaryCondition(Q, Bs, Bc)

!$omp end parallel

    time  = time + dt
    ntime = ntime + 1

    if (mod(ntime, 10) == 0) then
      write(*,'(A,I0,A,ES12.5,A,ES12.5)') "ntime = ", ntime, " time = ", time, " dt = ", dt
      call RealtimeAnalysis(xv, yv, Q, Bc, Bs, phys_evo)
      write(unitevo,'(*(1X,ES24.16E3))') time, phys_evo(1:nevo)
    endif
    call Output(time, .FALSE., dirname, xv, yv, Q, Bc)

    if (time >= timemax) exit
  enddo
  close(unitevo)

  write(6,*) "program has been finished"

  ! suppress unused variable warnings
  if (.false.) then
    i = 0; j = 0
  endif

!=======================================================================
contains
!=======================================================================

!-------------------------------------------------------------------
!  Generate coordinates
!  xf,yf --> cell boundary xf(i) <==> x_{i-1/2}
!  xv,yv --> cell center   xv(i) <==> x_{i}
!-------------------------------------------------------------------
subroutine GenerateGrid(xf, xv, yf, yv)
use params, only : nxtot, nytot, nx, ny, ngh, xmax, xmin, ymax, ymin
implicit none
real(8), intent(out) :: xf(nxtot), xv(nxtot)
real(8), intent(out) :: yf(nytot), yv(nytot)
real(8) :: dx, dy
integer :: i, j

  dx = (xmax - xmin) / dble(nx)
  do i = 1, nxtot
    xf(i) = dx*(i-(ngh+1)) + xmin
  enddo
  do i = 1, nxtot-1
    xv(i) = 0.5d0*(xf(i+1) + xf(i))
  enddo

  dy = (ymax - ymin) / dble(ny)
  do j = 1, nytot
    yf(j) = dy*(j-(ngh+1)) + ymin
  enddo
  do j = 1, nytot-1
    yv(j) = 0.5d0*(yf(j+1) + yf(j))
  enddo

return
end subroutine GenerateGrid

!-------------------------------------------------------------------
!  Output snapshot files
!-------------------------------------------------------------------
subroutine Output(time, flag, dirname, xv, yv, Q, Bc)
use params, only : nx, ny, is, ie, js, je, NFLX, NVAR, unitsnap, unitbin, &
                  nxtot, nytot, IDN, IVX, IVY, IVZ, IPR, ISC, dtsnap, flag_binary
implicit none
real(8), intent(in)       :: time
logical, intent(in)       :: flag
character(20), intent(in) :: dirname
real(8), intent(in)       :: xv(nxtot), yv(nytot)
real(8), intent(in)       :: Q(NVAR,nxtot,nytot), Bc(3,nxtot,nytot)
integer :: i, j
character(100) :: filename
real(8), save :: tsnap = - dtsnap
integer, save :: nsnap = 0
real(8) :: Qout(NVAR,nxtot,nytot)
real(8) :: u2_loc, W_loc

  if (.not. flag) then
    if (time + 1.0d-14 < tsnap + dtsnap) return
  endif

  ! Convert 4-velocity u^i -> v^i for output (Q stores u^i = W*v^i internally)
  Qout(1:NVAR,is:ie,js:je) = Q(1:NVAR,is:ie,js:je)
  do j = js, je
  do i = is, ie
    u2_loc = Q(IVX,i,j)**2 + Q(IVY,i,j)**2 + Q(IVZ,i,j)**2
    W_loc  = dsqrt(1.0d0 + u2_loc)
    Qout(IVX,i,j) = Q(IVX,i,j) / W_loc
    Qout(IVY,i,j) = Q(IVY,i,j) / W_loc
    Qout(IVZ,i,j) = Q(IVZ,i,j) / W_loc
  enddo
  enddo

  write(filename,'(i5.5)') nsnap
  if (flag_binary) then
    filename = trim(dirname)//"/snap"//trim(filename)//".bin"
    open(unitsnap, file=filename, form='unformatted', access="stream", action="write")
    write(unitsnap) time
    write(unitsnap) nx
    write(unitsnap) ny
    write(unitsnap) NVAR
    write(unitsnap) NFLX-NVAR
    write(unitsnap) xv(is:ie)
    write(unitsnap) yv(js:je)
    write(unitsnap) real(Qout(1:NVAR,is:ie,js:je))  ! single precision, v^i
    write(unitsnap) real(Bc(1:3,is:ie,js:je))        ! single precision
    close(unitsnap)
  else
    filename = trim(dirname)//"/snap"//trim(filename)//".dat"
    open(unitsnap, file=filename, form='formatted', action="write")
    write(unitsnap,*) "# time = ", time
    write(unitsnap,*) "#nx, ny = ", nx, ny
    do j = js, je
    do i = is, ie
      write(unitsnap,'(1p,11(es24.16e3,1x))') xv(i), yv(j), Qout(IDN,i,j), Qout(IVX,i,j), Qout(IVY,i,j), Qout(IVZ,i,j), &
           Qout(IPR,i,j), Qout(ISC,i,j), Bc(1,i,j), Bc(2,i,j), Bc(3,i,j)
    enddo
    enddo
    close(unitsnap)
  endif

  write(6,'(A,A,1X,ES12.5)') "output file:", trim(filename), time

  nsnap = nsnap + 1
  tsnap = tsnap + dtsnap

return
end subroutine Output

!-------------------------------------------------------------------
!  Create directory
!-------------------------------------------------------------------
subroutine makedirs(outdir)
implicit none
character(len=*), intent(in) :: outdir
character(len=256) :: command
write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
call system(command)
end subroutine makedirs

!-------------------------------------------------------------------
!  Realtime analysis (v_y RMS and div B errors)
!-------------------------------------------------------------------
subroutine RealtimeAnalysis(xv, yv, Q, Bc, Bs, phys_evo)
use params, only : NVAR, nxtot, nytot, nx, ny, is, ie, js, je, &
                   IDN, IVX, IVY, IVZ, IPR, ISC, &
                   IMX, IMY, IMZ, IEN, gam, nevo
implicit none
real(8), intent(in)  :: xv(nxtot), yv(nytot)
real(8), intent(in)  :: Q(NVAR,nxtot,nytot), Bc(3,nxtot,nytot), Bs(3,nxtot,nytot)
real(8), intent(out) :: phys_evo(nevo)
integer :: i, j
real(8) :: dvy, er_divBc, er_divBs
real(8) :: u2_loc, vy_loc

  dvy      = 0.0d0
  er_divBc = 0.0d0
  er_divBs = 0.0d0
  do j = js, je
  do i = is, ie
    u2_loc = Q(IVX,i,j)**2 + Q(IVY,i,j)**2 + Q(IVZ,i,j)**2
    vy_loc = Q(IVY,i,j) / dsqrt(1.0d0 + u2_loc)   ! v^y = u^y / W
    dvy    = dvy + vy_loc**2
    er_divBs = er_divBs + ( Bs(1,i+1,j) - Bs(1,i,j) + Bs(2,i,j+1) - Bs(2,i,j) )**2
    er_divBc = er_divBc + 0.5d0*( Bc(1,i+1,j) - Bc(1,i-1,j) + Bc(2,i,j+1) - Bc(2,i,j-1) )**2
  enddo
  enddo
  phys_evo(1) = sqrt(dvy / dble(nx*ny))
  phys_evo(2) = sqrt(er_divBc / dble(nx*ny))
  phys_evo(3) = sqrt(er_divBs / dble(nx*ny))

  ! suppress unused warnings
  if (.false.) then
    i = 0; j = 0
    vy_loc = xv(1); vy_loc = yv(1)
  endif

return
end subroutine RealtimeAnalysis

!-------------------------------------------------------------------
!  Generate initial condition: SRMHD Kelvin-Helmholtz instability
!  Mattia & Mignone 2022, Fig.12 / Beckwith & Stone 2011
!-------------------------------------------------------------------
subroutine GenerateProblem(xv, yv, Q, Bs, Bc)
use params, only : NVAR, nxtot, nytot, &
                   is, ie, js, je, IDN, IVX, IVY, IVZ, IPR, ISC
use prim_mod, only : CellCenterMagneticField
implicit none
integer :: i, j
real(8), intent(in)  :: xv(nxtot), yv(nytot)
real(8), intent(out) :: Q(NVAR,nxtot,nytot)
real(8), intent(out) :: Bs(3,nxtot,nytot)
real(8), intent(out) :: Bc(3,nxtot,nytot)
real(8) :: pi, vsh, a_shear, A0, sigma, rho_h, rho_l, B0
real(8) :: sgny, arg_shear, arg_sig, vx, vy, rho, W_loc

  pi      = dacos(-1.0d0)
  vsh     = 0.5d0
  a_shear = 0.01d0
  A0      = 0.1d0
  sigma   = 0.1d0
  rho_h   = 1.0d0
  rho_l   = 0.01d0
  B0      = 1.0d-3

  do j = js, je
  do i = is, ie
    sgny      = dsign(1.0d0, yv(j))
    arg_shear = (2.0d0*yv(j) - sgny) / (2.0d0*a_shear)
    arg_sig   = (2.0d0*yv(j) - sgny) / (2.0d0*sigma)

    vx  = sgny * vsh * dtanh(arg_shear)
    vy  = sgny * A0 * vsh * dsin(2.0d0*pi*xv(i)) * dexp(-arg_sig**2)
    rho = 0.5d0*(rho_l + rho_h) + 0.5d0*(rho_h - rho_l)*(vx/vsh)
    W_loc = 1.0d0 / dsqrt(1.0d0 - vx**2 - vy**2)

    Q(IDN,i,j) = rho
    Q(IVX,i,j) = W_loc * vx    ! u^x = W*v^x
    Q(IVY,i,j) = W_loc * vy    ! u^y = W*v^y
    Q(IVZ,i,j) = 0.0d0
    Q(IPR,i,j) = 1.0d0
    Q(ISC,i,j) = (rho - rho_l) / (rho_h - rho_l)
  enddo
  enddo

  do j = js, je
  do i = is, ie+1
    Bs(1,i,j) = B0
  enddo
  enddo

  do j = js, je+1
  do i = is, ie
    Bs(2,i,j) = 0.0d0
  enddo
  enddo

  do j = js, je
  do i = is, ie
    Bs(3,i,j) = 0.0d0
  enddo
  enddo

  call CellCenterMagneticField(is, ie, js, je, Bs, Bc)

return
end subroutine GenerateProblem

end program main2d
