!---------------------------------------------------------------------
!  1D SRMHD test driver
!  Mattia & Mignone 2022, Table 1 test problems
!
!  Usage: ./main1d.x <flag_problem>
!    1 = CW  isolated contact wave        (gamma=5/3, t=1.0,  N=40)
!    2 = RW  isolated rotational wave     (gamma=5/3, t=1.0,  N=40)
!    3 = ST1 Brio-Wu shock tube           (gamma=2,   t=0.4,  N=400)
!    4 = ST2 generic shock tube 2         (gamma=5/3, t=0.55, N=800)
!    5 = ST3 highly relativistic          (gamma=5/3, t=0.4,  N=400)
!    6 = ST4 generic shock tube 4         (gamma=5/3, t=0.55, N=800)
!
!  Output: <casename>/snapNNNNN.dat
!    Columns: x, rho, vx, vy, vz, p, Bx, By, Bz, sc
!
!  Numerics: PLM + van Leer limiter + HLL + RK2 (midpoint method)
!  Q(IVX/Y/Z) stores 4-velocity u^i = W*v^i internally (same as 2D)
!---------------------------------------------------------------------
program main1d
use params,      only: NVAR, NFLX, IDN, IVX, IVY, IVZ, IPR, ISC, &
                       IBX, IBY, IBZ, IEN, IMX, IMY, IMZ
use flux_mod,    only: compute_dt, apply_bc_1d, reconstruct_1d, compute_flux_1d, update_1d
use prim_mod,    only: prim2consv_cell, consv2prim_cell
implicit none

integer, parameter :: ngh = 2   ! number of ghost cells

!--- Runtime parameters (set from flag_problem) ---
integer       :: flag_problem
integer       :: n1d, ntot, is, ie
real(8)       :: gam1d, tmax, dtsnap, cfl
real(8)       :: xmin, xmax, xdis, dx
character(20) :: outdir

!--- Left/right initial condition states ---
real(8) :: rho_L, p_L, vx_L, vy_L, vz_L
real(8) :: rho_R, p_R, vx_R, vy_R, vz_R
real(8) :: Bx_val, By_L, Bz_L, By_R, Bz_R

!--- Allocatable arrays ---
real(8), allocatable :: Q(:,:)      ! primitive (u^i for velocity)
real(8), allocatable :: Uo(:,:)     ! saved primitives for RK2
real(8), allocatable :: U(:,:)      ! conservative (workspace)
real(8), allocatable :: Bf(:,:)     ! cell-center B: (1)=Bx,(2)=By,(3)=Bz
real(8), allocatable :: Bfo(:,:)    ! saved B for RK2
real(8), allocatable :: Ql_arr(:,:) ! PLM left  state at each face
real(8), allocatable :: Qr_arr(:,:) ! PLM right state at each face
real(8), allocatable :: Flux(:,:)   ! numerical flux at face i (at x_{i-1/2})
real(8), allocatable :: xf(:)       ! face positions
real(8), allocatable :: xv(:)       ! cell-center positions

!--- Time integration ---
real(8)  :: dt, time, tsnap_next
integer  :: ntime, nsnap, i

!--- Command-line argument ---
character(20) :: argstr

!======================================================================
!  Parse command-line argument
!======================================================================
if (command_argument_count() < 1) then
  write(6,'(A)') 'Usage: ./main1d.x <problem> [outdir]'
  write(6,'(A)') '  1=CW  2=RW  3=ST1  4=ST2  5=ST3  6=ST4'
  stop
endif
call get_command_argument(1, argstr)
read(argstr,*) flag_problem

!======================================================================
!  Set problem parameters (Mattia & Mignone 2022, Table 1)
!======================================================================
cfl  = 0.4d0
xmin = 0.0d0;  xmax = 1.0d0;  xdis = 0.5d0

select case(flag_problem)

case(1)   !--- CW: isolated contact wave ---
  n1d    = 40;   gam1d = 5.0d0/3.0d0;  tmax = 1.0d0;  dtsnap = 0.5d0
  outdir = 'cw'
  rho_L = 10.0d0;  p_L = 1.0d0
  vx_L  =  0.0d0;  vy_L =  0.7d0;  vz_L = 0.2d0
  rho_R =  1.0d0;  p_R = 1.0d0
  vx_R  =  0.0d0;  vy_R =  0.7d0;  vz_R = 0.2d0
  Bx_val = 5.0d0
  By_L =  1.0d0;  Bz_L = 0.5d0
  By_R =  1.0d0;  Bz_R = 0.5d0

case(2)   !--- RW: isolated rotational wave ---
  n1d    = 40;   gam1d = 5.0d0/3.0d0;  tmax = 1.0d0;  dtsnap = 0.5d0
  outdir = 'rw'
  rho_L = 1.0d0;  p_L = 1.0d0
  vx_L  =  0.4d0;      vy_L = -0.3d0;      vz_L =  0.5d0
  rho_R = 1.0d0;  p_R = 1.0d0
  vx_R  =  0.377237d0; vy_R = -0.482389d0;  vz_R =  0.424190d0
  Bx_val = 2.4d0
  By_L =  1.0d0;  Bz_L = -1.6d0
  By_R = -0.1d0;  Bz_R = -2.178213d0

case(3)   !--- ST1: Brio-Wu relativistic shock tube (gamma=2) ---
  n1d    = 400;  gam1d = 2.0d0;          tmax = 0.4d0;  dtsnap = 0.1d0
  outdir = 'st1'
  rho_L = 1.0d0;    p_L = 1.0d0
  vx_L  = 0.0d0;  vy_L = 0.0d0;  vz_L = 0.0d0
  rho_R = 0.125d0;  p_R = 0.1d0
  vx_R  = 0.0d0;  vy_R = 0.0d0;  vz_R = 0.0d0
  Bx_val = 0.5d0
  By_L =  1.0d0;  Bz_L = 0.0d0
  By_R = -1.0d0;  Bz_R = 0.0d0

case(4)   !--- ST2: generic shock tube 2 ---
  n1d    = 800;  gam1d = 5.0d0/3.0d0;  tmax = 0.55d0; dtsnap = 0.1d0
  outdir = 'st2'
  rho_L = 1.08d0;  p_L = 0.95d0
  vx_L  =  0.4d0;  vy_L =  0.3d0;  vz_L = 0.2d0
  rho_R =  1.0d0;  p_R =  1.0d0
  vx_R  = -0.45d0; vy_R = -0.2d0;  vz_R = 0.2d0
  Bx_val = 2.0d0
  By_L =  0.3d0;  Bz_L = 0.3d0
  By_R = -0.7d0;  Bz_R = 0.5d0

case(5)   !--- ST3: highly relativistic ---
  n1d    = 400;  gam1d = 5.0d0/3.0d0;  tmax = 0.4d0;  dtsnap = 0.1d0
  outdir = 'st3'
  rho_L = 1.0d0;  p_L = 0.1d0
  vx_L  =  0.999d0;  vy_L = 0.0d0;  vz_L = 0.0d0
  rho_R = 1.0d0;  p_R = 0.1d0
  vx_R  = -0.999d0;  vy_R = 0.0d0;  vz_R = 0.0d0
  Bx_val = 10.0d0
  By_L =  7.0d0;  Bz_L = 7.0d0
  By_R = -7.0d0;  Bz_R = 7.0d0

case(6)   !--- ST4: generic shock tube 4 ---
  n1d    = 800;  gam1d = 5.0d0/3.0d0;  tmax = 0.55d0; dtsnap = 0.1d0
  outdir = 'st4'
  rho_L = 1.0d0;  p_L = 5.0d0
  vx_L  = 0.0d0;  vy_L = 0.3d0;  vz_L = 0.4d0
  rho_R = 0.9d0;  p_R = 5.3d0
  vx_R  = 0.0d0;  vy_R = 0.0d0;  vz_R = 0.0d0
  Bx_val = 1.0d0
  By_L = 6.0d0;  Bz_L = 2.0d0
  By_R = 5.0d0;  Bz_R = 2.0d0

case default
  write(6,*) 'Unknown flag_problem:', flag_problem
  stop
end select

! Optional second argument overrides the output directory
if (command_argument_count() >= 2) then
  call get_command_argument(2, outdir)
endif

ntot = n1d + 2*ngh
is   = ngh + 1
ie   = n1d + ngh

!======================================================================
!  Allocate arrays
!======================================================================
allocate(Q(NVAR,ntot), Uo(NVAR,ntot), U(NVAR,ntot))
allocate(Bf(3,ntot), Bfo(3,ntot))
allocate(Ql_arr(NFLX,ntot), Qr_arr(NFLX,ntot), Flux(NFLX,ntot))
allocate(xf(ntot), xv(ntot))

!======================================================================
!  Grid setup
!======================================================================
dx = (xmax - xmin) / n1d
do i = 1, ntot
  xf(i) = xmin + (i - ngh - 1) * dx
  xv(i) = xf(i) + 0.5d0 * dx
enddo

!======================================================================
!  Initial conditions and first output
!======================================================================
call setup_ic()
call system('mkdir -p '//trim(outdir))
nsnap      = 0
tsnap_next = 0.0d0
call output_1d(nsnap, 0.0d0)
nsnap      = nsnap + 1
tsnap_next = tsnap_next + dtsnap

!======================================================================
!  Time integration: RK2 (midpoint method)
!======================================================================
time  = 0.0d0
ntime = 0

do
  dt = compute_dt(Q, Bf, is, ie, dx, gam1d, cfl)
  if (time + dt > tmax) dt = tmax - time

  !--- Save state ---
  Uo  = Q
  Bfo = Bf

  !--- Stage 1: flux from Q0 -> update by dt/2 -> Q_half ---
  call reconstruct_1d(Q, Bf, Ql_arr, Qr_arr, is, ie)
  if (ntime == 0) call diag_compare_flux(Ql_arr, Qr_arr, is, ie, Bx_val, dt, dx, gam1d)
  call compute_flux_1d(Ql_arr, Qr_arr, Flux, is, ie, Bx_val, dt, dx, gam1d)
  call update_1d(Q, Uo, Bf, Bfo, Flux, is, ie, dx, Bx_val, gam1d, 0.5d0*dt)
  call apply_bc_1d(Q, Bf, ngh, is, ie)

  !--- Stage 2: flux from Q_half -> update Uo by dt -> Q_new ---
  call reconstruct_1d(Q, Bf, Ql_arr, Qr_arr, is, ie)
  call compute_flux_1d(Ql_arr, Qr_arr, Flux, is, ie, Bx_val, dt, dx, gam1d)
  call update_1d(Q, Uo, Bf, Bfo, Flux, is, ie, dx, Bx_val, gam1d, dt)
  call apply_bc_1d(Q, Bf, ngh, is, ie)

  time  = time  + dt
  ntime = ntime + 1

  if (mod(ntime,50) == 0) write(6,'(A,I0,A,ES12.5,A,ES12.5)') &
    'ntime=', ntime, '  time=', time, '  dt=', dt

  if (time >= tsnap_next - 1.0d-14) then
    call output_1d(nsnap, time)
    nsnap      = nsnap + 1
    tsnap_next = tsnap_next + dtsnap
  endif

  if (time >= tmax) exit
enddo

write(6,'(A,I0,A,F20.17)') ' main1d finished, ntime=', ntime, '  time=', time

deallocate(Q, Uo, U, Bf, Bfo, Ql_arr, Qr_arr, Flux, xf, xv)

!=======================================================================
contains
!=======================================================================

!-----------------------------------------------------------------------
!  Diagnostic: compare HLL vs HLLC flux at each face (ntime==0 only)
!-----------------------------------------------------------------------
subroutine diag_compare_flux(Ql_arr, Qr_arr, is_d, ie_d, Bx_val_d, dt_d, dx_d, gam_d)
use params,   only: NFLX, IDN
use riemann_mod, only: HLL, HLLC
implicit none
real(8), intent(in) :: Ql_arr(NFLX,*), Qr_arr(NFLX,*)
integer, intent(in) :: is_d, ie_d
real(8), intent(in) :: Bx_val_d, dt_d, dx_d, gam_d
real(8) :: flx_hll(NFLX), flx_hllc(NFLX), wgt
integer :: k
real(8) :: dFD, dFmx, dFE

write(6,'(A)') '=== diag: HLL vs HLLC flux comparison (step 1, stage 1) ==='
write(6,'(A6,4A14)') 'face', 'FD_hll', 'FD_hllc', 'dFD', 'dFmx'
do k = is_d, ie_d+1
  call HLL (1, Ql_arr(1,k), Qr_arr(1,k), Bx_val_d, dt_d, dx_d, flx_hll,  wgt, gam_d)
  call HLLC(1, Ql_arr(1,k), Qr_arr(1,k), Bx_val_d, dt_d, dx_d, flx_hllc, wgt, gam_d)
  dFD  = flx_hllc(IDN)   - flx_hll(IDN)
  dFmx = flx_hllc(IDN+1) - flx_hll(IDN+1)
  if (abs(dFD) > 1.0d-6 .or. abs(dFmx) > 1.0d-6) then
    write(6,'(I6,4ES14.5)') k, flx_hll(IDN), flx_hllc(IDN), dFD, dFmx
  endif
enddo
write(6,'(A)') '=== end diag ==='
end subroutine diag_compare_flux

!-----------------------------------------------------------------------
!  Initial conditions
!-----------------------------------------------------------------------
subroutine setup_ic()
implicit none
integer :: k
real(8) :: W_loc

Q  = 0.0d0
Bf = 0.0d0

do k = is, ie
  Bf(1,k) = Bx_val
  if (xv(k) < xdis) then
    W_loc     = 1.0d0 / dsqrt(1.0d0 - vx_L**2 - vy_L**2 - vz_L**2)
    Q(IDN,k)  = rho_L
    Q(IVX,k)  = W_loc * vx_L
    Q(IVY,k)  = W_loc * vy_L
    Q(IVZ,k)  = W_loc * vz_L
    Q(IPR,k)  = p_L
    Q(ISC,k)  = 1.0d0
    Bf(2,k)   = By_L;  Bf(3,k) = Bz_L
  else
    W_loc     = 1.0d0 / dsqrt(1.0d0 - vx_R**2 - vy_R**2 - vz_R**2)
    Q(IDN,k)  = rho_R
    Q(IVX,k)  = W_loc * vx_R
    Q(IVY,k)  = W_loc * vy_R
    Q(IVZ,k)  = W_loc * vz_R
    Q(IPR,k)  = p_R
    Q(ISC,k)  = 0.0d0
    Bf(2,k)   = By_R;  Bf(3,k) = Bz_R
  endif
enddo
call apply_bc_1d(Q, Bf, ngh, is, ie)
end subroutine setup_ic

!-----------------------------------------------------------------------
!  Snapshot output (converts u^i -> v^i for writing)
!-----------------------------------------------------------------------
subroutine output_1d(snap_idx, time_out)
implicit none
integer, intent(in) :: snap_idx
real(8), intent(in) :: time_out
integer  :: k, unit_out
character(100) :: filename
real(8)  :: u2_loc, W_loc

write(filename,'(A,A,I5.5,A)') trim(outdir), '/snap', snap_idx, '.dat'
unit_out = 20
open(unit_out, file=trim(filename), form='formatted', action='write')
write(unit_out,'(A,ES14.6)') '# time = ', time_out
write(unit_out,'(A,I0)')     '# nx = ', n1d
write(unit_out,'(A)')        '# x  rho  vx  vy  vz  p  Bx  By  Bz  sc'
do k = is, ie
  u2_loc = Q(IVX,k)**2 + Q(IVY,k)**2 + Q(IVZ,k)**2
  W_loc  = dsqrt(1.0d0 + u2_loc)
  write(unit_out,'(1p,10(ES22.12E3,1x))') &
    xv(k), Q(IDN,k), Q(IVX,k)/W_loc, Q(IVY,k)/W_loc, Q(IVZ,k)/W_loc, &
    Q(IPR,k), Bf(1,k), Bf(2,k), Bf(3,k), Q(ISC,k)
enddo
close(unit_out)
write(6,'(A,A,A,ES12.5)') 'output: ', trim(filename), '  t=', time_out
end subroutine output_1d

end program main1d
