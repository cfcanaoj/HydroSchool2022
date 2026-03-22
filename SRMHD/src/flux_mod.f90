module flux_mod
implicit none
contains

!-------------------------------------------------------------------
!  2D: timestep control
!-------------------------------------------------------------------
real(8) function TimestepControl(xf, yf, Q, Bc)
use params, only : nxtot, nytot, NVAR, is, ie, js, je, &
                   IDN, IVX, IVY, IVZ, IPR, gam, cfl_number
implicit none
real(8), intent(in) :: xf(nxtot), yf(nytot)
real(8), intent(in) :: Q(NVAR, nxtot, nytot), Bc(3, nxtot, nytot)
real(8) :: dtl1, dtl2, dtmin
real(8) :: rho, vx, vy, vz, p_prim, Bx, By, Bz
real(8) :: v2, W2, rho_h, h, b2, vB, cs2, vA2, cf2, cf
real(8) :: lam_x, lam_y
integer :: i, j

dtmin = 1.0d90

!$omp parallel do default(none) collapse(2) schedule(static) &
!$omp shared(xf, yf, Q, Bc) &
!$omp private(i, j, rho, vx, vy, vz, p_prim, Bx, By, Bz, &
!$omp         v2, W2, rho_h, h, b2, vB, cs2, vA2, cf2, cf, &
!$omp         lam_x, lam_y, dtl1, dtl2) &
!$omp reduction(min:dtmin)
do j = js, je
do i = is, ie
  rho    = Q(IDN, i, j)
  vx     = Q(IVX, i, j)    ! u^x = W*v^x
  vy     = Q(IVY, i, j)    ! u^y = W*v^y
  vz     = Q(IVZ, i, j)    ! u^z = W*v^z
  p_prim = Q(IPR, i, j)
  Bx     = Bc(1, i, j)
  By     = Bc(2, i, j)
  Bz     = Bc(3, i, j)

  v2  = vx**2 + vy**2 + vz**2   ! u^2
  W2  = 1.0d0 + v2               ! W^2 = 1 + u^2
  vx  = vx / dsqrt(W2);  vy = vy / dsqrt(W2);  vz = vz / dsqrt(W2)
  v2  = v2 / W2                  ! now v^2
  h   = 1.0d0 + gam / (gam - 1.0d0) * p_prim / rho
  rho_h = rho * h

  vB  = vx*Bx + vy*By + vz*Bz
  b2  = (Bx**2 + By**2 + Bz**2) / W2 + vB**2

  cs2 = gam * p_prim / rho_h
  vA2 = b2 / (rho_h + b2)

  cf2 = cs2 + vA2*(1.0d0 - cs2)
  cf  = dsqrt(cf2)

  lam_x = (dabs(vx) + cf) / (1.0d0 + dabs(vx)*cf)
  lam_y = (dabs(vy) + cf) / (1.0d0 + dabs(vy)*cf)

  dtl1 = (xf(i+1) - xf(i)) / lam_x
  dtl2 = (yf(j+1) - yf(j)) / lam_y
  dtmin = min(dtl1, dtl2, dtmin)
enddo
enddo
!$omp end parallel do

TimestepControl = cfl_number * dtmin

return
end function TimestepControl

!-------------------------------------------------------------------
!  2D: PLM reconstruction + Riemann fluxes (all faces, both directions)
!-------------------------------------------------------------------
subroutine NumericalFlux(dt, xf, yf, Q, Bc, Bs, F, G, E)
use params, only : nxtot, nytot, NFLX, NVAR, &
                   is, ie, js, je, &
                   IDN, IPR, IVX, IVY, IVZ, IBX, IBY, IBZ, flag_flux, gam
use riemann_mod, only : HLL, HLLC, HLLD, vanLeer
implicit none
real(8), intent(in)  :: dt, xf(nxtot), yf(nytot)
real(8), intent(in)  :: Q(NVAR,nxtot,nytot)
real(8), intent(in)  :: Bc(3,nxtot,nytot)
real(8), intent(in)  :: Bs(3,nxtot,nytot)
real(8), intent(out) :: F(NVAR,nxtot,nytot)
real(8), intent(out) :: G(NVAR,nxtot,nytot)
real(8), intent(out) :: E(3,nxtot,nytot)

! the following arrays are shared between threads
real(8), save :: Ql(NFLX,nxtot,nytot), Qr(NFLX,nxtot,nytot)
real(8), save :: weight1(nxtot,nytot), weight2(nxtot,nytot)
real(8), save :: e2_xf(nxtot,nytot), e3_xf(nxtot,nytot)
real(8), save :: e1_yf(nxtot,nytot), e3_yf(nxtot,nytot)
real(8), save :: Etmp(nxtot,nytot)

integer :: i, j
real(8) :: flx(NFLX)
real(8) :: dQm(NVAR), dQp(NVAR), dQmon(NVAR)
real(8) :: wghtCT
real(8) :: de3_l1, de3_r1, de3_l2, de3_r2
real(8) :: Ql_c(NFLX), Qr_c(NFLX)

! numerical flux in the x direction

! hydro part
!$omp do collapse(2) schedule(static) private(i,j,dQp,dQm,dQmon)
  do j = js-1, je+1
  do i = is-1, ie+1
    dQp(1:NVAR) = Q(1:NVAR,i+1,j) - Q(1:NVAR,i,j)
    dQm(1:NVAR) = Q(1:NVAR,i,j)   - Q(1:NVAR,i-1,j)

    call vanLeer(NVAR, dQp, dQm, dQmon)

    Ql(1:NVAR,i+1,j) = Q(1:NVAR,i,j) + 0.5d0*dQmon(1:NVAR)
    Qr(1:NVAR,i  ,j) = Q(1:NVAR,i,j) - 0.5d0*dQmon(1:NVAR)
  enddo
  enddo
!$omp end do

! B field part
!$omp do collapse(2) schedule(static) private(i,j,dQp,dQm,dQmon)
  do j = js-1, je+1
  do i = is-1, ie+1
    dQp(1:3) = Bc(1:3,i+1,j) - Bc(1:3,i  ,j)
    dQm(1:3) = Bc(1:3,i  ,j) - Bc(1:3,i-1,j)

    call vanLeer(3, dQp, dQm, dQmon)

    Ql(NVAR+1:NFLX,i+1,j) = Bc(1:3,i,j) + 0.5d0*dQmon(1:3)
    Qr(NVAR+1:NFLX,i  ,j) = Bc(1:3,i,j) - 0.5d0*dQmon(1:3)
  enddo
  enddo
!$omp end do

!$omp do collapse(2) schedule(static) private(i,j,flx,wghtCT,Ql_c,Qr_c)
  do j = js-1, je+1
  do i = is, ie+1
    Ql_c(1:NFLX) = Ql(1:NFLX, i, j)
    Qr_c(1:NFLX) = Qr(1:NFLX, i, j)
    if (flag_flux == 1) then
      call HLL (1, Ql_c, Qr_c, Bs(1,i,j), dt, xf(i+1)-xf(i), flx, wghtCT, gam)
    else if (flag_flux == 2) then
      call HLLC(1, Ql_c, Qr_c, Bs(1,i,j), dt, xf(i+1)-xf(i), flx, wghtCT, gam)
    else
      call HLLD(1, Ql_c, Qr_c, Bs(1,i,j), dt, xf(i+1)-xf(i), flx, wghtCT, gam)
    endif

    F(1:NVAR,i,j)  = flx(1:NVAR)
    e3_xf(i,j) =  -flx(IBY)
    e2_xf(i,j) =  +flx(IBZ)

    weight1(i,j) = wghtCT
  enddo
  enddo
!$omp end do


! numerical flux in the y direction
!$omp do collapse(2) schedule(static) private(i,j,dQp,dQm,dQmon)
  do j = js-1, je+1
  do i = is-1, ie+1
    dQp(1:NVAR) = Q(1:NVAR,i,j+1) - Q(1:NVAR,i,j)
    dQm(1:NVAR) = Q(1:NVAR,i,j)   - Q(1:NVAR,i,j-1)

    call vanLeer(NVAR, dQp, dQm, dQmon)

    Ql(1:NVAR,i,j+1) = Q(1:NVAR,i,j) + 0.5d0*dQmon(1:NVAR)
    Qr(1:NVAR,i,j  ) = Q(1:NVAR,i,j) - 0.5d0*dQmon(1:NVAR)
  enddo
  enddo
!$omp end do

! B field part
!$omp do collapse(2) schedule(static) private(i,j,dQp,dQm,dQmon)
  do j = js-1, je+1
  do i = is-1, ie+1
    dQp(1:3) = Bc(1:3,i,j+1) - Bc(1:3,i,j  )
    dQm(1:3) = Bc(1:3,i,j  ) - Bc(1:3,i,j-1)

    call vanLeer(3, dQp, dQm, dQmon)

    Ql(NVAR+1:NFLX,i,j+1) = Bc(1:3,i,j) + 0.5d0*dQmon(1:3)
    Qr(NVAR+1:NFLX,i,j  ) = Bc(1:3,i,j) - 0.5d0*dQmon(1:3)
  enddo
  enddo
!$omp end do

!$omp do collapse(2) schedule(static) private(i,j,flx,wghtCT,Ql_c,Qr_c)
  do j = js, je+1
  do i = is-1, ie+1
    Ql_c(1:NFLX) = Ql(1:NFLX, i, j)
    Qr_c(1:NFLX) = Qr(1:NFLX, i, j)
    if (flag_flux == 1) then
      call HLL (2, Ql_c, Qr_c, Bs(2,i,j), dt, yf(j+1)-yf(j), flx, wghtCT, gam)
    else if (flag_flux == 2) then
      call HLLC(2, Ql_c, Qr_c, Bs(2,i,j), dt, yf(j+1)-yf(j), flx, wghtCT, gam)
    else
      call HLLD(2, Ql_c, Qr_c, Bs(2,i,j), dt, yf(j+1)-yf(j), flx, wghtCT, gam)
    endif

    G(1:NVAR,i,j) = flx(1:NVAR)
    e1_yf(i,j) =  - flx(IBZ)
    e3_yf(i,j) =    flx(IBX)

    weight2(i,j) = wghtCT
  enddo
  enddo
!$omp end do


!$omp do collapse(2) schedule(static) private(i,j)
  do j = js-1, je+1
  do i = is-1, ie+1
    Etmp(i,j) = Q(IVY,i,j)*Bc(1,i,j) - Q(IVX,i,j)*Bc(2,i,j)
  enddo
  enddo
!$omp end do


!$omp do collapse(2) schedule(static) private(i,j)
  do j = js, je
  do i = is, ie+1
    E(2,i,j) = e2_xf(i,j)
  enddo
  enddo
!$omp end do

!$omp do collapse(2) schedule(static) private(i,j)
  do j = js, je+1
  do i = is, ie
    E(1,i,j) = e1_yf(i,j)
  enddo
  enddo
!$omp end do

!$omp do collapse(2) schedule(static) private(de3_l2, de3_r2, de3_l1, de3_r1, i, j)
  do j = js, je+1
  do i = is, ie+1
    de3_l2 = (1.0d0-weight1(i,j-1))*(e3_yf(i  ,j) - Etmp(i  ,j-1)) + &
             (      weight1(i,j-1))*(e3_yf(i-1,j) - Etmp(i-1,j-1))

    de3_r2 = (1.0d0-weight1(i,j  ))*(e3_yf(i  ,j) - Etmp(i  ,j  )) + &
             (      weight1(i,j  ))*(e3_yf(i-1,j) - Etmp(i-1,j  ))

    de3_l1 = (1.0d0-weight2(i-1,j))*(e3_xf(i,j  ) - Etmp(i-1,j  )) + &
             (      weight2(i-1,j))*(e3_xf(i,j-1) - Etmp(i-1,j-1))

    de3_r1 = (1.0d0-weight2(i  ,j))*(e3_xf(i,j  ) - Etmp(i  ,j  )) + &
             (      weight2(i  ,j))*(e3_xf(i,j-1) - Etmp(i  ,j-1))

    E(3,i,j) = 0.25d0*( de3_l1 + de3_r1 + de3_l2 + de3_r2 + &
                         e3_yf(i-1,j) + e3_yf(i,j) + e3_xf(i,j-1) + e3_xf(i,j))
  enddo
  enddo
!$omp end do

return
end subroutine NumericalFlux

!-------------------------------------------------------------------
!  2D: boundary conditions (periodic, all sides)
!-------------------------------------------------------------------
subroutine BoundaryCondition(Q, Bs, Bc)
use params, only : NVAR, nxtot, nytot, is, ie, js, je, ngh, &
                   IDN, IVX, IVY, IVZ, IPR
use prim_mod, only : CellCenterMagneticField
implicit none
real(8), intent(inout) :: Q(NVAR,nxtot,nytot)
real(8), intent(inout) :: Bs(3,nxtot,nytot)
real(8), intent(inout) :: Bc(3,nxtot,nytot)
integer :: i, j

!$omp do collapse(2) schedule(static) private(i,j)
  do j = js-ngh, je+ngh
  do i = 1, ngh
    Q(:,is-i,j)  = Q(:,ie+1-i,j)
    Q(:,ie+i,j)  = Q(:,is+i-1,j)
  enddo
  enddo
!$omp end do

!$omp do collapse(2) schedule(static) private(i,j)
  do j = js-ngh, je+ngh
  do i = 1, ngh
    Bs(1,is-i,j)   = Bs(1,ie+1-i,j)
    Bs(1,ie+i+1,j) = Bs(1,is+i,j)
  enddo
  enddo
!$omp end do

!$omp do collapse(2) schedule(static) private(i,j)
  do j = js-ngh, je+ngh+1
  do i = 1, ngh
    Bs(2,is-i,j) = Bs(2,ie+1-i,j)
    Bs(2,ie+i,j) = Bs(2,is+i-1,j)
  enddo
  enddo
!$omp end do

!$omp do collapse(2) schedule(static) private(i,j)
  do j = js-ngh, je+ngh
  do i = 1, ngh
    Bs(3,is-i,j) = Bs(3,ie+1-i,j)
    Bs(3,ie+i,j) = Bs(3,is+i-1,j)
  enddo
  enddo
!$omp end do

!$omp do collapse(2) schedule(static) private(i,j)
  ! y inner boundary
  do j = 1, ngh
  do i = is-ngh, ie+ngh
    Q(:,i,js-j) = Q(:,i,je+1-j)
    Q(:,i,je+j) = Q(:,i,js+j-1)
  enddo
  enddo
!$omp end do

!$omp do collapse(2) schedule(static) private(i,j)
  do j = 1, ngh
  do i = is-ngh, ie+ngh+1
    Bs(1,i,js-j) = Bs(1,i,je+1-j)
    Bs(1,i,je+j) = Bs(1,i,js+j-1)
  enddo
  enddo
!$omp end do

!$omp do collapse(2) schedule(static) private(i,j)
  do j = 1, ngh
  do i = is-ngh, ie+ngh
    Bs(2,i,js-j)   = Bs(2,i,je+1-j)
    Bs(2,i,je+j+1) = Bs(2,i,js+j)
  enddo
  enddo
!$omp end do

!$omp do collapse(2) schedule(static) private(i,j)
  do j = 1, ngh
  do i = is-ngh, ie+ngh
    Bs(3,i,js-j) = Bs(3,i,je+1-j)
    Bs(3,i,je+j) = Bs(3,i,js+j-1)
  enddo
  enddo
!$omp end do

  ! boundary condition for the cell centered B field
  call CellCenterMagneticField(is-ngh, is-1,   js-ngh, je+ngh, Bs, Bc)
  call CellCenterMagneticField(ie+1,   ie+ngh, js-ngh, je+ngh, Bs, Bc)
  call CellCenterMagneticField(is-ngh, ie+ngh, js-ngh, js-1,   Bs, Bc)
  call CellCenterMagneticField(is-ngh, ie+ngh, je+1,   je+ngh, Bs, Bc)

return
end subroutine BoundaryCondition

!-------------------------------------------------------------------
!  2D: update conservative variables U using numerical flux
!-------------------------------------------------------------------
subroutine UpdateConsv(dt1, xf, xv, yf, yv, F, G, E, Uo, Bso, U, Bs)
use params, only : NVAR, NFLX, nxtot, nytot, is, ie, js, je
implicit none
real(8), intent(in)  :: dt1
real(8), intent(in)  :: xf(nxtot), yf(nytot)
real(8), intent(in)  :: xv(nxtot), yv(nytot)
real(8), intent(in)  :: F(NVAR,nxtot,nytot), G(NVAR,nxtot,nytot), E(3,nxtot,nytot)
real(8), intent(in)  :: Uo(NVAR,nxtot,nytot), Bso(3,nxtot,nytot)
real(8), intent(out) :: U(NVAR,nxtot,nytot), Bs(3,nxtot,nytot)
integer :: i, j, ihy

!$omp do collapse(2) schedule(static) private(i,j,ihy)
  do j = js, je
  do i = is, ie
  do ihy = 1, NVAR
    U(ihy,i,j) = Uo(ihy,i,j) + dt1*(- F(ihy,i+1,j) + F(ihy,i,j))/(xf(i+1)-xf(i)) &
                              + dt1*(- G(ihy,i,j+1) + G(ihy,i,j))/(yf(j+1)-yf(j))
  enddo
  enddo
  enddo
!$omp end do

!$omp do collapse(2) schedule(static) private(i,j)
  do j = js, je
  do i = is, ie+1
    Bs(1,i,j) = Bso(1,i,j) &
              - dt1*(E(3,i,j+1) - E(3,i,j))/(yv(j+1) - yv(j))
  enddo
  enddo
!$omp end do

!$omp do collapse(2) schedule(static) private(i,j)
  do j = js, je+1
  do i = is, ie
    Bs(2,i,j) = Bso(2,i,j) &
              + dt1*(E(3,i+1,j) - E(3,i,j))/(xv(i+1) - xv(i))
  enddo
  enddo
!$omp end do

!$omp do collapse(2) schedule(static) private(i,j)
  do j = js, je
  do i = is, ie
    Bs(3,i,j) = Bso(3,i,j) &
              - dt1*(E(2,i+1,j) - E(2,i,j))/(xf(i+1) - xf(i)) &
              + dt1*(E(1,i,j+1) - E(1,i,j))/(yf(j+1) - yf(j))
  enddo
  enddo
!$omp end do

return
end subroutine UpdateConsv

!-------------------------------------------------------------------
!  2D: source terms (placeholder)
!-------------------------------------------------------------------
subroutine SrcTerms(dt1, Q, U)
use params, only : IDN, IVX, IVY, IVZ, IPR, IBX, IBY, IBZ, &
                   IMX, IMY, IMZ, IEN, NVAR, is, ie, js, je, gam
implicit none
real(8), intent(in)    :: dt1
real(8), intent(in)    :: Q(:,:,:)
real(8), intent(inout) :: U(:,:,:)
integer :: i, j

! suppress unused variable warnings
if (.false.) then
  i = 0; j = 0
endif

! Source term (currently empty placeholder)

end subroutine SrcTerms

!-------------------------------------------------------------------
!  1D: outflow boundary conditions
!-------------------------------------------------------------------
subroutine apply_bc_1d(Q, Bf, ngh, is, ie)
implicit none
integer, intent(in)    :: ngh, is, ie
real(8), intent(inout) :: Q(:,:)
real(8), intent(inout) :: Bf(:,:)
integer :: k

do k = 1, ngh
  Q(:, is-k)  = Q(:, is)
  Q(:, ie+k)  = Q(:, ie)
  Bf(:, is-k) = Bf(:, is)
  Bf(:, ie+k) = Bf(:, ie)
enddo
end subroutine apply_bc_1d

!-------------------------------------------------------------------
!  1D: CFL timestep
!-------------------------------------------------------------------
real(8) function compute_dt(Q, Bf, is, ie, dx, gam1d, cfl)
use params, only: IDN, IVX, IVY, IVZ, IPR
implicit none
integer, intent(in) :: is, ie
real(8), intent(in) :: Q(:,:), Bf(:,:)
real(8), intent(in) :: dx, gam1d, cfl
integer :: k
real(8) :: u2_loc, W_loc, vx_loc, rho_loc, p_loc
real(8) :: rho_h, vB_loc, b2_loc, cs2, vA2, cf, lam, dt_min

dt_min = 1.0d90
do k = is, ie
  u2_loc  = Q(IVX,k)**2 + Q(IVY,k)**2 + Q(IVZ,k)**2
  W_loc   = dsqrt(1.0d0 + u2_loc)
  vx_loc  = Q(IVX,k) / W_loc
  rho_loc = Q(IDN,k);  p_loc = Q(IPR,k)

  rho_h  = rho_loc * (1.0d0 + gam1d/(gam1d-1.0d0)*p_loc/rho_loc)
  vB_loc = (Q(IVX,k)*Bf(1,k) + Q(IVY,k)*Bf(2,k) + Q(IVZ,k)*Bf(3,k)) / W_loc
  b2_loc = (Bf(1,k)**2 + Bf(2,k)**2 + Bf(3,k)**2)/(1.0d0+u2_loc) + vB_loc**2
  cs2    = gam1d*p_loc / rho_h
  vA2    = b2_loc / (rho_h + b2_loc)
  cf     = dsqrt(cs2 + vA2*(1.0d0-cs2))
  lam    = (dabs(vx_loc)+cf) / (1.0d0+dabs(vx_loc)*cf)
  dt_min = min(dt_min, dx/lam)
enddo
compute_dt = cfl * dt_min
end function compute_dt

!-------------------------------------------------------------------
!  1D: PLM reconstruction (van Leer limiter)
!  Ql_arr(:,i) = left  state at face x_{i-1/2} (from cell i-1)
!  Qr_arr(:,i) = right state at face x_{i-1/2} (from cell i)
!-------------------------------------------------------------------
subroutine reconstruct_1d(Q, Bf, Ql_arr, Qr_arr, is, ie)
use params, only: NVAR, flag_recon
use riemann_mod, only: vanLeer
implicit none
integer, intent(in)  :: is, ie
real(8), intent(in)  :: Q(:,:), Bf(:,:)
real(8), intent(out) :: Ql_arr(:,:), Qr_arr(:,:)
integer :: k
real(8) :: dQp(NVAR), dQm(NVAR), dQmon(NVAR)

Ql_arr = 0.0d0
Qr_arr = 0.0d0

do k = is-1, ie+1
  if (flag_recon == 2) then
    ! PLM (2nd-order) with van Leer limiter
    dQp(1:NVAR) = Q(1:NVAR, k+1) - Q(1:NVAR, k  )
    dQm(1:NVAR) = Q(1:NVAR, k  ) - Q(1:NVAR, k-1)
    call vanLeer(NVAR, dQp, dQm, dQmon)
    Ql_arr(1:NVAR, k+1) = Q(1:NVAR, k) + 0.5d0*dQmon(1:NVAR)
    Qr_arr(1:NVAR, k  ) = Q(1:NVAR, k) - 0.5d0*dQmon(1:NVAR)

    dQp(1:3) = Bf(1:3, k+1) - Bf(1:3, k  )
    dQm(1:3) = Bf(1:3, k  ) - Bf(1:3, k-1)
    call vanLeer(3, dQp, dQm, dQmon)
    Ql_arr(NVAR+1:NVAR+3, k+1) = Bf(1:3, k) + 0.5d0*dQmon(1:3)
    Qr_arr(NVAR+1:NVAR+3, k  ) = Bf(1:3, k) - 0.5d0*dQmon(1:3)
  else
    ! Flat (1st-order) reconstruction
    Ql_arr(1:NVAR,       k+1) = Q(1:NVAR, k)
    Qr_arr(1:NVAR,       k  ) = Q(1:NVAR, k)
    Ql_arr(NVAR+1:NVAR+3,k+1) = Bf(1:3, k)
    Qr_arr(NVAR+1:NVAR+3,k  ) = Bf(1:3, k)
  endif
enddo
end subroutine reconstruct_1d

!-------------------------------------------------------------------
!  1D: Riemann flux at all faces (HLL or HLLC depending on flag_flux)
!-------------------------------------------------------------------
subroutine compute_flux_1d(Ql_arr, Qr_arr, Flux, is, ie, Bx_val, dt_arg, dx, gam1d)
use params, only: NFLX, flag_flux
use riemann_mod, only: HLL, HLLC
implicit none
integer, intent(in)  :: is, ie
real(8), intent(in)  :: Ql_arr(:,:), Qr_arr(:,:)
real(8), intent(out) :: Flux(:,:)
real(8), intent(in)  :: Bx_val, dt_arg, dx, gam1d
integer :: k
real(8) :: Ql_c(NFLX), Qr_c(NFLX), flx(NFLX), wghtCT

Flux = 0.0d0
do k = is, ie+1
  Ql_c(1:NFLX) = Ql_arr(1:NFLX, k)
  Qr_c(1:NFLX) = Qr_arr(1:NFLX, k)
  if (flag_flux == 2) then
    call HLLC(1, Ql_c, Qr_c, Bx_val, dt_arg, dx, flx, wghtCT, gam1d)
  else
    call HLL (1, Ql_c, Qr_c, Bx_val, dt_arg, dx, flx, wghtCT, gam1d)
  endif
  Flux(1:NFLX, k) = flx(1:NFLX)
enddo
end subroutine compute_flux_1d

!-------------------------------------------------------------------
!  1D: update conservatives from Uo, Bfo + recover primitives Q
!-------------------------------------------------------------------
subroutine update_1d(Q, Uo, Bf, Bfo, Flux, is, ie, dx, Bx_val, gam1d, dt_stage)
use params, only: NVAR, IBY, IBZ
use prim_mod, only: prim2consv_cell, consv2prim_cell
implicit none
integer, intent(in)    :: is, ie
real(8), intent(inout) :: Q(:,:), Bf(:,:)
real(8), intent(in)    :: Uo(:,:), Bfo(:,:), Flux(:,:)
real(8), intent(in)    :: dx, Bx_val, gam1d, dt_stage
integer :: k
real(8) :: Q_c(NVAR), U_c(NVAR), Bc_c(3), dtodx

dtodx = dt_stage / dx

do k = is, ie
  Q_c(1:NVAR) = Uo(1:NVAR, k)
  Bc_c(1) = Bx_val;  Bc_c(2) = Bfo(2,k);  Bc_c(3) = Bfo(3,k)
  call prim2consv_cell(Q_c, Bc_c, U_c, gam1d)

  U_c(1:NVAR) = U_c(1:NVAR) - dtodx*(Flux(1:NVAR,k+1) - Flux(1:NVAR,k))

  Bf(2,k) = Bfo(2,k) - dtodx*(Flux(IBY,k+1) - Flux(IBY,k))
  Bf(3,k) = Bfo(3,k) - dtodx*(Flux(IBZ,k+1) - Flux(IBZ,k))

  Bc_c(2) = Bf(2,k);  Bc_c(3) = Bf(3,k)
  call consv2prim_cell(U_c, Bc_c, Q_c, gam1d)
  Q(1:NVAR, k) = Q_c(1:NVAR)
enddo
end subroutine update_1d

end module flux_mod
