module prim_mod
implicit none
contains

!-------------------------------------------------------------------
!       Cell-level Primitive -> Conservative (SRMHD)
!       Q_c(IVX/Y/Z) stores u^i = W*v^i
!-------------------------------------------------------------------
subroutine prim2consv_cell(Q_c, Bc_c, U_c, gam_arg)
use params, only : NVAR, IDN, IVX, IVY, IVZ, IPR, ISC, &
                   IMX, IMY, IMZ, IEN
implicit none
real(8), intent(in)  :: Q_c(NVAR), Bc_c(3), gam_arg
real(8), intent(out) :: U_c(NVAR)
real(8) :: rho, ux, uy, uz, vx, vy, vz, p_prim, Bx, By, Bz
real(8) :: u2, W2, W, h, xi, vB, B2

  rho    = Q_c(IDN)
  ux     = Q_c(IVX)    ! u^x = W*v^x
  uy     = Q_c(IVY)    ! u^y = W*v^y
  uz     = Q_c(IVZ)    ! u^z = W*v^z
  p_prim = Q_c(IPR)
  Bx     = Bc_c(1)
  By     = Bc_c(2)
  Bz     = Bc_c(3)

  u2 = ux**2 + uy**2 + uz**2
  W2 = 1.0d0 + u2
  W  = dsqrt(W2)
  vx = ux / W;  vy = uy / W;  vz = uz / W
  h  = 1.0d0 + gam_arg / (gam_arg - 1.0d0) * p_prim / rho

  B2 = Bx**2 + By**2 + Bz**2
  vB = vx*Bx + vy*By + vz*Bz

  xi = rho*h*W2 + B2    ! = rho*h*W^2 + B^2

  U_c(IDN) = rho * W
  U_c(IMX) = xi*vx - vB*Bx
  U_c(IMY) = xi*vy - vB*By
  U_c(IMZ) = xi*vz - vB*Bz
  U_c(IEN) = rho*h*W2 - p_prim - rho*W &
           + 0.5d0*B2*(1.0d0 + vx**2+vy**2+vz**2) - 0.5d0*vB**2
  U_c(ISC) = rho*W * Q_c(ISC)    ! D * sc

return
end subroutine prim2consv_cell

!-------------------------------------------------------------------
!       Cell-level Conservative -> Primitive (SRMHD, Newton-Raphson)
!       Q_c(IVX/Y/Z) stores u^i = W*v^i
!-------------------------------------------------------------------
subroutine consv2prim_cell(U_c, Bc_c, Q_c, gam_arg)
use params, only : NVAR, IDN, IVX, IVY, IVZ, IPR, ISC, &
                   IMX, IMY, IMZ, IEN
implicit none
real(8), intent(in)  :: U_c(NVAR), Bc_c(3), gam_arg
real(8), intent(out) :: Q_c(NVAR)
integer :: iter
real(8) :: D_con, Sx, Sy, Sz, E_con
real(8) :: Bx, By, Bz
real(8) :: Smag2, B2, SB
real(8) :: xi, aa, omv2, omv, vB, p_new
real(8) :: nW2, dW2, W2_loc              ! numerator/denominator/value of W²
real(8) :: dv2dxi, dpdxi, f_xi, df_xi
integer, parameter :: itermax = 50
real(8), parameter :: tol = 1.0d-10

  D_con = U_c(IDN)
  Sx    = U_c(IMX)
  Sy    = U_c(IMY)
  Sz    = U_c(IMZ)
  E_con = U_c(IEN)
  Bx    = Bc_c(1)
  By    = Bc_c(2)
  Bz    = Bc_c(3)

  Smag2 = Sx**2 + Sy**2 + Sz**2
  B2    = Bx**2 + By**2 + Bz**2
  SB    = Sx*Bx + Sy*By + Sz*Bz

  ! Initial guess: E+D+B²/2 avoids the superluminal region
  ! for magnetised states where E+D alone undershoots xi_min
  xi = E_con + D_con + 0.5d0*B2

  ! Newton-Raphson: solve f(xi) = 0
  do iter = 1, itermax
    aa = xi - B2

    ! Compute W directly:  W² = (ξ·aa)² / [(ξ·aa)² − S²·aa² − (S·B)²·(aa+ξ)]
    ! This avoids ever taking sqrt(1-v²); W ≥ 1 is enforced naturally.
    nW2 = (xi * aa)**2
    dW2 = nW2 - Smag2*aa**2 - SB**2*(aa + xi)

    if (dW2 <= 0.0d0) then
      ! Superluminal: ξ is below the physical minimum.  Push upward.
      xi = xi * 1.5d0 + B2 * 0.01d0
      xi = max(xi, B2*(1.0d0 + 1.0d-12))
      cycle
    endif

    W2_loc = nW2 / dW2
    W2_loc = max(W2_loc, 1.0d0 + 1.0d-15)   ! W ≥ 1
    omv2   = 1.0d0 / W2_loc                  ! = 1 − v²
    omv    = 1.0d0 / dsqrt(W2_loc)           ! = 1/W  (never NaN)

    p_new = (gam_arg - 1.0d0) / gam_arg * (aa*omv2 - D_con*omv)

    f_xi = xi - 0.5d0*B2*omv2 - 0.5d0*SB**2/aa**2 - p_new - D_con - E_con

    if (dabs(f_xi) < tol) exit

    dv2dxi = -2.0d0*(Smag2*aa**3 + SB**2*(aa**2 + aa*xi + xi**2)) / (xi*aa)**3

    ! D/(2·omv) = D·W/2; write as D·sqrt(W²)/2 to avoid dividing by small omv
    dpdxi = (gam_arg - 1.0d0) / gam_arg * &
            (omv2 + dv2dxi*(D_con*0.5d0*dsqrt(W2_loc) - aa))

    df_xi = 1.0d0 + 0.5d0*B2*dv2dxi + SB**2/aa**3 - dpdxi

    xi = xi - f_xi / df_xi
    xi = max(xi, B2*(1.0d0 + 1.0d-12))

  enddo

  ! Recover primitive variables from converged xi
  aa   = xi - B2
  nW2  = (xi * aa)**2
  dW2  = max(nW2 - Smag2*aa**2 - SB**2*(aa + xi), nW2*1.0d-15)
  W2_loc = max(nW2/dW2, 1.0d0 + 1.0d-15)
  omv2 = 1.0d0 / W2_loc
  omv  = 1.0d0 / dsqrt(W2_loc)
  vB   = SB / aa   ! v.B = (S.B)/(xi - B^2)

  Q_c(IDN) = D_con * omv                                     ! rho = D/W
  Q_c(IVX) = (Sx + vB*Bx) / (xi * omv)                     ! ux = W*vx
  Q_c(IVY) = (Sy + vB*By) / (xi * omv)                     ! uy = W*vy
  Q_c(IVZ) = (Sz + vB*Bz) / (xi * omv)                     ! uz = W*vz
  Q_c(IPR) = (gam_arg - 1.0d0) / gam_arg * (aa*omv2 - D_con*omv)    ! p
  Q_c(ISC) = U_c(ISC) / D_con                               ! passive scalar

return
end subroutine consv2prim_cell

!-------------------------------------------------------------------
!       Cell Surface B field ===> Cell Center B field
!-------------------------------------------------------------------
subroutine CellCenterMagneticField(ibeg, ifin, jbeg, jfin, Bs, Bc)
use params, only : nxtot, nytot
implicit none
integer, intent(in)  :: ibeg, ifin, jbeg, jfin
real(8), intent(in)  :: Bs(3,nxtot,nytot)
real(8), intent(out) :: Bc(3,nxtot,nytot)
integer :: i, j

!$omp do collapse(2) schedule(static) private(i,j)
  do j = jbeg, jfin
  do i = ibeg, ifin
    Bc(1,i,j) = 0.5d0*( Bs(1,i+1,j) + Bs(1,i,j) )
    Bc(2,i,j) = 0.5d0*( Bs(2,i,j+1) + Bs(2,i,j) )
    Bc(3,i,j) = Bs(3,i,j)
  enddo
  enddo
!$omp end do

return
end subroutine CellCenterMagneticField

!-------------------------------------------------------------------
!       Primitive -> Conservative wrapper (2D, NO OMP, called outside parallel)
!-------------------------------------------------------------------
subroutine Prim2Consv(Q, Bc, U)
use params, only : NVAR, nxtot, nytot, is, ie, js, je, gam
implicit none
real(8), intent(in)  :: Q(NVAR, nxtot, nytot)
real(8), intent(in)  :: Bc(3, nxtot, nytot)
real(8), intent(out) :: U(NVAR, nxtot, nytot)
integer :: i, j
real(8) :: Q_c(NVAR), Bc_c(3), U_c(NVAR)

do j = js, je
do i = is, ie
  Q_c(1:NVAR)  = Q(1:NVAR, i, j)
  Bc_c(1:3)    = Bc(1:3, i, j)
  call prim2consv_cell(Q_c, Bc_c, U_c, gam)
  U(1:NVAR, i, j) = U_c(1:NVAR)
enddo
enddo

return
end subroutine Prim2Consv

!-------------------------------------------------------------------
!       Conservative -> Primitive wrapper (2D, HAS !$omp do, called inside parallel)
!-------------------------------------------------------------------
subroutine Consv2Prim(U, Bs, Q, Bc)
use params, only : NVAR, nxtot, nytot, is, ie, js, je, gam
implicit none
real(8), intent(in)  :: U(NVAR, nxtot, nytot), Bs(3, nxtot, nytot)
real(8), intent(out) :: Q(NVAR, nxtot, nytot), Bc(3, nxtot, nytot)
integer :: i, j
real(8) :: U_c(NVAR), Bc_c(3), Q_c(NVAR)

call CellCenterMagneticField(is, ie, js, je, Bs, Bc)

!$omp do collapse(2) schedule(static) private(i, j, U_c, Bc_c, Q_c)
do j = js, je
do i = is, ie
  U_c(1:NVAR)  = U(1:NVAR, i, j)
  Bc_c(1:3)    = Bc(1:3, i, j)
  call consv2prim_cell(U_c, Bc_c, Q_c, gam)
  Q(1:NVAR, i, j) = Q_c(1:NVAR)
enddo
enddo
!$omp end do

return
end subroutine Consv2Prim

end module prim_mod
