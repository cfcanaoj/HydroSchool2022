module riemann_mod
implicit none
contains

!---------------------------------------------------------------------
!     van Leer monotonicity limiter
!---------------------------------------------------------------------
subroutine vanLeer(n, dvp, dvm, dv)
use params, only : NVAR
implicit none
integer,  intent(in)  :: n
real(8),  intent(in)  :: dvp(NVAR), dvm(NVAR)
real(8),  intent(out) :: dv(NVAR)
integer :: i

  do i = 1, n
    if (dvp(i)*dvm(i) > 0.0d0) then
      dv(i) = 2.0d0*dvp(i)*dvm(i) / (dvp(i) + dvm(i))
    else
      dv(i) = 0.0d0
    endif
  enddo

return
end subroutine vanLeer

!---------------------------------------------------------------------
!     SRMHD HLL Riemann Solver
!---------------------------------------------------------------------
subroutine HLL(idir, Ql, Qr, b1, dt, dx, flx, wghtCT, gam_arg)
use params, only : IDN, IVX, IVY, IVZ, IPR, ISC, IBX, IBY, IBZ, IEN, NFLX
implicit none
integer, intent(in)  :: idir
real(8), intent(in)  :: Ql(NFLX), Qr(NFLX)
real(8), intent(in)  :: b1, dt, dx
real(8), intent(out) :: flx(NFLX), wghtCT
real(8), intent(in)  :: gam_arg
integer :: IVpara, IVperp1, IVperp2, IBpara, IBperp1, IBperp2
real(8) :: Ul(NFLX), Ur(NFLX), Fl(NFLX), Fr(NFLX)
real(8) :: sl, sr, v_over_c
! Left state primitives and derived quantities
real(8) :: rho_l, vpara_l, vperp1_l, vperp2_l, p_l, Bp1_l, Bp2_l
real(8) :: v2_l, W2_l, W_l, h_l, xi_l, vB_l, Btot2_l, b2_l, ptot_l, Spara_l
real(8) :: rho_h_l, cs2_l, vA2_l, cf_l, lam_lm, lam_lp
! Right state primitives and derived quantities
real(8) :: rho_r, vpara_r, vperp1_r, vperp2_r, p_r, Bp1_r, Bp2_r
real(8) :: v2_r, W2_r, W_r, h_r, xi_r, vB_r, Btot2_r, b2_r, ptot_r, Spara_r
real(8) :: rho_h_r, cs2_r, vA2_r, cf_r, lam_rm, lam_rp

if (idir == 1) then
  IVpara = IVX;  IVperp1 = IVY;  IVperp2 = IVZ
  IBpara = IBX;  IBperp1 = IBY;  IBperp2 = IBZ
else
  IVpara = IVY;  IVperp1 = IVZ;  IVperp2 = IVX
  IBpara = IBY;  IBperp1 = IBZ;  IBperp2 = IBX
endif

!--- Left state ---
rho_l    = Ql(IDN);    vpara_l  = Ql(IVpara)
vperp1_l = Ql(IVperp1); vperp2_l = Ql(IVperp2); p_l = Ql(IPR)
Bp1_l    = Ql(IBperp1); Bp2_l    = Ql(IBperp2)

! Ql stores 4-velocity u^i; convert to v^i
v2_l  = vpara_l**2 + vperp1_l**2 + vperp2_l**2   ! u^2
W2_l  = 1.0d0 + v2_l                               ! W^2 = 1 + u^2
W_l   = dsqrt(W2_l)
vpara_l  = vpara_l  / W_l
vperp1_l = vperp1_l / W_l
vperp2_l = vperp2_l / W_l
v2_l  = v2_l / W2_l    ! now v^2
h_l   = 1.0d0 + gam_arg / (gam_arg - 1.0d0) * p_l / rho_l
Btot2_l = b1**2 + Bp1_l**2 + Bp2_l**2
vB_l    = vpara_l*b1 + vperp1_l*Bp1_l + vperp2_l*Bp2_l
b2_l    = Btot2_l / W2_l + vB_l**2
ptot_l  = p_l + 0.5d0*b2_l
xi_l    = rho_l*h_l*W2_l + Btot2_l
Spara_l = xi_l*vpara_l - vB_l*b1

Ul(IDN)    = rho_l * W_l
Ul(IVpara) = Spara_l
Ul(IVperp1)= xi_l*vperp1_l - vB_l*Bp1_l
Ul(IVperp2)= xi_l*vperp2_l - vB_l*Bp2_l
Ul(IEN)    = rho_l*h_l*W2_l - p_l - rho_l*W_l &
           + 0.5d0*Btot2_l*(1.0d0 + v2_l) - 0.5d0*vB_l**2
Ul(IBperp1)= Bp1_l
Ul(IBperp2)= Bp2_l

Fl(IDN)    = Ul(IDN) * vpara_l
Fl(IVpara) = Spara_l*vpara_l - b1*(vB_l*vpara_l + b1/W2_l) + ptot_l
Fl(IVperp1)= Ul(IVperp1)*vpara_l - b1*(vB_l*vperp1_l + Bp1_l/W2_l)
Fl(IVperp2)= Ul(IVperp2)*vpara_l - b1*(vB_l*vperp2_l + Bp2_l/W2_l)
Fl(IEN)    = Spara_l - Ul(IDN)*vpara_l  ! = (tau + ptot)*vpara - vB*b1
Fl(IBperp1)= Bp1_l*vpara_l - b1*vperp1_l
Fl(IBperp2)= Bp2_l*vpara_l - b1*vperp2_l

!--- Right state ---
rho_r    = Qr(IDN);    vpara_r  = Qr(IVpara)
vperp1_r = Qr(IVperp1); vperp2_r = Qr(IVperp2); p_r = Qr(IPR)
Bp1_r    = Qr(IBperp1); Bp2_r    = Qr(IBperp2)

! Qr stores 4-velocity u^i; convert to v^i
v2_r  = vpara_r**2 + vperp1_r**2 + vperp2_r**2   ! u^2
W2_r  = 1.0d0 + v2_r                               ! W^2 = 1 + u^2
W_r   = dsqrt(W2_r)
vpara_r  = vpara_r  / W_r
vperp1_r = vperp1_r / W_r
vperp2_r = vperp2_r / W_r
v2_r  = v2_r / W2_r    ! now v^2
h_r   = 1.0d0 + gam_arg / (gam_arg - 1.0d0) * p_r / rho_r
Btot2_r = b1**2 + Bp1_r**2 + Bp2_r**2
vB_r    = vpara_r*b1 + vperp1_r*Bp1_r + vperp2_r*Bp2_r
b2_r    = Btot2_r / W2_r + vB_r**2
ptot_r  = p_r + 0.5d0*b2_r
xi_r    = rho_r*h_r*W2_r + Btot2_r
Spara_r = xi_r*vpara_r - vB_r*b1

Ur(IDN)    = rho_r * W_r
Ur(IVpara) = Spara_r
Ur(IVperp1)= xi_r*vperp1_r - vB_r*Bp1_r
Ur(IVperp2)= xi_r*vperp2_r - vB_r*Bp2_r
Ur(IEN)    = rho_r*h_r*W2_r - p_r - rho_r*W_r &
           + 0.5d0*Btot2_r*(1.0d0 + v2_r) - 0.5d0*vB_r**2
Ur(IBperp1)= Bp1_r
Ur(IBperp2)= Bp2_r

Fr(IDN)    = Ur(IDN) * vpara_r
Fr(IVpara) = Spara_r*vpara_r - b1*(vB_r*vpara_r + b1/W2_r) + ptot_r
Fr(IVperp1)= Ur(IVperp1)*vpara_r - b1*(vB_r*vperp1_r + Bp1_r/W2_r)
Fr(IVperp2)= Ur(IVperp2)*vpara_r - b1*(vB_r*vperp2_r + Bp2_r/W2_r)
Fr(IEN)    = Spara_r - Ur(IDN)*vpara_r  ! = (tau + ptot)*vpara - vB*b1
Fr(IBperp1)= Bp1_r*vpara_r - b1*vperp1_r
Fr(IBperp2)= Bp2_r*vpara_r - b1*vperp2_r

!--- Wave speed estimates: relativistic velocity addition ---
rho_h_l = rho_l * h_l
cs2_l   = gam_arg * p_l / rho_h_l
vA2_l   = b2_l / (rho_h_l + b2_l)
cf_l    = dsqrt(cs2_l + vA2_l*(1.0d0 - cs2_l))
lam_lm  = (vpara_l - cf_l) / (1.0d0 - vpara_l*cf_l)
lam_lp  = (vpara_l + cf_l) / (1.0d0 + vpara_l*cf_l)

rho_h_r = rho_r * h_r
cs2_r   = gam_arg * p_r / rho_h_r
vA2_r   = b2_r / (rho_h_r + b2_r)
cf_r    = dsqrt(cs2_r + vA2_r*(1.0d0 - cs2_r))
lam_rm  = (vpara_r - cf_r) / (1.0d0 - vpara_r*cf_r)
lam_rp  = (vpara_r + cf_r) / (1.0d0 + vpara_r*cf_r)

sl = min(lam_lm, lam_rm)
sr = max(lam_lp, lam_rp)

!--- HLL flux ---
if (sl > 0.0d0) then
  flx(:) = Fl(:)
else if (sr <= 0.0d0) then
  flx(:) = Fr(:)
else
  flx(:) = (sr*Fl(:) - sl*Fr(:) + sl*sr*(Ur(:) - Ul(:))) / (sr - sl)
endif

flx(IBpara) = 0.0d0   ! parallel B handled by constrained transport

! CT weight: upwind estimate for electric field interpolation
v_over_c = 1024.0d0*dt*flx(IDN) / (dx*(Ul(IDN) + Ur(IDN)))
wghtCT   = 0.5d0 + max(-0.5d0, min(0.5d0, v_over_c))

! Passive scalar: upwind advection from mass flux direction
if (flx(IDN) >= 0.0d0) then
  flx(ISC) = flx(IDN) * Ql(ISC)
else
  flx(ISC) = flx(IDN) * Qr(ISC)
endif

return
end subroutine HLL

!---------------------------------------------------------------------
!     SRMHD HLLC Riemann Solver (Mignone & Bodo 2006, HLLC-MB)
!---------------------------------------------------------------------
subroutine HLLC(idir, Ql, Qr, b1, dt, dx, flx, wghtCT, gam_arg)
use params, only : IDN, IVX, IVY, IVZ, IPR, ISC, IBX, IBY, IBZ, IEN, NFLX
implicit none
integer, intent(in)  :: idir
real(8), intent(in)  :: Ql(NFLX), Qr(NFLX)
real(8), intent(in)  :: b1, dt, dx
real(8), intent(out) :: flx(NFLX), wghtCT
real(8), intent(in)  :: gam_arg
integer :: IVpara, IVperp1, IVperp2, IBpara, IBperp1, IBperp2
real(8) :: Ul(NFLX), Ur(NFLX), Fl(NFLX), Fr(NFLX)
real(8) :: sl, sr, v_over_c
! Left state
real(8) :: rho_l, vpara_l, vperp1_l, vperp2_l, p_l, Bp1_l, Bp2_l
real(8) :: v2_l, W2_l, W_l, h_l, xi_l, vB_l, Btot2_l, b2_l, ptot_l, Spara_l
real(8) :: rho_h_l, cs2_l, vA2_l, cf_l, lam_lm, lam_lp
! Right state
real(8) :: rho_r, vpara_r, vperp1_r, vperp2_r, p_r, Bp1_r, Bp2_r
real(8) :: v2_r, W2_r, W_r, h_r, xi_r, vB_r, Btot2_r, b2_r, ptot_r, Spara_r
real(8) :: rho_h_r, cs2_r, vA2_r, cf_r, lam_rm, lam_rp
! HLL single state and flux
real(8) :: D_hll, mx_hll, my_hll, mz_hll, E_hll, Bp1_hll, Bp2_hll
real(8) :: FD_hll, Fmx_hll, Fmy_hll, Fmz_hll, FE_hll, FBp1_hll, FBp2_hll
real(8) :: Eprime_hll
! Star region
real(8) :: aa_q, bb_q, cc_q, discr
real(8) :: vxstar, vystar, vzstar, vBstar, vstar2, omv2star
real(8) :: Bp1star, Bp2star, pstar
real(8) :: Dl_star, El_star, mxl_star, myl_star, mzl_star
real(8) :: Dr_star, Er_star, mxr_star, myr_star, mzr_star
real(8) :: Ul_star(NFLX), Ur_star(NFLX)

if (idir == 1) then
  IVpara = IVX;  IVperp1 = IVY;  IVperp2 = IVZ
  IBpara = IBX;  IBperp1 = IBY;  IBperp2 = IBZ
else
  IVpara = IVY;  IVperp1 = IVZ;  IVperp2 = IVX
  IBpara = IBY;  IBperp1 = IBZ;  IBperp2 = IBX
endif

!--- Left state: primitives -> U, F ---
rho_l    = Ql(IDN);     vpara_l  = Ql(IVpara)
vperp1_l = Ql(IVperp1); vperp2_l = Ql(IVperp2); p_l = Ql(IPR)
Bp1_l    = Ql(IBperp1); Bp2_l    = Ql(IBperp2)

! Ql stores 4-velocity u^i; convert to v^i
v2_l    = vpara_l**2 + vperp1_l**2 + vperp2_l**2   ! u^2
W2_l    = 1.0d0 + v2_l;  W_l = dsqrt(W2_l)         ! W^2 = 1 + u^2
vpara_l  = vpara_l  / W_l
vperp1_l = vperp1_l / W_l
vperp2_l = vperp2_l / W_l
v2_l    = v2_l / W2_l    ! now v^2
h_l     = 1.0d0 + gam_arg/(gam_arg-1.0d0)*p_l/rho_l
Btot2_l = b1**2 + Bp1_l**2 + Bp2_l**2
vB_l    = vpara_l*b1 + vperp1_l*Bp1_l + vperp2_l*Bp2_l
b2_l    = Btot2_l/W2_l + vB_l**2;  ptot_l = p_l + 0.5d0*b2_l
xi_l    = rho_l*h_l*W2_l + Btot2_l;  Spara_l = xi_l*vpara_l - vB_l*b1

Ul(IDN)     = rho_l*W_l;    Ul(IVpara)  = Spara_l
Ul(IVperp1) = xi_l*vperp1_l - vB_l*Bp1_l
Ul(IVperp2) = xi_l*vperp2_l - vB_l*Bp2_l
Ul(IEN)     = rho_l*h_l*W2_l - p_l - rho_l*W_l &
            + 0.5d0*Btot2_l*(1.0d0+v2_l) - 0.5d0*vB_l**2
Ul(IBperp1) = Bp1_l;         Ul(IBperp2) = Bp2_l

Fl(IDN)     = Ul(IDN)*vpara_l
Fl(IVpara)  = Spara_l*vpara_l - b1*(vB_l*vpara_l + b1/W2_l) + ptot_l
Fl(IVperp1) = Ul(IVperp1)*vpara_l - b1*(vB_l*vperp1_l + Bp1_l/W2_l)
Fl(IVperp2) = Ul(IVperp2)*vpara_l - b1*(vB_l*vperp2_l + Bp2_l/W2_l)
Fl(IEN)     = Spara_l - Ul(IDN)*vpara_l  ! F_tau = S^x - D*vx
Fl(IBperp1) = Bp1_l*vpara_l - b1*vperp1_l
Fl(IBperp2) = Bp2_l*vpara_l - b1*vperp2_l

!--- Right state: primitives -> U, F ---
rho_r    = Qr(IDN);     vpara_r  = Qr(IVpara)
vperp1_r = Qr(IVperp1); vperp2_r = Qr(IVperp2); p_r = Qr(IPR)
Bp1_r    = Qr(IBperp1); Bp2_r    = Qr(IBperp2)

! Qr stores 4-velocity u^i; convert to v^i
v2_r    = vpara_r**2 + vperp1_r**2 + vperp2_r**2   ! u^2
W2_r    = 1.0d0 + v2_r;  W_r = dsqrt(W2_r)         ! W^2 = 1 + u^2
vpara_r  = vpara_r  / W_r
vperp1_r = vperp1_r / W_r
vperp2_r = vperp2_r / W_r
v2_r    = v2_r / W2_r    ! now v^2
h_r     = 1.0d0 + gam_arg/(gam_arg-1.0d0)*p_r/rho_r
Btot2_r = b1**2 + Bp1_r**2 + Bp2_r**2
vB_r    = vpara_r*b1 + vperp1_r*Bp1_r + vperp2_r*Bp2_r
b2_r    = Btot2_r/W2_r + vB_r**2;  ptot_r = p_r + 0.5d0*b2_r
xi_r    = rho_r*h_r*W2_r + Btot2_r;  Spara_r = xi_r*vpara_r - vB_r*b1

Ur(IDN)     = rho_r*W_r;    Ur(IVpara)  = Spara_r
Ur(IVperp1) = xi_r*vperp1_r - vB_r*Bp1_r
Ur(IVperp2) = xi_r*vperp2_r - vB_r*Bp2_r
Ur(IEN)     = rho_r*h_r*W2_r - p_r - rho_r*W_r &
            + 0.5d0*Btot2_r*(1.0d0+v2_r) - 0.5d0*vB_r**2
Ur(IBperp1) = Bp1_r;         Ur(IBperp2) = Bp2_r

Fr(IDN)     = Ur(IDN)*vpara_r
Fr(IVpara)  = Spara_r*vpara_r - b1*(vB_r*vpara_r + b1/W2_r) + ptot_r
Fr(IVperp1) = Ur(IVperp1)*vpara_r - b1*(vB_r*vperp1_r + Bp1_r/W2_r)
Fr(IVperp2) = Ur(IVperp2)*vpara_r - b1*(vB_r*vperp2_r + Bp2_r/W2_r)
Fr(IEN)     = Spara_r - Ur(IDN)*vpara_r  ! F_tau = S^x - D*vx
Fr(IBperp1) = Bp1_r*vpara_r - b1*vperp1_r
Fr(IBperp2) = Bp2_r*vpara_r - b1*vperp2_r

!--- Fast wave speeds (same as HLL) ---
rho_h_l = rho_l*h_l
cs2_l   = gam_arg*p_l/rho_h_l;  vA2_l = b2_l/(rho_h_l+b2_l)
cf_l    = dsqrt(cs2_l + vA2_l*(1.0d0-cs2_l))
lam_lm  = (vpara_l-cf_l)/(1.0d0-vpara_l*cf_l)
lam_lp  = (vpara_l+cf_l)/(1.0d0+vpara_l*cf_l)

rho_h_r = rho_r*h_r
cs2_r   = gam_arg*p_r/rho_h_r;  vA2_r = b2_r/(rho_h_r+b2_r)
cf_r    = dsqrt(cs2_r + vA2_r*(1.0d0-cs2_r))
lam_rm  = (vpara_r-cf_r)/(1.0d0-vpara_r*cf_r)
lam_rp  = (vpara_r+cf_r)/(1.0d0+vpara_r*cf_r)

sl = min(lam_lm, lam_rm)
sr = max(lam_lp, lam_rp)

!--- HLL single state (Eq. 15) and flux (Eq. 16) ---
D_hll    = (sr*Ur(IDN)    -sl*Ul(IDN)    +Fl(IDN)    -Fr(IDN)   )/(sr-sl)
mx_hll   = (sr*Ur(IVpara) -sl*Ul(IVpara) +Fl(IVpara) -Fr(IVpara))/(sr-sl)
my_hll   = (sr*Ur(IVperp1)-sl*Ul(IVperp1)+Fl(IVperp1)-Fr(IVperp1))/(sr-sl)
mz_hll   = (sr*Ur(IVperp2)-sl*Ul(IVperp2)+Fl(IVperp2)-Fr(IVperp2))/(sr-sl)
E_hll    = (sr*Ur(IEN)    -sl*Ul(IEN)    +Fl(IEN)    -Fr(IEN)   )/(sr-sl)
Bp1_hll  = (sr*Ur(IBperp1)-sl*Ul(IBperp1)+Fl(IBperp1)-Fr(IBperp1))/(sr-sl)
Bp2_hll  = (sr*Ur(IBperp2)-sl*Ul(IBperp2)+Fl(IBperp2)-Fr(IBperp2))/(sr-sl)

FD_hll   = (sr*Fl(IDN)    -sl*Fr(IDN)    +sl*sr*(Ur(IDN)   -Ul(IDN)   ))/(sr-sl)
Fmx_hll  = (sr*Fl(IVpara) -sl*Fr(IVpara) +sl*sr*(Ur(IVpara)-Ul(IVpara)))/(sr-sl)
Fmy_hll  = (sr*Fl(IVperp1)-sl*Fr(IVperp1)+sl*sr*(Ur(IVperp1)-Ul(IVperp1)))/(sr-sl)
Fmz_hll  = (sr*Fl(IVperp2)-sl*Fr(IVperp2)+sl*sr*(Ur(IVperp2)-Ul(IVperp2)))/(sr-sl)
FE_hll   = (sr*Fl(IEN)    -sl*Fr(IEN)    +sl*sr*(Ur(IEN)   -Ul(IEN)   ))/(sr-sl)
FBp1_hll = (sr*Fl(IBperp1)-sl*Fr(IBperp1)+sl*sr*(Ur(IBperp1)-Ul(IBperp1)))/(sr-sl)
FBp2_hll = (sr*Fl(IBperp2)-sl*Fr(IBperp2)+sl*sr*(Ur(IBperp2)-Ul(IBperp2)))/(sr-sl)

Eprime_hll = E_hll + D_hll   ! E' = E + D  (energy including rest mass)

! suppress unused warnings
if (.false.) then
  my_hll = my_hll; mz_hll = mz_hll; Fmy_hll = Fmy_hll; Fmz_hll = Fmz_hll
endif

!--- HLLC-MB: require Bx large enough relative to total B (else HLL) ---
! When b1 is tiny compared to transverse B, Step 3 division by b1 is ill-conditioned
if (b1**2 < 1.0d-4 * (b1**2 + Bp1_l**2 + Bp2_l**2 + Bp1_r**2 + Bp2_r**2 + 1.0d-30)) goto 100

!--- Step 1: B* = B^hll (Eq. 26) ---
Bp1star = Bp1_hll
Bp2star = Bp2_hll

!--- Step 2: v*_x from quadratic a*vx^2+b*vx+c=0 (Eq. 27-28) ---
aa_q = FE_hll + FD_hll - Bp1_hll*FBp1_hll - Bp2_hll*FBp2_hll
bb_q = -(Fmx_hll + Eprime_hll) + Bp1_hll**2 + Bp2_hll**2 &
       + FBp1_hll**2 + FBp2_hll**2
cc_q = mx_hll - Bp1_hll*FBp1_hll - Bp2_hll*FBp2_hll

discr = bb_q**2 - 4.0d0*aa_q*cc_q
if (discr < 0.0d0 .or. dabs(aa_q) < 1.0d-14) goto 100

vxstar = (-bb_q - dsqrt(discr)) / (2.0d0*aa_q)   ! negative branch (Mignone & Bodo 2006)
if (dabs(vxstar) >= 1.0d0) goto 100               ! causality check

!--- Step 3: v*_t from Eq. (29): b1*vt* = Bt**v*_x - F^hll_Bt ---
vystar = (Bp1star*vxstar - FBp1_hll) / b1
vzstar = (Bp2star*vxstar - FBp2_hll) / b1

vstar2 = vxstar**2 + vystar**2 + vzstar**2
if (vstar2 >= 1.0d0) goto 100                      ! superluminal: fallback

omv2star = 1.0d0 - vstar2     ! 1/W*^2 = 1 - v*^2
vBstar   = b1*vxstar + Bp1star*vystar + Bp2star*vzstar

!--- Step 4: p* from Eq. (30) ---
pstar = Fmx_hll - (FE_hll + FD_hll - b1*vBstar)*vxstar + b1*b1*omv2star
if (pstar < 0.0d0) goto 100                        ! unphysical pressure

!--- Step 5: D*, E*, m*_t, m*_x for each side (Eqs. 31, 24) ---
! Left (K=L, lambda_K=sl)
Dl_star  = Ul(IDN) * (sl - vpara_l) / (sl - vxstar)
El_star  = (sl*Ul(IEN) - Fl(IEN) + pstar*vxstar - vBstar*b1) / (sl - vxstar)
myl_star = (sl*Ul(IVperp1) - Fl(IVperp1) &
             - b1*(Bp1star*omv2star + vBstar*vystar)) / (sl - vxstar)
mzl_star = (sl*Ul(IVperp2) - Fl(IVperp2) &
             - b1*(Bp2star*omv2star + vBstar*vzstar)) / (sl - vxstar)
mxl_star = (El_star + Dl_star + pstar)*vxstar - vBstar*b1   ! Eq. 24

! Right (K=R, lambda_K=sr)
Dr_star  = Ur(IDN) * (sr - vpara_r) / (sr - vxstar)
Er_star  = (sr*Ur(IEN) - Fr(IEN) + pstar*vxstar - vBstar*b1) / (sr - vxstar)
myr_star = (sr*Ur(IVperp1) - Fr(IVperp1) &
             - b1*(Bp1star*omv2star + vBstar*vystar)) / (sr - vxstar)
mzr_star = (sr*Ur(IVperp2) - Fr(IVperp2) &
             - b1*(Bp2star*omv2star + vBstar*vzstar)) / (sr - vxstar)
mxr_star = (Er_star + Dr_star + pstar)*vxstar - vBstar*b1   ! Eq. 24

!--- Build star state arrays ---
Ul_star(IDN)     = Dl_star;   Ul_star(IVpara)  = mxl_star
Ul_star(IVperp1) = myl_star;  Ul_star(IVperp2) = mzl_star
Ul_star(IEN)     = El_star
Ul_star(IBperp1) = Bp1star;   Ul_star(IBperp2) = Bp2star
Ul_star(IBpara)  = 0.0d0;     Ul_star(ISC)     = 0.0d0

Ur_star(IDN)     = Dr_star;   Ur_star(IVpara)  = mxr_star
Ur_star(IVperp1) = myr_star;  Ur_star(IVperp2) = mzr_star
Ur_star(IEN)     = Er_star
Ur_star(IBperp1) = Bp1star;   Ur_star(IBperp2) = Bp2star
Ur_star(IBpara)  = 0.0d0;     Ur_star(ISC)     = 0.0d0

!--- Flux selection (Eq. 19, lambda* = v*_x) ---
if (sl > 0.0d0) then
  flx(:) = Fl(:)
else if (vxstar >= 0.0d0) then
  flx(:) = Fl(:) + sl*(Ul_star(:) - Ul(:))
else if (sr >= 0.0d0) then
  flx(:) = Fr(:) + sr*(Ur_star(:) - Ur(:))
else
  flx(:) = Fr(:)
endif

flx(IBpara) = 0.0d0   ! parallel B handled by constrained transport

v_over_c = 1024.0d0*dt*flx(IDN) / (dx*(Ul(IDN)+Ur(IDN)))
wghtCT   = 0.5d0 + max(-0.5d0, min(0.5d0, v_over_c))

if (flx(IDN) >= 0.0d0) then
  flx(ISC) = flx(IDN) * Ql(ISC)
else
  flx(ISC) = flx(IDN) * Qr(ISC)
endif

return

! Fallback to HLL (Bx=0, non-physical star state, or superluminal velocity)
100 continue
call HLL(idir, Ql, Qr, b1, dt, dx, flx, wghtCT, gam_arg)

return
end subroutine HLLC

!---------------------------------------------------------------------
!     SRMHD HLLD Riemann Solver (placeholder - delegates to HLL)
!---------------------------------------------------------------------
subroutine HLLD(idir, Ql, Qr, b1, dt, dx, flx, wghtCT, gam_arg)
use params, only : NFLX
implicit none
integer, intent(in)  :: idir
real(8), intent(in)  :: Ql(NFLX), Qr(NFLX)
real(8), intent(in)  :: b1, dt, dx
real(8), intent(out) :: flx(NFLX), wghtCT
real(8), intent(in)  :: gam_arg

! SRMHD HLLD (Mignone et al. 2009) is not yet implemented.
! Delegate to SRMHD HLL.
call HLL(idir, Ql, Qr, b1, dt, dx, flx, wghtCT, gam_arg)

return
end subroutine HLLD

end module riemann_mod
