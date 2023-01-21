module commons
implicit none
integer::ntime                        ! counter of the timestep
integer,parameter::ntimemax=200000     ! the maximum timesteps
real(8)::time,dt                    ! time, timewidth
data time / 0.0d0 /
real(8),parameter:: timemax=0.1d0   
real(8),parameter:: dtout=5.0d-3

integer,parameter::ngrid=256       ! the number of grids in the simulation box
integer,parameter::mgn=2            ! the number of ghost cells
integer,parameter::in=ngrid+2*mgn+1 ! the total number of grids including ghost cells
integer,parameter::is=mgn+1         ! the index of the leftmost grid
integer,parameter::ie=ngrid+mgn     ! the index of the rightmost grid
integer,parameter::nflux = 3
real(8),parameter:: x1min=-0.5d0,x1max=0.5d0
real(8),parameter:: Bx=0.75

integer, parameter :: IDN = 1
integer, parameter :: IM1 = 2
integer, parameter :: IM2 = 3
integer, parameter :: IM3 = 4
integer, parameter :: IPR = 5
integer, parameter :: IB2 = 6
integer, parameter :: IB3 = 7
integer, parameter :: NVAR = 7
integer, parameter :: NFLX = 7

integer, parameter :: IV1 = 2
integer, parameter :: IV2 = 3
integer, parameter :: IV3 = 4
integer, parameter :: IEN = 5

end module commons
     
      module eosmod
      implicit none
! adiabatic
      real(8),parameter::gam=5.0d0/3.0d0 !! adiabatic index
!      real(8),parameter::gam=1.4d0 !! adiabatic index
! isothermal
!      real(8)::csiso  !! isothemal sound speed
      end module eosmod

      program main
      use commons
      implicit none

      real(8),dimension(in)::x1a,x1b
      real(8),dimension(in,NVAR) :: U
      real(8),dimension(in,NVAR) :: Q
      real(8),dimension(in,NFLX) :: flux

!      write(6,*) "setup grids and initial condition"
      call GenerateGrid(x1a, x1b)
      call GenerateProblem(x1b, Q)
      call ConsvVariable(Q, U)
      call Output( x1a, x1b, Q, .FALSE. )

! main loop
      mloop: do ntime=1,ntimemax
         call TimestepControl(x1a, Q)
         if( time + dt > timemax ) dt = timemax - time
         call BoundaryCondition( Q)
         call NumericalFlux( Q, flux )
         call UpdateConsv( flux, U)
         call PrimVariable( U, Q )
         time=time+dt
         call Output( x1a, x1b, Q, .TRUE.)

         if(time >= timemax) exit mloop
      enddo mloop
      call Output( x1a, x1b, Q, .FALSE.)

!      write(6,*) "program has been finished"
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

      subroutine GenerateProblem(x1b, Q)
      use commons
      use eosmod
      implicit none
      integer::i
      real(8), intent(in ) :: x1b(:)
      real(8), intent(out) :: Q(:,:)
      real(8) :: rho1,rho2,Lsm,u1,u2


      do i=is,ie
         if( x1b(i) < 0.0d0 ) then 
             Q(i,IDN) = 1.0d0
             Q(i,IV1) = 0.0d0
             Q(i,IV2) = 0.0d0
             Q(i,IV3) = 0.0d0
             Q(i,IPR) = 1.0d0
             Q(i,IB2) = 1.0d0
             Q(i,IB3) = 0.0d0
         else 
             Q(i,IDN) = 0.125d0
             Q(i,IV1) = 0.0d0
             Q(i,IV2) = 0.0d0
             Q(i,IV3) = 0.0d0
             Q(i,IPR) = 0.1d0
             Q(i,IB2) = -1.0d0
             Q(i,IB3) = 0.0d0
         endif
      enddo

      
!      call BoundaryCondition

      return
      end subroutine GenerateProblem

      subroutine BoundaryCondition(Q)
      use commons
      implicit none
      real(8), intent(inout) :: Q(:,:)
      integer::i


      do i=1,mgn
          Q(is-i,IDN)  = Q(is-1+i,IDN)
          Q(is-i,IV1)  = Q(is-1+i,IV1)
          Q(is-i,IV2)  = Q(is-1+i,IV2)
          Q(is-i,IV3)  = Q(is-1+i,IV3)
          Q(is-i,IPR)  = Q(is-1+i,IPR)
          Q(is-i,IB2)  = Q(is-1+i,IB2)
          Q(is-i,IB3)  = Q(is-1+i,IB3)
      enddo

      do i=1,mgn
          Q(ie+i,IDN) = Q(ie-i+1,IDN)
          Q(ie+i,IV1) = Q(ie-i+1,IV1)
          Q(ie+i,IV2) = Q(ie-i+1,IV2)
          Q(ie+i,IV3) = Q(ie-i+1,IV3)
          Q(ie+i,IPR) = Q(ie-i+1,IPR)
          Q(ie+i,IB2) = Q(ie-i+1,IB2)
          Q(ie+i,IB3) = Q(ie-i+1,IB3)
      enddo

      return
      end subroutine BoundaryCondition
!
      subroutine ConsvVariable(Q, U)
      use commons
      use eosmod
      implicit none
      real(8), intent(in) :: Q(:,:)
      real(8), intent(out) :: U(:,:)
      integer::i

      do i=is,ie
          U(i,IDN) = Q(i,IDN)
          U(i,IM1) = Q(i,IDN)*Q(i,IV1)
          U(i,IM2) = Q(i,IDN)*Q(i,IV2)
          U(i,IM3) = Q(i,IDN)*Q(i,IV3)
          U(i,IEN) = 0.5d0*Q(i,IDN)*( Q(i,IV1)**2 + Q(i,IV2)**2 + Q(i,IV3)**2 ) &
                   + 0.5d0*( Bx**2 + Q(i,IB2)**2 + Q(i,IB3)**2 ) &
                   + Q(i,IPR)/(gam - 1.0d0)
          U(i,IB2) = Q(i,IB2)
          U(i,IB3) = Q(i,IB3)
      enddo
      
      return
      end subroutine Consvvariable

      subroutine PrimVariable( U, Q )
      use commons
      use eosmod
      implicit none
      real(8), intent(in) :: U(:,:)
      real(8), intent(out) :: Q(:,:)
      integer::i
      real(8) :: inv_d;

      do i=is,ie
           Q(i,IDN) = U(i,IDN)
           inv_d = 1.0d0/U(i,IDN)
           Q(i,IV1) = U(i,IM1)*inv_d
           Q(i,IV2) = U(i,IM2)*inv_d
           Q(i,IV3) = U(i,IM3)*inv_d
           Q(i,IPR) = ( U(i,IEN) &
                    - 0.5d0*(U(i,IM1)**2 + U(i,IM2)**2 + U(i,IM3)**2)*inv_d  &
                    - 0.5d0*(Bx**2 + U(i,IB2)**2 + U(i,IB3)**2) )*(gam-1.0d0)
           Q(i,IB2) = U(i,IB2)
           Q(i,IB3) = U(i,IB3)
      enddo

      return
      end subroutine PrimVariable

      subroutine TimestepControl(x1a, Q)
      use commons
      use eosmod
      implicit none
      real(8), intent(in) :: x1a(:), Q(:,:)
      real(8)::dtl1
      real(8)::dtl2
      real(8)::dtl3
      real(8)::dtlocal
      real(8)::dtmin,cf
      integer::i

      dtmin=1.0d90

      do i=is,ie
         cf = dsqrt( (gam*Q(i,IPR) + Bx**2 + Q(i,IB2)**2 + Q(i,IB3)**2)/Q(i,IDN))
         dtlocal =(x1a(i+1)-x1a(i))/(abs(Q(i,IV1)) + cf)
         if(dtlocal .lt. dtmin) dtmin = dtlocal
      enddo

      dt = 0.05d0 * dtmin
!      write(6,*)"dt",dt
      return
      end subroutine TimestepControl

!---------------------------------------------------------------------
!     van Leer monotonicity limiter 
!---------------------------------------------------------------------
      subroutine vanLeer(n,dvp,dvm,dv)
      implicit none
      real(8),intent(in)::dvp(:),dvm(:)
      integer,intent(in) :: n
      real(8),intent(out)::dv(:)
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
      use commons !, only: is, ie, in
      implicit none
      integer::i
      real(8), intent(in) :: Q(:,:)
      real(8), intent(out) :: flux(:,:)
      real(8),dimension(in,NFLX):: Ql,Qr
      real(8),dimension(NFLX):: flx
      real(8) :: dQm(NFLX), dQp(NFLX), dQmon(NFLX)
      real(8) :: ddmon, dvmon, dpmon

      do i=is-1,ie+1
         dQp(1:NVAR) = Q(i+1,1:NVAR) - Q(i  ,1:NVAR)
         dQm(1:NVAR) = Q(i  ,1:NVAR) - Q(i-1,1:NVAR)

         call vanLeer(NFLX, dQp, dQm, dQmon)

         Ql(i+1,1:NVAR) = Q(i,1:NVAR) + 0.5d0*dQmon(1:NVAR)*1.0d0
         Qr(i  ,1:NVAR) = Q(i,1:NVAR) - 0.5d0*dQmon(1:NVAR)*1.0d0
      enddo

      do i=is,ie+1
!         call HLL(Ql(i,:),Qr(i,:),flx)
         call HLLD(Ql(i,:),Qr(i,:),flx)
         flux(i,:)  = flx(:)
      enddo


      return
      end subroutine Numericalflux

!---------------------------------------------------------------------
!     HLL Riemann Solver
!---------------------------------------------------------------------
!     solve the HLL Riemann solver 
!
!     Input: Ql, Qr: primitive variables containing the perpendicular B fields 
!                    at the left and right states
!            1D array (IDN, IV1, IV2, IV3, IPR, IB2, IB3)
!                                 |
!                                 |
!                           Ql    |    Qr
!                                 |
!                                -->
!                                flx
!
!     Input: b1    : magnetic field perpendicular to the initial discontinuity
!
!     Output: flx  : flux estimated at the initial discontinuity
!            index: (IDN, IV1, IV2, IV3, IPR, IB2, IB3)
!---------------------------------------------------------------------
      subroutine HLL(Ql,Qr,flx)
      use commons !, only : is, ie, NVAR
      use eosmod
      implicit none
      real(8),intent(in)::Ql(:), Qr(:)
      real(8),intent(out) :: flx(:)
      real(8):: Ul(NFLX), Ur(NFLX)
      real(8):: Fl(NFLX), Fr(NFLX)
      real(8):: Ust(NFLX)
      real(8):: Fst(NFLX)
      real(8):: cfl,cfr
      real(8):: sl, sr
      real(8):: pbl, pbr, ptotl, ptotr
      integer :: i, n

          pbl = 0.5d0*(Bx**2 + Ql(IB2)**2 + Ql(IB3)**2)
          pbr = 0.5d0*(Bx**2 + Qr(IB2)**2 + Qr(IB3)**2)
          ptotl = Ql(IPR) + pbl
          ptotr = Qr(IPR) + pbr

          ! conserved variables in the left and right states
          Ul(IDN) = Ql(IDN)
          Ul(IM1) = Ql(IDN)*Ql(IV1)
          Ul(IM2) = Ql(IDN)*Ql(IV2)
          Ul(IM3) = Ql(IDN)*Ql(IV3)
          Ul(IEN) = 0.5d0*Ql(IDN)*( Ql(IV1)**2 + Ql(IV2)**2 + Ql(IV3)**2) & 
                  + pbl + Ql(IPR)/(gam - 1.0d0)
          Ul(IB2) = Ql(IB2)
          Ul(IB3) = Ql(IB3)

          Ur(IDN) = Qr(IDN)
          Ur(IM1) = Qr(IDN)*Qr(IV1)
          Ur(IM2) = Qr(IDN)*Qr(IV2)
          Ur(IM3) = Qr(IDN)*Qr(IV3)
          Ur(IEN) = 0.5d0*Qr(IDN)*( Qr(IV1)**2 + Qr(IV2)**2 + Qr(IV3)**2) & 
                  + pbr + Qr(IPR)/(gam - 1.0d0)
          Ur(IB2) = Qr(IB2)
          Ur(IB3) = Qr(IB3)

          ! flux in the left and right states
          Fl(IDN) = Ul(IM1)
          Fl(IM1) = Ul(IM1)*Ql(IV1) + ptotl - Bx**2
          Fl(IM2) = Ul(IM2)*Ql(IV1) - Bx*Ql(IB2)
          Fl(IM3) = Ul(IM3)*Ql(IV1) - Bx*Ql(IB3)
          Fl(IEN) = ( Ul(IEN) + ptotl )*Ql(IV1) &
                  - Bx*( Bx*Ql(IV1) + Ql(IB2)*Ql(IV2) + Ql(IB3)*Ql(IV3) )
          Fl(IB2) = Ql(IB2)*Ql(IV1) - Bx*Ql(IV2)
          Fl(IB3) = Ql(IB3)*Ql(IV1) - Bx*Ql(IV3)

          Fr(IDN) = Ur(IM1)
          Fr(IM1) = Ur(IM1)*Qr(IV1) + ptotr - Bx**2
          Fr(IM2) = Ur(IM2)*Qr(IV1) - Bx*Qr(IB2)
          Fr(IM3) = Ur(IM3)*Qr(IV1) - Bx*Qr(IB3)
          Fr(IEN) = ( Ur(IEN) + ptotr )*Qr(IV1) &
                  - Bx*( Bx*Qr(IV1) + Qr(IB2)*Qr(IV2) + Qr(IB3)*Qr(IV3) )
          Fr(IB2) = Qr(IB2)*Qr(IV1) - Bx*Qr(IV2)
          Fr(IB3) = Qr(IB3)*Qr(IV1) - Bx*Qr(IV3)

          cfl = dsqrt( (gam*Ql(IPR) + Ql(IB2)**2 + Ql(IB3)**2 + Bx**2)/Ql(IDN))
          cfr = dsqrt( (gam*Qr(IPR) + Qr(IB2)**2 + Qr(IB3)**2 + Bx**2)/Qr(IDN))

!          sl = min(Ql(IV1),Qr(IV1)) - max(cfl,cfr)
!          sr = max(Ql(IV1),Qr(IV1)) + max(cfl,cfr)
          sl = min(Ql(IV1) - cfl,Qr(IV1) - cfr)
          sr = max(Ql(IV1) + cfl,Qr(IV1) + cfr)
!          Fst(:)  = (sr*Fl(:) - sl*Fr(:) + sl*sr*( Ur(:) - Ul(:) ))/(sr - sl)
!          Ust(:) = ( sr*Ur(:) - sl*Ul(:) - Fr(:) + Fl(:) )/(sr - sl)


          if( sl > 0.0d0 ) then
               flx(:)  = Fl(:) 
          else if (sr <= 0.0d0 ) then
               flx(:)  = Fr(:) 
          else 
               flx(:)  = (sr*Fl(:) - sl*Fr(:) + sl*sr*( Ur(:) - Ul(:) ))/(sr - sl)
          endif

      return
      end subroutine HLL

!---------------------------------------------------------------------
!     HLLD Riemann Solver
!---------------------------------------------------------------------
!     solve the HLL Riemann solver 
!
!     Input: Ql, Qr: primitive variables containing the perpendicular B fields 
!                    at the left and right states
!            1D array (IDN, IV1, IV2, IV3, IPR, IBperp1, IBperp2)
!                                 |
!                                 |
!                           Ql    |    Qr
!                                 |
!                                -->
!                                flx
!
!     Input: b1    : magnetic field perpendicular to the initial discontinuity
!
!     Output: flx  : flux estimated at the initial discontinuity
!            index: (IDN, IV1, IV2, IV3, IPR, IBperp1, IBperp2)
!---------------------------------------------------------------------
      subroutine HLLD(Ql,Qr,flx)
      use commons !, only : is, ie, NVAR
      use eosmod
      implicit none
      real(8),intent(in)  ::Ql(:), Qr(:)
      real(8),intent(out) :: flx(:)
      real(8):: b1
      real(8):: Ul(NFLX), Ur(NFLX)
      real(8):: Ulst(NFLX), Urst(NFLX)
      real(8):: Uldst(NFLX), Urdst(NFLX)
      real(8):: Fl(NFLX), Fr(NFLX)
      real(8):: test(NFLX)
      real(8):: cfl,cfr
      real(8):: S0, S1, S2, S3, S4
      real(8):: pbl, pbr, ptotl, ptotr
      real(8) :: sqrtdl, sqrtdr, v_dot_B_stl, v_dot_B_str
      real(8) :: Ulst_d_inv, Urst_d_inv, sum_sqrtd_inv, tmp
      real(8) :: ptot_stl, ptot_str,ptot_st, Cl, Cr, Cml, Cmr, Cml_inv, Cmr_inv, bxsgn
      integer :: i, n

          pbl = 0.5d0*(Bx**2 + Ql(IB2)**2 + Ql(IB3)**2)
          pbr = 0.5d0*(Bx**2 + Qr(IB2)**2 + Qr(IB3)**2)
          ptotl = Ql(IPR) + pbl
          ptotr = Qr(IPR) + pbr

          cfl = dsqrt( 0.5d0*( 2.0d0*pbl + gam*Ql(IPR) &
                     + dsqrt( (2.0d0*pbl - gam*Ql(IPR))**2 &
                     + 4.0d0*gam*Ql(IPR)*( Ql(IB2)**2 + Ql(IB3)**2 ) ) )/Ql(IDN) )
          cfr = dsqrt( 0.5d0*( 2.0d0*pbr + gam*Qr(IPR) &
                     + dsqrt( (2.0d0*pbr - gam*Qr(IPR))**2 &
                     + 4.0d0*gam*Qr(IPR)*( Qr(IB2)**2 + Qr(IB3)**2 ) ) )/Qr(IDN) )

          S0 = min( Ql(IV1) - cfl, Qr(IV1) - cfr)
          S4 = max( Ql(IV1) + cfl, Qr(IV1) + cfr)

          ! conserved variables in the left and right states
          Ul(IDN) = Ql(IDN)
          Ul(IV1) = Ql(IDN)*Ql(IV1)
          Ul(IV2) = Ql(IDN)*Ql(IV2)
          Ul(IV3) = Ql(IDN)*Ql(IV3)
          Ul(IEN) = 0.5d0*Ql(IDN)*( Ql(IV1)**2 + Ql(IV2)**2 + Ql(IV3)**2) & 
                  + pbl + Ql(IPR)/(gam - 1.0d0)
          Ul(IB2) = Ql(IB2)
          Ul(IB3) = Ql(IB3)

          Ur(IDN) = Qr(IDN)
          Ur(IV1) = Qr(IDN)*Qr(IV1)
          Ur(IV2) = Qr(IDN)*Qr(IV2)
          Ur(IV3) = Qr(IDN)*Qr(IV3)
          Ur(IEN) = 0.5d0*Qr(IDN)*( Qr(IV1)**2 + Qr(IV2)**2 + Qr(IV3)**2) & 
                  + pbr + Qr(IPR)/(gam - 1.0d0)
          Ur(IB2) = Qr(IB2)
          Ur(IB3) = Qr(IB3)

    !--- Step 3.  Compute L/R fluxes
          Fl(IDN) = Ul(IV1)
          Fl(IV1) = Ul(IV1)*Ql(IV1) + ptotl - Bx**2
          Fl(IV2) = Ul(IV2)*Ql(IV1) - Bx*Ql(IB2)
          Fl(IV3) = Ul(IV3)*Ql(IV1) - Bx*Ql(IB3)
          Fl(IEN) = ( Ul(IEN) + ptotl - Bx**2 )*Ql(IV1) &
                  - Bx*( Ql(IB2)*Ql(IV2) + Ql(IB3)*Ql(IV3) )
          Fl(IB2) = Ql(IB2)*Ql(IV1) - Bx*Ql(IV2)
          Fl(IB3) = Ql(IB3)*Ql(IV1) - Bx*Ql(IV3)

          Fr(IDN) = Ur(IV1)
          Fr(IV1) = Ur(IV1)*Qr(IV1) + ptotr - Bx**2
          Fr(IV2) = Ur(IV2)*Qr(IV1) - Bx*Qr(IB2)
          Fr(IV3) = Ur(IV3)*Qr(IV1) - Bx*Qr(IB3)
          Fr(IEN) = ( Ur(IEN) + ptotr - Bx**2 )*Qr(IV1) &
                  - Bx*( Qr(IB2)*Qr(IV2) + Qr(IB3)*Qr(IV3) )
          Fr(IB2) = Qr(IB2)*Qr(IV1) - Bx*Qr(IV2)
          Fr(IB3) = Qr(IB3)*Qr(IV1) - Bx*Qr(IV3)

    !--- Step 4.  Compute middle and Alfven wave speeds
          Cl = S0 - Ql(IV1)
          Cr = S4 - Qr(IV1)

          S2 = ( Cr*Ur(IV1) - Cl*Ul(IV1) + (ptotl - ptotr) ) &
                 /( Cr*Ur(IDN) - Cl*Ul(IDN) )

          Cml = S0 - S2
          Cmr = S4 - S2
          Cml_inv = 1.0d0/Cml
          Cmr_inv = 1.0d0/Cmr

          Ulst(IDN) = Ul(IDN)*Cl*Cml_inv
          Urst(IDN) = Ur(IDN)*Cr*Cmr_inv
          Ulst_d_inv = 1.0d0/Ulst(IDN)
          Urst_d_inv = 1.0d0/Urst(IDN)
          sqrtdl = dsqrt(Ulst(IDN))
          sqrtdr = dsqrt(Urst(IDN))

!          if( sqrtdr .ne. sqrtdr ) then
!              print*, "sqrtdr",sqrtdr, Cr, Cmr
!             print*,"S", S0,S2,S4
!              stop
!          endif

          S1 = S2 - dabs(Bx)/sqrtdl
          S3 = S2 + dabs(Bx)/sqrtdr

    !--- Step 5.  Compute intermediate states
         ptot_stl = ptotl + Ul(IDN)*Cl*(S2 - Ql(IV1))
         ptot_str = ptotr + Ur(IDN)*Cr*(S2 - Qr(IV1))

         ptot_st = 0.5d0*(ptot_stl + ptot_str)

         Ulst(IV1) = Ulst(IDN)*S2
         if( dabs( Ul(IDN)*Cl*Cml-Bx**2) < 1.0d-8*ptot_st ) then
             Ulst(IV2) = Ulst(IDN)*Ql(IV2)
             Ulst(IV3) = Ulst(IDN)*Ql(IV3)

             Ulst(IB2) = Ul(IB2)
             Ulst(IB3) = Ul(IB3)
         else 
             tmp = Bx*( Cl - Cml )/(Ul(IDN)*Cl*Cml - Bx**2)
             Ulst(IV2) = Ulst(IDN)*( Ql(IV2) - Ul(IB2)*tmp )
             Ulst(IV3) = Ulst(IDN)*( Ql(IV3) - Ul(IB3)*tmp )

             tmp = (Ul(IDN)*Cl**2 - Bx**2)/( Ul(IDN)*Cl*Cml - Bx**2)
             Ulst(IB2) = Ul(IB2)*tmp
             Ulst(IB3) = Ul(IB3)*tmp
         endif

         v_dot_B_stl = ( Ulst(IV1)*Bx + Ulst(IV2)*Ulst(IB2) + Ulst(IV3)*Ulst(IB3) )*Ulst_d_inv
         Ulst(IEN) = ( Cl*Ul(IEN) - ptotl*Ql(IV1) + ptot_st*S2 &
                     + Bx*( Ql(IV1)*Bx + Ql(IV2)*Ul(IB2) + Ql(IV3)*Ul(IB3) - v_dot_B_stl) )*Cml_inv

         Urst(IV1) = Urst(IDN)*S2
         if( dabs( Ur(IDN)*Cr*Cmr-Bx**2) < 1.0d-8*ptot_st ) then
             Urst(IV2) = Urst(IDN)*Qr(IV2)
             Urst(IV3) = Urst(IDN)*Qr(IV3)

             Urst(IB2) = Ur(IB2)
             Urst(IB3) = Ur(IB3)
         else 
             tmp = Bx*( Cr - Cmr )/(Ur(IDN)*Cr*Cmr - Bx**2)
             Urst(IV2) = Urst(IDN)*( Qr(IV2) - Ur(IB2)*tmp )
             Urst(IV3) = Urst(IDN)*( Qr(IV3) - Ur(IB3)*tmp )

             tmp = (Ur(IDN)*Cr**2 - Bx**2)/( Ur(IDN)*Cr*Cmr - Bx**2)
             Urst(IB2) = Ur(IB2)*tmp
             Urst(IB3) = Ur(IB3)*tmp
         endif

         v_dot_B_str = ( Urst(IV1)*Bx + Urst(IV2)*Urst(IB2) + Urst(IV3)*Urst(IB3) )*Urst_d_inv
         Urst(IEN) = ( Cr*Ur(IEN) - ptotr*Qr(IV1) + ptot_st*S2 &
                     + Bx*( Qr(IV1)*Bx + Qr(IV2)*Ur(IB2) + Qr(IV3)*Ur(IB3) - v_dot_B_str) )*Cmr_inv
       
         if( 0.5d0*Bx**2 < 1.0d-8*ptot_st )then
             Uldst(:) = Ulst(:)
             Urdst(:) = Urst(:)
         else 
             sum_sqrtd_inv = 1.0d0/(sqrtdl + sqrtdr) 
             if (Bx > 0.0d0 ) then 
                 bxsgn = 1.0d0
             else 
                 bxsgn = -1.0d0
             endif

             Uldst(IDN) = Ulst(IDN)
             Urdst(IDN) = Urst(IDN)

             Uldst(IV1) = Ulst(IV1)
             Urdst(IV1) = Urst(IV1)

             tmp = sum_sqrtd_inv*(  sqrtdl*(Ulst(IV2)*Ulst_d_inv) + sqrtdr*(Urst(IV2)*Urst_d_inv) &
                                  + bxsgn*(Urst(IB2) - Ulst(IB2)) )
             Uldst(IV2) = Uldst(IDN)*tmp
             Urdst(IV2) = Urdst(IDN)*tmp
!

             tmp = sum_sqrtd_inv*(  sqrtdl*(Ulst(IV3)*Ulst_d_inv) + sqrtdr*(Urst(IV3)*Urst_d_inv) &
                                  + bxsgn*(Urst(IB3) - Ulst(IB3)) )
             Uldst(IV3) = Uldst(IDN)*tmp
             Urdst(IV3) = Urdst(IDN)*tmp

             tmp = sum_sqrtd_inv*(  sqrtdl*Urst(IB2) + sqrtdr*Ulst(IB2) &
                       + bxsgn*sqrtdl*sqrtdr*( (Urst(IV2)*Urst_d_inv) - (Ulst(IV2)*Ulst_d_inv) ) )
             Uldst(IB2) = tmp
             Urdst(IB2) = tmp

             tmp = sum_sqrtd_inv*(  sqrtdl*Urst(IB3) + sqrtdr*Ulst(IB3) &
                       + bxsgn*sqrtdl*sqrtdr*( (Urst(IV3)*Urst_d_inv) - (Ulst(IV3)*Ulst_d_inv) ) )
             Uldst(IB3) = tmp
             Urdst(IB3) = tmp
!
             tmp = S2*Bx + (Uldst(IV2)*Uldst(IB2) + Uldst(IV3)*Uldst(IB3))/Uldst(IDN)
             Uldst(IEN) = Ulst(IEN) - sqrtdl*bxsgn*(v_dot_B_stl - tmp)
             Urdst(IEN) = Urst(IEN) + sqrtdr*bxsgn*(v_dot_B_str - tmp)
         endif

!         test = (S4 - S3)*Urst + (S3 - S2)*Urdst + (S2 - S1)*Uldst + (S1 - S0)*Ulst - S4*Ur + S0*Ul + Fr - Fl
!         print*,test(IVperp2)
         

    !--- Step 6.  Compute flux
          if( S0 >= 0.0d0 ) then
               flx(:) = Fl(:)
          else if (S4 <= 0.0d0 ) then
               flx(:) = Fr(:)
          else  if  (S1 >= 0.0d0) then
               flx(:) = Fl(:) + S0*(Ulst(:) - Ul(:))
           else if (S2 >= 0.0d0) then
               flx(:) = Fl(:) + S0*(Ulst(:) - Ul(:)) + S1*(Uldst(:) - Ulst(:))
           else if (S3 > 0.0d0 ) then
               flx(:) = Fr(:) + S4*(Urst(:) - Ur(:)) + S3*(Urdst(:) - Urst(:))
           else 
               flx(:) = Fr(:) + S4*(Urst(:) - Ur(:)) 
           endif

      return
      end subroutine HLLD

!      subroutine HLLC(dl,vl,pl,dr,vr,pr)
!!=====================================================================
!!
!! HLLC Scheme
!!
!! Purpose
!! Calculation of Numerical Flux by HLLC method
!!
!! Reference
!!  Toro EF, Spruce M, Speares Q. (1992,1994)
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
!!----- Qave speed ---
!! sl ::  left-going fastest signal velocity
!! sr :: right-going fastest signal velocity
!! sm :: contact discontinuity velocity
!! slst ::  left-going alfven velocity
!! srst :: right-going alfven velocity
!      real(8) :: sm,sl,sr
!
!! cfl :: left-state Fast Qave velocity
!! cfr :: right-sate Fast Qave velocity
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
!        rzl = leftst(muvQ)
!        vxl = leftst(muvu)/leftst(mudn)
!        vyl = leftst(muvv)/leftst(mudn)
!        vzl = leftst(muvQ)/leftst(mudn)
!        ptl = leftst(mpre)
!
!!---- Right state
!        
!        ror = rigtst(mudn)
!        eer = rigtst(muet)
!        rxr = rigtst(muvu)
!        ryr = rigtst(muvv)
!        rzr = rigtst(muvQ)
!        vxr = rigtst(muvu)/rigtst(mudn)
!        vyr = rigtst(muvv)/rigtst(mudn)
!        vzr = rigtst(muvQ)/rigtst(mudn)
!        ptr = rigtst(mpre)
!!----- Step 1. ----------------------------------------------------------|
!! Compute Qave left & right Qave speed
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
!        frzl = leftst(mfvQ)
!
!! Right value
!! Left value
!        fror = rigtst(mfdn)
!        feer = rigtst(mfet)
!        frxr = rigtst(mfvu)
!        fryr = rigtst(mfvv) 
!        frzr = rigtst(mfvQ)
!
!!----- Step 4. ----------------------------------------------------------|
!! compute middle and alfven Qave
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
!        nflux(mrvQ) = (frzl+msl*(rzlst-rzl))*maxs1 &
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

      subroutine Output( x1a, x1b, Q, flag )
      use commons
      implicit none
      real(8), intent(in) :: x1a(:), x1b(:), Q(:,:)
      logical, intent(in) :: flag
      integer::i
      character(20),parameter::dirname="snap/"
      character(20),parameter::base="bw2nd_hlld"
      character(20),parameter::suffix=".dat"
      character(100)::filename
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

      if( time + 1.0d-14.lt. tout+dtout) return

      write(filename,'(i5.5)') nout

      filename = trim(dirname)//trim(base)//trim(filename)//trim(suffix)
      open(unitbin,file=filename,form='formatted',action="write")
      write(unitbin,"(a2,e14.8)") "# ",time
      do i=1,in-1
          write(unitbin,*) x1b(i), Q(i,IDN), Q(i,IV1), Q(i,IV2), Q(i,IV3), Q(i,IPR), Bx, &
          Q(i,IB2), Q(i,IB3)
      enddo
      close(unitbin)

      print*, "output time=", time, trim(filename)

      nout=nout+1
      tout=tout + dtout

      return
      end subroutine Output

      subroutine makedirs(outdir)
      implicit none
      character(len=*), intent(in) :: outdir
      character(len=256) command
      write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
!      write(*, *) trim(command)
      call system(command)
      end subroutine makedirs
end program main
