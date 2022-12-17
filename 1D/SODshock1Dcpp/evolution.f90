      program riemann
      implicit none

! time evolution of the exact riemann problem
! original code from bruce fryxell

! declare
      character*80       :: string
      integer            :: npts,itmax,iter,i
      integer, parameter :: n = 1000
      real*8             :: x(n),rho(n),u(n),p(n), &
                            rhol,pl,ul,rhor,pr,ur,xi,xl,xr, &
                            rho1,p1,u1,rho5,p5,u5,p40,p41,f0,tube,eps, &
                            f1,p4,error,z,c5,gmfac1,gmfac2,fact, &
                            u4,rho4,w,p3,u3,rho3,c1,c3,xsh,xcd,xft,xhd,dx

      real*8             :: energy,temper,cv,const1,const2,entropy

! for the time evoluton
      integer            :: j, nstep
      real*8             :: time, tlo, thi, tstep

! for communicating the equation of state
      real*8             :: gamma, gm1, gp1, g2
      common /csmall/       gamma, gm1, gp1, g2

! common formats
10    format('snap/t',i5.5,'.dat')
11    format ('#',1x, 10x, 'x', 12x, 'density', 8x, 'pressure', 8x, 'velocity'/)
12    format (12(2x, 1pe10.2))


! set initial conditions
! state at left of discontinuity
      rhol = 1.0d0
      pl   = 1.0d0
      ul   = 0.d0

! state at right of discontinuity
      rhor = 0.125d0
      pr   = 0.1d0
      ur   = 0.0d0

      if (ul .ne. 0. .or. ur .ne. 0.) stop 'must have ul = ur = 0.'

!  equation of state
      gamma = 1.4d0
!      gamma = 5.0d0/3.0d0
      gm1   = gamma - 1.0d0
      gp1   = gamma + 1.0d0
      g2    = 2.0d0 * gamma

! location of discontinuity at t = 0
      xi = 0.5d0

! time at which solution is desired
      tlo   = 0.0d0
      thi   = 0.4d0
      nstep = 401
      tstep = 0.0d0
      if (nstep .ne. 1) tstep = (thi - tlo) / float(nstep - 1)


! number of points in solution
      npts = 500
      if (npts .gt. n) stop 'number of points npts exceeds array size n'

! spatial domain and grid
      xl = 0.0d0
      xr = 1.0d0
      if (xr .lt. xl) stop 'xr must be greater than xl'
      dx = (xr - xl) / (npts - 1)
      do i = 1, npts
       x(i) = xl + dx * float(i - 1)
      enddo

! shock solution
       if (pl .gt. pr) then
        rho1 = rhol
        p1   = pl
        u1   = ul
        rho5 = rhor
        p5   = pr
        u5   = ur
       else
        rho1 = rhor
        p1   = pr
        u1   = ur
        rho5 = rhol
        p5   = pl
        u5   = ul
       endif


! solve for post-shock pressure by secant method
! initial guesses

      p40 = p1
      p41 = p5
      f0  = tube(p40, p1, p5, rho1, rho5)

! maximum number of iterations and maximum allowable relative error
      itmax = 20
      eps   = 1.0d-5

      do iter = 1, itmax
       f1 = tube(p41, p1, p5, rho1, rho5)
       if (f1 .eq. f0) go to 05

       p4 = p41 - (p41 - p40) * f1 / (f1 - f0)

       error = abs (p4 - p41) / p41
       if (error .lt. eps) goto 05

       p40 = p41
       p41 = p4
       f0  = f1
      enddo
      write (6,*) 'post-shock pressure iteration failed to converge'
      stop 'abnormal termination'

05    continue


! compute post-shock density and velocity
      z  = (p4 / p5 - 1.0d0)
      c5 = sqrt (gamma * p5 / rho5)

      gmfac1 = 0.5d0 * gm1 / gamma
      gmfac2 = 0.5d0 * gp1 / gamma

      fact = sqrt (1.0d0 + gmfac2 * z)

      u4 = c5 * z / (gamma * fact)
      rho4 = rho5 * (1.0d0 + gmfac2 * z) / (1.0d0 + gmfac1 * z)

! shock speed
      w = c5 * fact


! compute values at foot of rarefaction
      p3 = p4
      u3 = u4
      rho3 = rho1 * (p3 / p1)**(1.0d0 /gamma)



! begin time evolution 
      do j = 1, nstep
       time = tlo + tstep * float(j - 1)


! compute positions of waves
       if (pl .gt. pr) then
        c1 = sqrt (gamma * p1 / rho1)
        c3 = sqrt (gamma * p3 / rho3)

        xsh = xi + w * time
        xcd = xi + u3 * time
        xft = xi + (u3 - c3) * time
        xhd = xi - c1 * time

        do i = 1, npts
         if (x(i) .lt. xhd) then
          rho(i) = rho1
          p(i)   = p1
          u(i)   = u1
         else if (x(i) .lt. xft) then
          u(i)   = 2.0d0 / gp1 * (c1 + (x(i) - xi) / time)
          fact   = 1.0d0 - 0.5d0 * gm1 * u(i) / c1
          rho(i) = rho1 * fact ** (2.0d0 / gm1)
          p(i)   = p1 * fact ** (2.0d0 * gamma / gm1)
         else if (x(i) .lt. xcd) then
          rho(i) = rho3
          p(i)   = p3
          u(i)   = u3
         else if (x(i) .lt. xsh) then
          rho(i) = rho4
          p(i)   = p4
          u(i)   = u4
         else
          rho(i) = rho5
          p(i)   = p5
          u(i)   = u5
         endif
        enddo
       endif

! if pr > pl, reverse solution
       if (pr .gt. pl) then
        c1 = sqrt (gamma * p1 / rho1)
        c3 = sqrt (gamma * p3 / rho3)

        xsh = xi - w * time
        xcd = xi - u3 * time
        xft = xi - (u3 - c3) * time
        xhd = xi + c1 * time

        do i = 1, npts
         if (x(i) .lt. xsh) then
          rho(i) = rho5
          p(i)   = p5
          u(i)   = -u5
         else if (x(i) .lt. xcd) then
         rho(i) = rho4
          p(i)   = p4
          u(i)   = -u4
         else if (x(i) .lt. xft) then
          rho(i) = rho3
          p(i)   = p3
          u(i)   = -u3
         else if (x(i) .lt. xhd) then
          u(i)   = -2.0d0 / gp1 * (c1 + (xi - x(i)) / time)
         fact   = 1.0d0 + 0.5d0 * gm1 * u(i) / c1
          rho(i) = rho1 * fact ** (2.0d0 / gm1)
          p(i)   = p1 * fact ** (2.0d0 * gamma / gm1)
         else
          rho(i) = rho1
          p(i)   = p1
          u(i)   = -u1
         endif
        enddo
       endif


! write solution at this time
       write(string,10) j
       open (unit=1, file=string,status='unknown')
       write(1,'(a,1f6.3)') '# ',time
       write (1, 11)
       do i = 1, npts
        energy = p(i)/((gamma - 1.0d0)*rho(i))
        temper = energy*(gamma - 1.0d0)
        cv     = energy/temper
!        const1 = gamma/(gamma-1.0d0)
        const1 = (energy + p(i)/rho(i))/temper
        const2 = p(i)/rho(i)**gamma
        entropy = cv*log(const2)

        write (1,12) x(i), rho(i), p(i), u(i),energy, temper, cv, const1, const2, entropy

       enddo
       close (unit=1)

! end of time evolution loop
      end do

      stop 'normal termination'
      end program riemann




      real*8 function tube(p4, p1, p5, rho1, rho5)
      implicit none

! shock tube equation

! declare the pass
      real*8  :: p4,p1,p5,rho1,rho5


! for communicating the equation of state
      real*8             :: gamma, gm1, gp1, g2
      common /csmall/       gamma, gm1, gp1, g2


! local variables
      real*8  :: z,c1,c5,fact

      z = (p4 / p5 - 1.0d0)
      c1 = sqrt (gamma * p1 / rho1)
      c5 = sqrt (gamma * p5 / rho5)

      fact = gm1 / g2 * (c5 / c1) * z / sqrt (1.0d0 + gp1 / g2 * z)
      fact = (1.0d0 - fact) ** (g2 / gm1)

      tube = p1 * fact - p4

      return
      end function tube

