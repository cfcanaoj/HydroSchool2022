      module commons
      implicit none
      integer::nhy
      real(8)::time,dt
      data time / 0.0d0 /
      integer,parameter::ngridi=150
      integer,parameter::ngridj=50
      integer,parameter::mgn=3
      integer,parameter::in=ngridi+2*mgn+1
     &                  ,jn=ngridj+2*mgn+1
     &                  ,kn=1
      integer,parameter::is=mgn+1
     &                  ,js=mgn+1
     &                  ,ks=1
      integer,parameter::ie=ngridi+mgn
     &                  ,je=ngridj+mgn
     &                  ,ke=1

      real(8),parameter:: x1min=-0.75d0,x1max=0.75d0
      real(8),parameter:: x2min=-0.25d0,x2max=0.25d0
      real(8),dimension(in)::x1a,x1b
      real(8),dimension(jn)::x2a,x2b
      real(8),dimension(kn)::x3a,x3b

      real(8),dimension(in,jn,kn)::d,et,mv1,mv2,mv3
      real(8),dimension(in,jn,kn)::p,ei,v1,v2,v3,cs
      real(8),dimension(in,jn,kn)::b1,b2,b3,bp
      real(8),dimension(in,jn,kn)::gp,gp1a,gp2a

      real(8),parameter::gam=1.4d0

      end module commons
      
      module fluxmod
      use commons, only : in,jn,kn
      implicit none
      real(8):: chg

      integer,parameter::nden=1,nve1=2,nve2=3,nve3=4,nene=5,npre=6
     &                         ,nbm1=7,nbm2=8,nbm3=9,nbps=10
      integer,parameter::nhyd=10
      real(8),dimension(nhyd,in,jn,kn):: svc

      integer,parameter::mudn= 1,muvu= 2,muvv= 3,muvw= 4,muet= 5
     &                          ,mubu= 6,mubv= 7,mubw= 8,mubp= 9
     &                  ,mfdn=10,mfvu=11,mfvv=12,mfvw=13,mfet=14
     &                          ,mfbu=15,mfbv=16,mfbw=17,mfbp=18
     &                          ,mcsp=19,mvel=20,mpre=21
      integer,parameter:: mflx=9,madd=3

      integer,parameter:: mden=1,mrv1=2,mrv2=3,mrv3=4,meto=5
     &                          ,mrvu=muvu,mrvv=muvv,mrvw=muvw
     &                          ,mbm1=6,mbm2=7,mbm3=8,mbps=9
     &                          ,mbmu=mubu,mbmv=mubv,mbmw=mubw
      real(8),dimension(mflx,in,jn,kn):: nflux1,nflux2,nflux3
      real(8),dimension(in,jn,kn):: grvsrc1,grvsrc2,grvsrc3

      logical:: MP5on
      data MP5on / .true. /

      real(8),parameter :: Alpha = 1.0d0, BC2 = 1.0d0/3d0

      end module fluxmod

      program main
      use commons
      implicit none
      logical::is_final
      data is_final /.false./

      write(6,*) "setup grids and fields"
      call GenerateGrid
      call GenerateProblem
      call ConsvVariable
      write(6,*) "entering main loop"
! main loop
      mainloop: do nhy=1,100000
         if(mod(nhy,100) .eq. 0 )write(6,*)nhy,time,dt
         call TimestepControl
         call BoundaryCondition
         call StateVevtor
         call EvaulateCh
         call GravForce
         call NumericalFlux1
         call NumericalFlux2
         call UpdateConsv
         call DampPsi
         call PrimVariable
         time=time+dt
         call Output(is_final)
         if (time .gt. 20.0d0 ) exit mainloop
      enddo mainloop
      is_final = .true.
      call Output(is_final)

      write(6,*) "program has been finished"
      end program main

      subroutine GenerateGrid
      use commons
      implicit none
      real(8)::dx,dy
      integer::i,j,k
! x coordinates
      dx=(x1max-x1min)/dble(ngridi)
      do i=1,in
         x1a(i) = dx*(i-(mgn+1))+x1min
      enddo
      do i=1,in-1
         x1b(i) = 0.5d0*(x1a(i+1)+x1a(i))
      enddo
 
! y coordinates
      dy=(x2max-x2min)/dble(ngridj)
      do j=1,jn
         x2a(j) = dy*(j-(mgn+1))+x2min
      enddo

      do j=1,jn-1
         x2b(j) = 0.5d0*(x2a(j+1)+x2a(j))
      enddo

      return
      end subroutine GenerateGrid

      subroutine GenerateProblem
      use commons
      use fluxmod, only:MP5on
      implicit none
      integer::i,j,k

      real(8) :: RU,RD,Gra,P0,Rcenter
      data RU  / 2.0d0 /
      data RD  / 1.0d0 /
      data GRA / 0.1d0 /
      data P0  / 2.5d0 /
      data Rcenter / 0.0d0 /

      real(8)::d1d(in),pa1d(in),pb1d(in)

      real(8)::pi
      pi=acos(-1.0d0)

      d(:,:,:) = 1.0d0

      do i=1,in-1
         if(x1b(i) .lt. Rcenter) then
            d1d(i) = RD   
         else
            d1d(i) = RU
         endif
      enddo

      pa1d(1:is) = P0 - GRA * RD * x1a(1:is) ! x1a(1) is negative

      do i=is,in-1
          pa1d(i+1) = pa1d(i) 
     &       - d1d(i)*GRA*(x1a(i+1)-x1a(i))
      enddo

      if(.not. MP5on) then

         do i=1,in-1
            pb1d(i) = 0.5d0*(pa1d(i+1)+pa1d(i))
         enddo

      else ! MP5on is true 
         do i=1,mgn
            pb1d(i) = 0.5d0*(pa1d(i+1)+pa1d(i))
         enddo

         do i=is,ie
            pb1d(i) =  1.0d0/30.0d0*pa1d(i-2)
     &               -13.0d0/60.0d0*pa1d(i-1)
     &               +47.0d0/60.0d0*pa1d(i  )
     &               + 9.0d0/20.0d0*pa1d(i+1)
     &               - 1.0d0/20.0d0*pa1d(i+2)
         enddo

         do i=ie+1,in-1
            pb1d(i) = 0.5d0*(pa1d(i+1)+pa1d(i))
         enddo
      endif ! MP5on

      do k=ks,ke
      do j=js,je
      do i=is,ie
          d(i,j,k) = d1d(i)
          p(i,j,k) = pb1d(i)
         v1(i,j,k) = 0.0d0
         v2(i,j,k) = 0.0d0
         v3(i,j,k) = 0.0d0
         b1(i,j,k) = 0.0d0
         b2(i,j,k) = 0.0d0
         b3(i,j,k) = 0.0d0
         bp(i,j,k) = 0.0d0
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=js,je
      do i=1,in-1
         gp(i,j,k) = -GRA*x1b(i)
      enddo
      enddo
      enddo

! pert
      do k=ks,ke
      do j=js,je
      do i=is,ie
          v1(i,j,k)= 0.01d0/4.0d0
     & *(1.0d0+cos(2.0d0*pi*(x2b(j)-(x2max+x2min)/2.0d0)/(x2max-x2min)))
     & *(1.0d0+cos(2.0d0*pi*(x1b(i)-(x1max+x1min)/2.0d0)/(x1max-x1min)))
       enddo
      enddo
      enddo

      do k=ks,ke
      do j=js,je
      do i=is,ie
          ei(i,j,k) = p(i,j,k)/(gam-1.0d0)
          cs(i,j,k) = sqrt(gam*p(i,j,k)/d(i,j,k))
      enddo
      enddo
      enddo
      

      call BoundaryCondition

      return
      end subroutine GenerateProblem

      subroutine BoundaryCondition
      use commons
      implicit none
      integer::i,j,k

!reflection
      k=ks
      do j=1,jn-1
      do i=1,mgn
           d(is-i,j,k) =  d(is+i-1,j,k)
          ei(is-i,j,k) = ei(is+i-1,j,k)
          v1(is-i,j,k) = -v1(is+i-1,j,k)
          v2(is-i,j,k) = v2(is+i-1,j,k)
          v3(is-i,j,k) = v3(is+i-1,j,k)
          gp(is-i,j,k) = gp(is+i-1,j,k)
          b1(is-i,j,k) = b1(is+i-1,j,k)
          b2(is-i,j,k) = b2(is+i-1,j,k)
          b3(is-i,j,k) = b3(is+i-1,j,k)
          bp(is-i,j,k) = bp(is+i-1,j,k)
      enddo
      enddo

!reflection
      k=ks
      do j=1,jn-1
      do i=1,mgn
           d(ie+i,j,k) =  d(ie-i+1,j,k)
          ei(ie+i,j,k) = ei(ie-i+1,j,k)
          v1(ie+i,j,k) = -v1(ie-i+1,j,k)
          v2(ie+i,j,k) = v2(ie-i+1,j,k)
          v3(ie+i,j,k) = v3(ie-i+1,j,k)
          gp(ie+i,j,k) = gp(ie-i+1,j,k)
          b1(ie+i,j,k) = b1(ie-i+1,j,k)
          b2(ie+i,j,k) = b2(ie-i+1,j,k)
          b3(ie+i,j,k) = b3(ie-i+1,j,k)
          bp(ie+i,j,k) = bp(ie-i+1,j,k)
      enddo
      enddo

! periodic
      k=ks
      do i=1,in-1
      do j=1,mgn
           d(i,js-j,k) =  d(i,je-j+1,k)
          ei(i,js-j,k) = ei(i,je-j+1,k)
          v1(i,js-j,k) = v1(i,je-j+1,k)
          v2(i,js-j,k) = v2(i,je-j+1,k)
          v3(i,js-j,k) = v3(i,je-j+1,k)
          gp(i,js-j,k) = gp(i,je-j+1,k)
          b1(i,js-j,k) = b1(i,je-j+1,k)
          b2(i,js-j,k) = b2(i,je-j+1,k)
          b3(i,js-j,k) = b3(i,je-j+1,k)
          bp(i,js-j,k) = bp(i,je-j+1,k)
      enddo
      enddo

! periodic
      k=ks
      do i=1,in-1
      do j=1,mgn
           d(i,je+j,k) =  d(i,js+j-1,k)
          ei(i,je+j,k) = ei(i,js+j-1,k)
          v1(i,je+j,k) = v1(i,js+j-1,k)
          v2(i,je+j,k) = v2(i,js+j-1,k)
          v3(i,je+j,k) = v3(i,js+j-1,k)
          gp(i,je+j,k) = gp(i,js+j-1,k)
          b1(i,je+j,k) = b1(i,js+j-1,k)
          b2(i,je+j,k) = b2(i,js+j-1,k)
          b3(i,je+j,k) = b3(i,js+j-1,k)
          bp(i,je+j,k) = bp(i,js+j-1,k)
      enddo
      enddo


      return
      end subroutine BoundaryCondition

      subroutine ConsvVariable
      use commons
      implicit none
      integer::i,j,k
      do k=ks,ke
      do j=js,je
      do i=is,ie
          et(i,j,k) = 0.5d0*d(i,j,k)*(
     &                    +v1(i,j,k)**2
     &                    +v2(i,j,k)**2
     &                    +v3(i,j,k)**2)
     &               +0.5d0*(
     &                    +b1(i,j,k)**2
     &                    +b2(i,j,k)**2
     &                    +b3(i,j,k)**2)
     &                    +ei(i,j,k)
          mv1(i,j,k) =d(i,j,k)*v1(i,j,k)
          mv2(i,j,k) =d(i,j,k)*v2(i,j,k)
          mv3(i,j,k) =d(i,j,k)*v3(i,j,k)
      enddo
      enddo
      enddo
      
      return
      end subroutine Consvvariable

      subroutine PrimVariable
      use commons
      implicit none
      integer::i,j,k
      do k=ks,ke
      do j=js,je
      do i=is,ie
          v1(i,j,k) = mv1(i,j,k)/d(i,j,k)
          v2(i,j,k) = mv2(i,j,k)/d(i,j,k)
          v3(i,j,k) = mv3(i,j,k)/d(i,j,k)

          ei(i,j,k) =  et(i,j,k)
     &          -0.5d0*d(i,j,k)*(
     &                    +v1(i,j,k)**2
     &                    +v2(i,j,k)**2
     &                    +v3(i,j,k)**2)
     &          -0.5d0*(
     &                    +b1(i,j,k)**2
     &                    +b2(i,j,k)**2
     &                    +b3(i,j,k)**2)

           p(i,j,k) =  ei(i,j,k)*(gam-1.0d0)
          cs(i,j,k) =  sqrt(gam*p(i,j,k)/d(i,j,k))
      enddo
      enddo
      enddo

      return
      end subroutine PrimVariable

      subroutine TimestepControl
      use commons
      implicit none
      real(8)::dtl1
      real(8)::dtl2
      real(8)::dtl3
      real(8)::dtlocal
      real(8)::dtmin
      real(8)::ctot
      integer::i,j,k
      dtmin=1.0d90
      do k=ks,ke
      do j=js,je
      do i=is,ie
         ctot = sqrt(cs(i,j,k)**2
     &            +( b1(i,j,k)**2
     &              +b2(i,j,k)**2
     &              +b3(i,j,k)**2
     &                           )/d(i,j,k))
         dtl1 =(x1a(i+1)-x1a(i))/(abs(v1(i,j,k)) + ctot)
         dtl2 =(x2a(j+1)-x2a(j))/(abs(v2(i,j,k)) + ctot)
!         dtl3 =(x3a(k+1)-x3a(k))/(abs(v3(i,j,k)) + ctot)
         dtlocal = min (dtl1,dtl2)
         if(dtlocal .lt. dtmin) dtmin = dtlocal
      enddo
      enddo
      enddo

      dt = 0.03d0 * dtmin

      return
      end subroutine TimestepControl

      subroutine StateVevtor
      use commons
      use fluxmod
      implicit none
      integer::i,j,k

      k=ks
      do j=1,jn-1
      do i=1,in-1
         svc(nden,i,j,k) =  d(i,j,k)
         svc(nve1,i,j,k) = v1(i,j,k)
         svc(nve2,i,j,k) = v2(i,j,k)
         svc(nve3,i,j,k) = v3(i,j,k)
         svc(nene,i,j,k) = ei(i,j,k)/d(i,j,k)
         svc(npre,i,j,k) = ei(i,j,k)*(gam-1.0d0)
         svc(nbm1,i,j,k) = b1(i,j,k)
         svc(nbm2,i,j,k) = b2(i,j,k)
         svc(nbm3,i,j,k) = b3(i,j,k)
         svc(nbps,i,j,k) = bp(i,j,k)
         p(i,j,k) =  ei(i,j,k)*(gam-1.0d0) ! for output boundary 
      enddo
      enddo

      return
      end subroutine StateVevtor

      subroutine minmod(a,b,d)
      use fluxmod, only : nhyd
      implicit none
      real(8),dimension(nhyd),intent(in)::a,b
      real(8),dimension(nhyd),intent(out)::d
      integer:: n

      do n=1,nhyd
         d(n) = sign(1.0d0,a(n))*max(0.0d0,min(abs(a(n))
     &                                        ,sign(1.0d0,a(n))*b(n)))
      enddo

      return
      end subroutine minmod


      subroutine vanLeer(dvp,dvm,dv)
      use fluxmod, only : nhyd
      implicit none
      real(8),dimension(nhyd),intent(in)::dvp,dvm
      real(8),dimension(nhyd),intent(out)::dv
      integer:: n

      do n=1,nhyd
         if(dvp(n)*dvm(n) .gt. 0.0d0)then
            dv(n) =2.0d0*dvp(n)*dvm(n)/(dvp(n)+dvm(n))
         else
            dv(n) = 0.0d0
         endif

      enddo

      return
      end subroutine vanLeer

      subroutine MClimiter(a,b,c,d)
      use fluxmod, only : nhyd
      implicit none
      real(8),dimension(nhyd),intent(in)::a,b,c
      real(8),dimension(nhyd),intent(out)::d
      integer:: n

      do n=1,nhyd
         d(n) = sign(1.0d0,a(n))*max(0.0d0,min(abs(a(n))
     &                                  ,sign(1.0d0,a(n))*b(n)
     &                                  ,sign(1.0d0,a(n))*c(n)))
      enddo

      return
      end subroutine MClimiter

      subroutine mc(dql,dqr,dqc,dq)
      real(8), intent(in) :: dql,dqr,dqc
      real(8), intent(out):: dq
      real(8):: tmp
      call minmod2(dql,dqr,tmp)
      call minmod2(tmp,dqc,dq)
      return
      end subroutine mc

      subroutine minmod2(dql,dqr,dq)
      real(8), intent(in) :: dql,dqr
      real(8), intent(out):: dq

      dq = 0.5d0*(sign(1.0d0,dql)+sign(1.0d0,dqr))*min(dabs(dql),dabs(dqr))
  
      return
      end subroutine minmod2

      subroutine minmod4(d1,d2,d3,d4,d)
      real(8), intent(in) :: d1,d2,d3,d4
      real(8), intent(out) :: d
      real(8) :: sign1 
      sign1 = sign(1d0,d1)
      d = 0.125d0*(sign1+sign(1.0d0,d2)) 
     & *dabs( (sign1+sign(1.0d0,d3))
     & *(sign1+sign(1.0d0,d4)))
     & *min(dabs(d1),dabs(d2),dabs(d3),dabs(d4))
      return
      end subroutine minmod4

      subroutine NumericalFlux1
      use commons, only: is,ie,in,js,je,jn,ks,ke,kn,gam,x1b,x1a
      use fluxmod
      implicit none
      integer::i,j,k,n
      real(8),dimension(nhyd):: dsvp,dsvm,dsvc,dsv
      real(8),dimension(nhyd,in,jn,kn):: leftpr,rigtpr
      real(8),dimension(2*mflx+madd,in,jn,kn):: leftco,rigtco
      real(8),dimension(2*mflx+madd):: leftst,rigtst
      real(8),dimension(mflx):: nflux
      real(8):: ptl,css,cts
! Variables for MP5
      real(8),dimension(nhyd,5) :: wwc
      real(8) :: wwor,djm1,dj,djp1,dm4jph,dm4jmh,qqul,qqmd,qqlc
      real(8) :: qqmin,qqmax,qqlr,qqll
      real(8) :: dqm,dqp,dqc,dq
      real(8),dimension(5,2) :: ccx
      
      ccx(1,1) =  -1.0d0/20.0d0 ! for i-1/2,R
      ccx(2,1) =   9.0d0/20.0d0
      ccx(3,1) =  47.0d0/60.0d0
      ccx(4,1) = -13.0d0/60.0d0
      ccx(5,1) =   1.0d0/30.0d0

      ccx(1,2) =   1.0d0/30.0d0 ! for i+1/2,L
      ccx(2,2) = -13.0d0/60.0d0
      ccx(3,2) =  47.0d0/60.0d0
      ccx(4,2) =   9.0d0/20.0d0
      ccx(5,2) =  -1.0d0/20.0d0

      if(.not. MP5on) then
      k=ks
      do j=js,je
      do i=is-1,ie+1
         dsvp(:) = (svc(:,i+1,j,k) -svc(:,i,j,k)                 )
         dsvm(:) = (                svc(:,i,j,k) - svc(:,i-1,j,k))

         call vanLeer(dsvp,dsvm,dsv)
!         call minmod(dsvp,dsvm,dsv)
         leftpr(:,i+1,j,k) = svc(:,i,j,k) + 0.5d0*dsv(:)
         rigtpr(:,i  ,j,k) = svc(:,i,j,k) - 0.5d0*dsv(:)
      enddo
      enddo

      else

      k=ks
      do j=js,je
      do i=is-1,ie+1
         wwc(1:nhyd,1:5) = svc(1:nhyd,i-2:i+2,j,k)
         do n=1,nhyd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! left 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         wwor = ccx(1,2)*wwc(n,1) + ccx(2,2)*wwc(n,2)
     &        + ccx(3,2)*wwc(n,3) + ccx(4,2)*wwc(n,4)
     &        + ccx(5,2)*wwc(n,5)
         dqp = (wwc(n,4)-wwc(n,3))
         dqm = (wwc(n,3)-wwc(n,2))
         dqc = (wwc(n,4)-wwc(n,2))/2.0d0
         
         call minmod2(dqp,dqm,dq)
         qqll = wwc(n,3) + 0.5d0*dq
         call mc(2.0d0*dqp,2.0d0*dqm,dqc,dq)
         qqul = wwc(n,3) + 0.5d0*dq
!         qqul = wwc(n,3) + Alpha*dqm

         djm1 = wwc(n,1) -2d0*wwc(n,2) + wwc(n,3)
         dj   = wwc(n,2) -2d0*wwc(n,3) + wwc(n,4)
         djp1 = wwc(n,3) -2d0*wwc(n,4) + wwc(n,5)           

         call minmod4(4.0d0*dj-djp1,4.0d0*djp1-dj,dj,djp1,dm4jph)
         call minmod4(4.0d0*dj-djm1,4.0d0*djm1-dj,dj,djm1,dm4jmh)

         qqmd = 0.5d0*(wwc(n,3)+wwc(n,4) - dm4jph)
         qqlc = wwc(n,3) + 0.5d0*(wwc(n,3)-wwc(n,2)) + BC2*dm4jmh

!         if(dm4jmh < 1d-20) qqlc = qqll
         qqmin = max(min(wwc(n,3),wwc(n,4),qqmd),min(qqll,qqul,qqlc))
         qqmax = min(max(wwc(n,3),wwc(n,4),qqmd),max(qqll,qqul,qqlc))
! left side
         call minmod2((qqmin-wwor),(qqmax-wwor),dq)
         leftpr(n,i+1,j,k) =  wwor + dq 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Right
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         wwor = ccx(5,1)*wwc(n,5) + ccx(4,1)*wwc(n,4)
     &        + ccx(3,1)*wwc(n,3) + ccx(2,1)*wwc(n,2)
     &        + ccx(1,1)*wwc(n,1)

         dqp = (wwc(n,4)-wwc(n,3))
         dqm = (wwc(n,3)-wwc(n,2))
         dqc = (wwc(n,4)-wwc(n,2))/2.0d0

         call minmod2(dqp,dqm,dq)
         qqll = wwc(n,3) - 0.5d0*dq 
         call mc(2.0d0*dqp,2.0d0*dqm,dqc,dq)
         qqlr = wwc(n,3) - 0.5d0* dq        ! the last term is netgative
!         qqlr = wwc(n,3) - Alpha*dqp            ! the last term is netgative

         qqmd = 0.5d0*(wwc(n,3)+wwc(n,2) - dm4jmh)
         qqlc = wwc(n,3) + 0.5d0*(wwc(n,3)-wwc(n,4)) + BC2*dm4jph
               
!         if(dm4jph < 1d-20) qqlc = qqll
         qqmin = max(min(wwc(n,3),wwc(n,2),qqmd),min(qqll,qqlr,qqlc))
         qqmax = min(max(wwc(n,3),wwc(n,2),qqmd),max(qqll,qqlr,qqlc))
         
         call minmod2((qqmin-wwor),(qqmax-wwor),dq)
         rigtpr(n,i  ,j,k) = wwor + dq
         enddo
      enddo
      enddo

      endif

      do j=js,je
      do i=is,ie+1
!====================
! Left
!====================

! Consvative variables
         leftco(mudn,i,j,k)=leftpr(nden,i,j,k) ! rho
         leftco(muvu,i,j,k)=leftpr(nve1,i,j,k)*leftpr(nden,i,j,k)  ! rho v_x
         leftco(muvv,i,j,k)=leftpr(nve2,i,j,k)*leftpr(nden,i,j,k)  ! rho v_y
         leftco(muvw,i,j,k)=leftpr(nve3,i,j,k)*leftpr(nden,i,j,k)  ! rho v_z
         leftco(muet,i,j,k)=leftpr(nene,i,j,k)*leftpr(nden,i,j,k)  ! e_i
     &               +0.5d0*leftpr(nden,i,j,k)*(
     &                     +leftpr(nve1,i,j,k)**2
     &                     +leftpr(nve2,i,j,k)**2
     &                     +leftpr(nve3,i,j,k)**2)                 ! + rho v^2/2
     &               +0.5d0*                    (
     &                     +leftpr(nbm1,i,j,k)**2
     &                     +leftpr(nbm2,i,j,k)**2
     &                     +leftpr(nbm3,i,j,k)**2)                 ! + B^2/2

         leftco(mubu,i,j,k)=leftpr(nbm1,i,j,k)  ! b_x
         leftco(mubv,i,j,k)=leftpr(nbm2,i,j,k)  ! b_y
         leftco(mubw,i,j,k)=leftpr(nbm3,i,j,k)  ! b_z
         leftco(mubp,i,j,k)=leftpr(nbps,i,j,k)  ! psi

! Flux
         ptl = leftpr(npre,i,j,k) + ( leftpr(nbm1,i,j,k)**2
     &                               +leftpr(nbm2,i,j,k)**2
     &                               +leftpr(nbm3,i,j,k)**2)/2.0d0 

         leftco(mfdn,i,j,k)=leftpr(nden,i,j,k)                   *leftpr(nve1,i,j,k)
         leftco(mfvu,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve1,i,j,k)*leftpr(nve1,i,j,k)
     &                     +ptl-leftpr(nbm1,i,j,k)**2
         leftco(mfvv,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve2,i,j,k)*leftpr(nve1,i,j,k)
     &                                        -leftpr(nbm2,i,j,k)*leftpr(nbm1,i,j,k)
         leftco(mfvw,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve3,i,j,k)*leftpr(nve1,i,j,k)
     &                                        -leftpr(nbm3,i,j,k)*leftpr(nbm1,i,j,k)
         leftco(mfet,i,j,k)=(leftco(muet,i,j,k)+ptl)*leftpr(nve1,i,j,k)
     &                     -( leftpr(nbm1,i,j,k)*leftpr(nve1,i,j,k)
     &                       +leftpr(nbm2,i,j,k)*leftpr(nve2,i,j,k)
     &                       +leftpr(nbm3,i,j,k)*leftpr(nve3,i,j,k))*leftpr(nbm1,i,j,k)

         leftco(mfbu,i,j,k) =  0.0d0
         leftco(mfbv,i,j,k) =  leftpr(nbm2,i,j,k)*leftpr(nve1,i,j,k)
     &                        -leftpr(nve2,i,j,k)*leftpr(nbm1,i,j,k)
         leftco(mfbw,i,j,k) =  leftpr(nbm3,i,j,k)*leftpr(nve1,i,j,k)
     &                        -leftpr(nve3,i,j,k)*leftpr(nbm1,i,j,k)
         leftco(mfbp,i,j,k) = 0.0d0  ! psi
     
         css =(gam*(gam-1.0d0)*leftpr(nene,i,j,k))
         cts =  css  !c_s^2*c_a^2
     &                       +( leftpr(nbm1,i,j,k)**2
     &                         +leftpr(nbm2,i,j,k)**2
     &                         +leftpr(nbm3,i,j,k)**2)/leftpr(nden,i,j,k) 

         leftco(mcsp,i,j,k)= sqrt((cts +sqrt(cts**2 
     &                             -4.0d0*css*leftpr(nbm1,i,j,k)**2
     &                                       /leftpr(nden,i,j,k))
     &                            )/2.0d0)
         leftco(mvel,i,j,k)= leftpr(nve1,i,j,k)
         leftco(mpre,i,j,k)= ptl
!====================
! Right
!====================
! Consvative variables
         rigtco(mudn,i,j,k)=rigtpr(nden,i,j,k) ! rho
         rigtco(muvu,i,j,k)=rigtpr(nve1,i,j,k)*rigtpr(nden,i,j,k)  ! rho v_x
         rigtco(muvv,i,j,k)=rigtpr(nve2,i,j,k)*rigtpr(nden,i,j,k)  ! rho v_y
         rigtco(muvw,i,j,k)=rigtpr(nve3,i,j,k)*rigtpr(nden,i,j,k)  ! rho v_z
         rigtco(muet,i,j,k)=rigtpr(nene,i,j,k)*rigtpr(nden,i,j,k)  ! e_i
     &               +0.5d0*rigtpr(nden,i,j,k)*(
     &                     +rigtpr(nve1,i,j,k)**2
     &                     +rigtpr(nve2,i,j,k)**2
     &                     +rigtpr(nve3,i,j,k)**2)                 ! + rho v^2/2
     &               +0.5d0*                    (
     &                     +rigtpr(nbm1,i,j,k)**2
     &                     +rigtpr(nbm2,i,j,k)**2
     &                     +rigtpr(nbm3,i,j,k)**2)                 ! + B^2/2

         rigtco(mubu,i,j,k)=rigtpr(nbm1,i,j,k)  ! b_x
         rigtco(mubv,i,j,k)=rigtpr(nbm2,i,j,k)  ! b_y
         rigtco(mubw,i,j,k)=rigtpr(nbm3,i,j,k)  ! b_z
         rigtco(mubp,i,j,k)=rigtpr(nbps,i,j,k)  ! psi

! Flux
         ptl = rigtpr(npre,i,j,k) + ( rigtpr(nbm1,i,j,k)**2
     &                               +rigtpr(nbm2,i,j,k)**2
     &                               +rigtpr(nbm3,i,j,k)**2)/2.0d0 

         rigtco(mfdn,i,j,k)=rigtpr(nden,i,j,k)                   *rigtpr(nve1,i,j,k)
         rigtco(mfvu,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve1,i,j,k)*rigtpr(nve1,i,j,k)
     &                     +ptl-rigtpr(nbm1,i,j,k)**2
         rigtco(mfvv,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve2,i,j,k)*rigtpr(nve1,i,j,k)
     &                                        -rigtpr(nbm2,i,j,k)*rigtpr(nbm1,i,j,k)
         rigtco(mfvw,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve3,i,j,k)*rigtpr(nve1,i,j,k)
     &                                        -rigtpr(nbm3,i,j,k)*rigtpr(nbm1,i,j,k)
         rigtco(mfet,i,j,k)=(rigtco(muet,i,j,k)+ptl)*rigtpr(nve1,i,j,k)
     &                     -( rigtpr(nbm1,i,j,k)*rigtpr(nve1,i,j,k)
     &                       +rigtpr(nbm2,i,j,k)*rigtpr(nve2,i,j,k)
     &                       +rigtpr(nbm3,i,j,k)*rigtpr(nve3,i,j,k))*rigtpr(nbm1,i,j,k)
     
         rigtco(mfbu,i,j,k) =  0.0d0
         rigtco(mfbv,i,j,k) =  rigtpr(nbm2,i,j,k)*rigtpr(nve1,i,j,k)
     &                        -rigtpr(nve2,i,j,k)*rigtpr(nbm1,i,j,k)
         rigtco(mfbw,i,j,k) =  rigtpr(nbm3,i,j,k)*rigtpr(nve1,i,j,k)
     &                        -rigtpr(nve3,i,j,k)*rigtpr(nbm1,i,j,k)
         rigtco(mfbp,i,j,k) = 0.0d0  ! b_z
         css =(gam*(gam-1.0d0)*rigtpr(nene,i,j,k))
         cts =  css  !c_s^2*c_a^2
     &                       +( rigtpr(nbm1,i,j,k)**2
     &                         +rigtpr(nbm2,i,j,k)**2
     &                         +rigtpr(nbm3,i,j,k)**2)/rigtpr(nden,i,j,k) 

         rigtco(mcsp,i,j,k)= sqrt((cts +sqrt(cts**2 
     &                             -4.0d0*css*rigtpr(nbm1,i,j,k)**2
     &                                       /rigtpr(nden,i,j,k))
     &                            )/2.0d0)
         rigtco(mvel,i,j,k)= rigtpr(nve1,i,j,k)
         rigtco(mpre,i,j,k)= ptl

      enddo
      enddo

      do j=js,je
      do i=is,ie+1
         leftst(:)=leftco(:,i,j,k)
         rigtst(:)=rigtco(:,i,j,k)
!         call HLLE(leftst,rigtst,nflux)
         call HLLC(leftst,rigtst,nflux)
         nflux1(mden,i,j,k)=nflux(mden)
         nflux1(mrv1,i,j,k)=nflux(mrvu)
         nflux1(mrv2,i,j,k)=nflux(mrvv)
         nflux1(mrv3,i,j,k)=nflux(mrvw)
         nflux1(meto,i,j,k)=nflux(meto)
         nflux1(mbm1,i,j,k)=nflux(mbmu)
         nflux1(mbm2,i,j,k)=nflux(mbmv)
         nflux1(mbm3,i,j,k)=nflux(mbmw)

         nflux1(mbm1,i,j,k) =  0.5d0*(leftst(mubp)+rigtst(mubp))
     &                    -0.5d0*chg*(rigtst(mubu)-leftst(mubu))        ! finite volume
         nflux1(mbps,i,j,k) = (0.5d0*(leftst(mubu)+rigtst(mubu))
     &                    -0.5d0/chg*(rigtst(mubp)-leftst(mubp)))*chg**2 ! finite volume

!         write(6,*) "bpf1",nflux1(mbps,i,j,k)

      enddo
      enddo

      return
      end subroutine Numericalflux1

      subroutine NumericalFlux2
      use commons, only: is,ie,in,js,je,jn,ks,ke,kn,gam,x2a,x2b
      use fluxmod
      implicit none
      integer::i,j,k,n
      real(8),dimension(nhyd):: dsvp,dsvm,dsvc,dsv
      real(8),dimension(nhyd,in,jn,kn):: leftpr,rigtpr
      real(8),dimension(2*mflx+madd,in,jn,kn):: leftco,rigtco
      real(8),dimension(2*mflx+madd):: leftst,rigtst
      real(8),dimension(mflx):: nflux
      real(8):: ptl,css,cts
! Variables for MP5
      real(8),dimension(nhyd,5) :: wwc
      real(8) :: wwor,djm1,dj,djp1,dm4jph,dm4jmh,qqul,qqmd,qqlc
      real(8) :: qqmin,qqmax,qqlr,qqll
      real(8) :: dqm,dqp,dqc,dq
      real(8),dimension(5,2) :: ccy
      
      ccy(1,1) =  -1.0d0/20.0d0
      ccy(2,1) =   9.0d0/20.0d0
      ccy(3,1) =  47.0d0/60.0d0
      ccy(4,1) = -13.0d0/60.0d0
      ccy(5,1) =   1.0d0/30.0d0

      ccy(1,2) =   1.0d0/30.0d0
      ccy(2,2) = -13.0d0/60.0d0
      ccy(3,2) =  47.0d0/60.0d0
      ccy(4,2) =   9.0d0/20.0d0
      ccy(5,2) =  -1.0d0/20.0d0

      if(.not. MP5on) then 

      k=ks
      do i=is,ie
      do j=js-1,je+1
         dsvp(:) = (svc(:,i,j+1,k) -svc(:,i,j,k)                 )
         dsvm(:) = (                svc(:,i,j,k) - svc(:,i,j-1,k))

         call vanLeer(dsvp,dsvm,dsv)
!         call minmod(dsvp,dsvm,dsv)
         leftpr(:,i,j+1,k) = svc(:,i,j,k) + 0.5d0*dsv(:)
         rigtpr(:,i,j  ,k) = svc(:,i,j,k) - 0.5d0*dsv(:)

!         leftpr(:,i,j,k) = svc(:,i,j-1,k)
!         rigtpr(:,i,j,k) = svc(:,i,j  ,k)

       enddo
       enddo  

      else

      k=ks
      do i=is,ie
      do j=js-1,je+1
         wwc(1:nhyd,1:5) = svc(1:nhyd,i,j-2:j+2,k)
         do n=1,nhyd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! left 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         wwor = ccy(1,2)*wwc(n,1) + ccy(2,2)*wwc(n,2)
     &        + ccy(3,2)*wwc(n,3) + ccy(4,2)*wwc(n,4)
     &        + ccy(5,2)*wwc(n,5)

         dqp = (wwc(n,4)-wwc(n,3))
         dqm = (wwc(n,3)-wwc(n,2))
         dqc = (wwc(n,4)-wwc(n,2))/2.0d0
         call minmod2(dqp,dqm,dq)
         qqll = wwc(n,3) + 0.5d0*dq
         call mc(2.0d0*dqp,2.0d0*dqm,dqc,dq)
         qqul = wwc(n,3) + 0.5d0*dq
!         qqul = wwc(n,3) + Alpha*dqm

         djm1 = wwc(n,1) -2d0*wwc(n,2) + wwc(n,3)
         dj   = wwc(n,2) -2d0*wwc(n,3) + wwc(n,4)
         djp1 = wwc(n,3) -2d0*wwc(n,4) + wwc(n,5)           

         call minmod4(4.0d0*dj-djp1,4.0d0*djp1-dj,dj,djp1,dm4jph)
         call minmod4(4.0d0*dj-djm1,4.0d0*djm1-dj,dj,djm1,dm4jmh)

         qqmd = 0.5d0*(wwc(n,3)+wwc(n,4) - dm4jph)
         qqlc = wwc(n,3) + 0.5d0*(wwc(n,3)-wwc(n,2)) + BC2*dm4jmh

!         if(dm4jmh < 1d-20) qqlc = qqll
         qqmin = max(min(wwc(n,3),wwc(n,4),qqmd),min(qqll,qqul,qqlc))
         qqmax = min(max(wwc(n,3),wwc(n,4),qqmd),max(qqll,qqul,qqlc))
! left side
         call minmod2((qqmin-wwor),(qqmax-wwor),dq) 
         leftpr(n,i,j+1,k) =  wwor + dq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Right
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         wwor = ccy(5,1)*wwc(n,5) + ccy(4,1)*wwc(n,4)
     &        + ccy(3,1)*wwc(n,3) + ccy(2,1)*wwc(n,2)
     &        + ccy(1,1)*wwc(n,1)

         dqp = (wwc(n,4)-wwc(n,3))
         dqm = (wwc(n,3)-wwc(n,2))
         dqc = (wwc(n,4)-wwc(n,2))/2.0d0
         call  minmod2(dqp,dqm,dq)
         qqll = wwc(n,3) - 0.5d0*dq
         call mc(2.0d0*dqp,2.0d0*dqm,dqc,dq)
         qqlr = wwc(n,3) - 0.5d0*dq
!         qqlr = wwc(n,3) - Alpha*dqp  ! the last term is negative

         qqmd = 0.5d0*(wwc(n,3)+wwc(n,2) - dm4jmh)
         qqlc = wwc(n,3) + 0.5d0*(wwc(n,3)-wwc(n,4)) + BC2*dm4jph
               
!         if(dm4jph < 1d-20) qqlc = qqll
         qqmin = max(min(wwc(n,3),wwc(n,2),qqmd),min(qqll,qqlr,qqlc))
         qqmax = min(max(wwc(n,3),wwc(n,2),qqmd),max(qqll,qqlr,qqlc))
         call minmod2((qqmin-wwor),(qqmax-wwor),dq)
         rigtpr(n,i  ,j,k) = wwor + dq
         enddo
      enddo
      enddo

      endif


      do i=is,ie
      do j=js,je+1
         leftco(mudn,i,j,k)=leftpr(nden,i,j,k)
         leftco(muvw,i,j,k)=leftpr(nve1,i,j,k)*leftpr(nden,i,j,k)
         leftco(muvu,i,j,k)=leftpr(nve2,i,j,k)*leftpr(nden,i,j,k) ! rho v
         leftco(muvv,i,j,k)=leftpr(nve3,i,j,k)*leftpr(nden,i,j,k)
         leftco(muet,i,j,k)=leftpr(nene,i,j,k)*leftpr(nden,i,j,k) ! internal
     &               +0.5d0*leftpr(nden,i,j,k)*(
     &                     +leftpr(nve1,i,j,k)**2
     &                     +leftpr(nve2,i,j,k)**2
     &                     +leftpr(nve3,i,j,k)**2) ! kinetic
     &               +0.5d0*                    (
     &                     +leftpr(nbm1,i,j,k)**2
     &                     +leftpr(nbm2,i,j,k)**2
     &                     +leftpr(nbm3,i,j,k)**2) ! magnetic

         leftco(mubw,i,j,k)=leftpr(nbm1,i,j,k)  ! b_x
         leftco(mubu,i,j,k)=leftpr(nbm2,i,j,k)  ! b_y
         leftco(mubv,i,j,k)=leftpr(nbm3,i,j,k)  ! b_z
         leftco(mubp,i,j,k)=leftpr(nbps,i,j,k)  ! psi

         ptl = leftpr(npre,i,j,k) + ( leftpr(nbm1,i,j,k)**2
     &                               +leftpr(nbm2,i,j,k)**2
     &                               +leftpr(nbm3,i,j,k)**2)/2.0d0 

         leftco(mfdn,i,j,k)=leftpr(nden,i,j,k)                   *leftpr(nve2,i,j,k) ! rho v
         leftco(mfvw,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve1,i,j,k)*leftpr(nve2,i,j,k)
     &                                        -leftpr(nbm1,i,j,k)*leftpr(nbm2,i,j,k)
         leftco(mfvu,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve2,i,j,k)*leftpr(nve2,i,j,k)
     &                     +ptl-leftpr(nbm2,i,j,k)**2
         leftco(mfvv,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve3,i,j,k)*leftpr(nve2,i,j,k)
     &                                        -leftpr(nbm3,i,j,k)*leftpr(nbm2,i,j,k)
         leftco(mfet,i,j,k)=(leftco(muet,i,j,k)+ptl)*leftpr(nve2,i,j,k)
     &                     -( leftpr(nbm1,i,j,k)*leftpr(nve1,i,j,k)
     &                       +leftpr(nbm2,i,j,k)*leftpr(nve2,i,j,k)
     &                       +leftpr(nbm3,i,j,k)*leftpr(nve3,i,j,k))*leftpr(nbm2,i,j,k)

         leftco(mfbw,i,j,k) =  leftpr(nbm1,i,j,k)*leftpr(nve2,i,j,k)
     &                        -leftpr(nve1,i,j,k)*leftpr(nbm2,i,j,k)
         leftco(mfbu,i,j,k) =  0.0d0
         leftco(mfbv,i,j,k) =  leftpr(nbm3,i,j,k)*leftpr(nve2,i,j,k)
     &                        -leftpr(nve3,i,j,k)*leftpr(nbm2,i,j,k)
         leftco(mfbp,i,j,k) = 0.0d0  ! psi
     
         css =(gam*(gam-1.0d0)*leftpr(nene,i,j,k))
         cts =  css  !c_s^2*c_a^2
     &                       +( leftpr(nbm1,i,j,k)**2
     &                         +leftpr(nbm2,i,j,k)**2
     &                         +leftpr(nbm3,i,j,k)**2)/leftpr(nden,i,j,k) 

         leftco(mcsp,i,j,k)= sqrt((cts +sqrt(cts**2 
     &                             -4.0d0*css*leftpr(nbm2,i,j,k)**2
     &                                          /leftpr(nden,i,j,k))
     &                            )/2.0d0)
         leftco(mvel,i,j,k)= leftpr(nve2,i,j,k)
         leftco(mpre,i,j,k)= ptl


         rigtco(mudn,i,j,k)=rigtpr(nden,i,j,k)
         rigtco(muvw,i,j,k)=rigtpr(nve1,i,j,k)*rigtpr(nden,i,j,k)
         rigtco(muvu,i,j,k)=rigtpr(nve2,i,j,k)*rigtpr(nden,i,j,k) ! rho v
         rigtco(muvv,i,j,k)=rigtpr(nve3,i,j,k)*rigtpr(nden,i,j,k)
         rigtco(muet,i,j,k)=rigtpr(nene,i,j,k)*rigtpr(nden,i,j,k) ! internal
     &               +0.5d0*rigtpr(nden,i,j,k)*(
     &                     +rigtpr(nve1,i,j,k)**2
     &                     +rigtpr(nve2,i,j,k)**2
     &                     +rigtpr(nve3,i,j,k)**2) ! kinetic
     &               +0.5d0*                    (
     &                     +rigtpr(nbm1,i,j,k)**2
     &                     +rigtpr(nbm2,i,j,k)**2
     &                     +rigtpr(nbm3,i,j,k)**2) ! magnetic

         rigtco(mubw,i,j,k)=rigtpr(nbm1,i,j,k)  ! b_x
         rigtco(mubu,i,j,k)=rigtpr(nbm2,i,j,k)  ! b_y
         rigtco(mubv,i,j,k)=rigtpr(nbm3,i,j,k)  ! b_z
         rigtco(mubp,i,j,k)=rigtpr(nbps,i,j,k)  ! psi

         ptl = rigtpr(npre,i,j,k) + ( rigtpr(nbm1,i,j,k)**2
     &                               +rigtpr(nbm2,i,j,k)**2
     &                               +rigtpr(nbm3,i,j,k)**2)/2.0d0 

         rigtco(mfdn,i,j,k)=rigtpr(nden,i,j,k)                   *rigtpr(nve2,i,j,k) ! rho v
         rigtco(mfvw,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve1,i,j,k)*rigtpr(nve2,i,j,k)
     &                     -rigtpr(nbm1,i,j,k)*rigtpr(nbm2,i,j,k)
         rigtco(mfvu,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve2,i,j,k)*rigtpr(nve2,i,j,k)
     &                     +ptl-rigtpr(nbm2,i,j,k)**2
         rigtco(mfvv,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve3,i,j,k)*rigtpr(nve2,i,j,k)
     &                     -rigtpr(nbm3,i,j,k)*rigtpr(nbm2,i,j,k)
         rigtco(mfet,i,j,k)=(rigtco(muet,i,j,k)+ptl)*rigtpr(nve2,i,j,k)
     &                     -( rigtpr(nbm1,i,j,k)*rigtpr(nve1,i,j,k)
     &                       +rigtpr(nbm2,i,j,k)*rigtpr(nve2,i,j,k)
     &                       +rigtpr(nbm3,i,j,k)*rigtpr(nve3,i,j,k))*rigtpr(nbm2,i,j,k)
         rigtco(mfbw,i,j,k) =  rigtpr(nbm1,i,j,k)*rigtpr(nve2,i,j,k)
     &                        -rigtpr(nve1,i,j,k)*rigtpr(nbm2,i,j,k)
         rigtco(mfbu,i,j,k) =  0.0d0
         rigtco(mfbv,i,j,k) =  rigtpr(nbm3,i,j,k)*rigtpr(nve2,i,j,k)
     &                        -rigtpr(nve3,i,j,k)*rigtpr(nbm2,i,j,k)
         rigtco(mfbp,i,j,k) = 0.0d0  ! psi
     
         css =(gam*(gam-1.0d0)*rigtpr(nene,i,j,k))
         cts =  css  !c_s^2*c_a^2
     &                       +( rigtpr(nbm1,i,j,k)**2
     &                         +rigtpr(nbm2,i,j,k)**2
     &                         +rigtpr(nbm3,i,j,k)**2)/rigtpr(nden,i,j,k) 

         rigtco(mcsp,i,j,k)= sqrt((cts +sqrt(cts**2 
     &                             -4.0d0*css*rigtpr(nbm2,i,j,k)**2
     &                                          /rigtpr(nden,i,j,k))
     &                            )/2.0d0)
         rigtco(mvel,i,j,k)= rigtpr(nve2,i,j,k)
         rigtco(mpre,i,j,k)= ptl

      enddo
      enddo

      do i=is,ie
      do j=js,je+1
         leftst(:)=leftco(:,i,j,k)
         rigtst(:)=rigtco(:,i,j,k)
!         call HLLE(leftst,rigtst,nflux)
         call HLLC(leftst,rigtst,nflux)
         nflux2(mden,i,j,k)=nflux(mden)
         nflux2(mrv1,i,j,k)=nflux(mrvw)
         nflux2(mrv2,i,j,k)=nflux(mrvu) ! mrv2=3, mrvu=2
         nflux2(mrv3,i,j,k)=nflux(mrvv)
         nflux2(meto,i,j,k)=nflux(meto)
         nflux2(mbm1,i,j,k)=nflux(mbmw)
         nflux2(mbm2,i,j,k)=nflux(mbmu)
         nflux2(mbm3,i,j,k)=nflux(mbmv)

         nflux2(mbm2,i,j,k) =  0.5d0*(leftst(mubp)+rigtst(mubp))
     &                    -0.5d0*chg*(rigtst(mubu)-leftst(mubu))        ! finite volume
         nflux2(mbps,i,j,k) = (0.5d0*(leftst(mubu)+rigtst(mubu))
     &                    -0.5d0/chg*(rigtst(mubp)-leftst(mubp)))*chg**2 ! finite volume
!         write(6,*) "bpf2",nflux2(mbps,i,j,k)

      enddo
      enddo

      return
      end subroutine Numericalflux2

      subroutine HLLE(leftst,rigtst,nflux)
      use fluxmod
      implicit none
      real(8),dimension(2*mflx+madd),intent(in)::leftst,rigtst
      real(8),dimension(mflx),intent(out)::nflux
      real(8),dimension(mflx)::ul,ur,fl,fr
      real(8)::csl,csr
      real(8):: vl, vr
      real(8):: sl, sr

      ul(1:mflx) = leftst(1:mflx)
      fl(1:mflx) = leftst(mflx+1:2*mflx)
      ur(1:mflx) = rigtst(1:mflx)
      fr(1:mflx) = rigtst(mflx+1:2*mflx)
      csl=leftst(mcsp)
      csr=rigtst(mcsp)
       vl=leftst(mvel)
       vr=rigtst(mvel)

       sl = min(vl,vr) - max(csl,csr)
       sl = min(0.0d0,sl)
       sr = max(vl,vr) + max(csl,csr)
       sr = max(0.0d0,sr)

       nflux(:) = (sr*fl(:)-sl*fr(:) +sl*sr*(ur(:)-ul(:)))/(sr-sl)

      return
      end subroutine HLLE

      subroutine HLLC(leftst,rigtst,nflux)
!=====================================================================
!
! HLLC Scheme
!
! Purpose
! Calculation of Numerical Flux by HLLC method
!
! Reference
!  Toro EF, Spruce M, Speares W. (1992,1994)
!
! Input
! Output
!=====================================================================
      use fluxmod, only: mflx,madd
     &                 , mudn,muvu,muvv,muvw,muet
     &                 , mfdn,mfvu,mfvv,mfvw,mfet
     &                 , mcsp,mvel,mpre
     &                 , mden,mrvu,mrvv,mrvw,meto

      implicit none
      real(8),dimension(2*mflx+madd),intent(in)::leftst,rigtst
      real(8),dimension(mflx),intent(out)::nflux

!----- U -----
! qql :: left state
! qqr :: right state
      real(8) :: rol,vxl,vyl,vzl,ptl,eel
      real(8) :: ror,vxr,vyr,vzr,ptr,eer
      real(8) :: rxl,ryl,rzl
      real(8) :: rxr,ryr,rzr
      real(8) :: ptst

!----- U* ----
! qqlst ::  left state
! qqrst :: right state
      real(8) :: rolst,vxlst,vylst,vzlst,eelst
      real(8) :: rorst,vxrst,vyrst,vzrst,eerst
      real(8) :: rxlst,rylst,rzlst
      real(8) :: rxrst,ryrst,rzrst

!----- flux ---
! fqql ::  left physical flux
! fqqr :: right physical flux
      real(8) :: frol,frxl,fryl,frzl,feel
      real(8) :: fror,frxr,fryr,frzr,feer

!----- wave speed ---
! sl ::  left-going fastest signal velocity
! sr :: right-going fastest signal velocity
! sm :: contact discontinuity velocity
! slst ::  left-going alfven velocity
! srst :: right-going alfven velocity
      real(8) :: sm,sl,sr

! cfl :: left-state Fast wave velocity
! cfr :: right-sate Fast wave velocity
      real(8) :: cfl,cfr

!--------------------
! temporary variables
      real(8) :: sdl,sdr,sdml,sdmr,isdml,isdmr,rosdl,rosdr
      real(8) :: temp
  
! no if
      real(8) :: sign1,maxs1,mins1
      real(8) :: msl,msr

!----- Step 0. ----------------------------------------------------------|

!---- Left state
        
        rol = leftst(mudn)
        eel = leftst(muet)
        rxl = leftst(muvu)
        ryl = leftst(muvv)
        rzl = leftst(muvw)
        vxl = leftst(muvu)/leftst(mudn)
        vyl = leftst(muvv)/leftst(mudn)
        vzl = leftst(muvw)/leftst(mudn)
        ptl = leftst(mpre)

!---- Right state
        
        ror = rigtst(mudn)
        eer = rigtst(muet)
        rxr = rigtst(muvu)
        ryr = rigtst(muvv)
        rzr = rigtst(muvw)
        vxr = rigtst(muvu)/rigtst(mudn)
        vyr = rigtst(muvv)/rigtst(mudn)
        vzr = rigtst(muvw)/rigtst(mudn)
        ptr = rigtst(mpre)
!----- Step 1. ----------------------------------------------------------|
! Compute wave left & right wave speed
!
         
        cfl = leftst(mcsp)
        cfr = rigtst(mcsp)

        sl = min(vxl,vxr)-max(cfl,cfr) ! note sl is negative
        sr = max(vxl,vxr)+max(cfl,cfr)
!----- Step 2. ----------------------------------------------------------|
! compute L/R fluxs
!
! Left value
        frol = leftst(mfdn)
        feel = leftst(mfet)
        frxl = leftst(mfvu)
        fryl = leftst(mfvv)
        frzl = leftst(mfvw)

! Right value
! Left value
        fror = rigtst(mfdn)
        feer = rigtst(mfet)
        frxr = rigtst(mfvu)
        fryr = rigtst(mfvv) 
        frzr = rigtst(mfvw)

!----- Step 4. ----------------------------------------------------------|
! compute middle and alfven wave
!
        sdl = sl - vxl
        sdr = sr - vxr
        rosdl = rol*sdl
        rosdr = ror*sdr

        temp = 1.0d0/(rosdr - rosdl)
! Eq. 45
        sm = (rosdr*vxr - rosdl*vxl - ptr + ptl)*temp
           
        sdml = sl - sm; isdml = 1.0d0/sdml
        sdmr = sr - sm; isdmr = 1.0d0/sdmr
        
!----- Step 5. ----------------------------------------------------------|
! compute intermediate states
!
! Eq. 49
        ptst = (rosdr*ptl-rosdl*ptr+rosdl*rosdr*(vxr-vxl))*temp

!----- Step 5A. ----------------------------------------------------------|
! compute Ul*
!

        rolst = rol*sdl   *isdml
        vxlst = sm
        rxlst = rolst*vxlst
           
        vylst = vyl
        rylst = rolst*vylst
        vzlst = vzl
        rzlst = rolst*vzlst

        eelst =(sdl*eel - ptl*vxl + ptst*sm  )*isdml

!----- Step 5B. ----------------------------------------------------------|
! compute Ur*
!

        rorst   = rosdr   *isdmr
        vxrst = sm
        rxrst = rorst*vxrst
        vyrst = vyr
        ryrst = rorst*vyrst
        vzrst = vzr
        rzrst = rorst*vzrst
           
        eerst = (sdr*eer - ptr*vxr  + ptst*sm  )*isdmr
              
!----- Step 6. ----------------------------------------------------------|
! compute flux
        sign1 = sign(1.0d0,sm)    ! 1 for sm>0, -1 for sm<0
        maxs1 =  max(0.0d0,sign1) ! 1 sm>0, 0 for sm<0
        mins1 = -min(0.0d0,sign1) ! 0 sm>0,-1 for sm<0

        msl   = min(sl  ,0.0d0)   ! 0 for sl > 0, sl for sl < 0
        msr   = max(sr  ,0.0d0)   ! S_R > 0

        nflux(mden) = (frol+msl*(rolst-rol))*maxs1
     &               +(fror+msr*(rorst-ror))*mins1
        nflux(meto) = (feel+msl*(eelst-eel))*maxs1 
     &               +(feer+msr*(eerst-eer))*mins1
        nflux(mrvu) = (frxl+msl*(rxlst-rxl))*maxs1 
     &               +(frxr+msr*(rxrst-rxr))*mins1
        nflux(mrvv) = (fryl+msl*(rylst-ryl))*maxs1 
     &               +(fryr+msr*(ryrst-ryr))*mins1
        nflux(mrvw) = (frzl+msl*(rzlst-rzl))*maxs1 
     &               +(frzr+msr*(rzrst-rzr))*mins1

      return
      end subroutine HLLC

      subroutine UpdateConsv
      use commons
      use fluxmod
      implicit none
      integer::i,j,k

      do k=ks,ke
      do j=js,je
      do i=is,ie
         
         d(i,j,k) = d(i,j,k) 
     & +dt*(
     &  (- nflux1(mden,i+1,j,k)
     &   + nflux1(mden,i  ,j,k))/(x1a(i+1)-x1a(i)) 
     & +(- nflux2(mden,i,j+1,k)
     &   + nflux2(mden,i,j  ,k))/(x2a(j+1)-x2a(j)) 
     &      )

         mv1(i,j,k) = mv1(i,j,k) 
     & +dt*(
     &    +  grvsrc1(i,j,k)
     & +(- nflux1(mrv1,i+1,j,k)
     &   + nflux1(mrv1,i  ,j,k))/(x1a(i+1)-x1a(i)) 
     & +(- nflux2(mrv1,i,j+1,k)
     &   + nflux2(mrv1,i,j  ,k))/(x2a(j+1)-x2a(j)) 
     &      )

         mv2(i,j,k) = mv2(i,j,k) 
     & +dt*(
     &    +  grvsrc2(i,j,k)
     & +(- nflux1(mrv2,i+1,j,k)
     &   + nflux1(mrv2,i  ,j,k))/(x1a(i+1)-x1a(i)) 
     & +(- nflux2(mrv2,i,j+1,k)
     &   + nflux2(mrv2,i,j  ,k))/(x2a(j+1)-x2a(j)) 
     &      )

         mv3(i,j,k) = mv3(i,j,k) 
     & +dt*(
     &    +  grvsrc3(i,j,k)
     & +(- nflux1(mrv3,i+1,j,k)
     &   + nflux1(mrv3,i  ,j,k))/(x1a(i+1)-x1a(i)) 
     & +(- nflux2(mrv3,i,j+1,k)
     &   + nflux2(mrv3,i,j  ,k))/(x2a(j+1)-x2a(j)) 
     &      )

          et(i,j,k) = et(i,j,k) 
     & +dt*(
     &  (- nflux1(meto,i+1,j,k)
     &   + nflux1(meto,i  ,j,k) 
     &   - nflux1(mden,i+1,j,k)*(gp(i,j,k)-gp1a(i+1,j,k))
     &   + nflux1(mden,i  ,j,k)*(gp(i,j,k)-gp1a(i  ,j,k))
     &  )/(x1a(i+1)-x1a(i)) 
     & +(- nflux2(meto,i,j+1,k)
     &   + nflux2(meto,i,j  ,k)
     &   - nflux2(mden,i,j+1,k)*(gp(i,j,k)-gp2a(i,j+1,k))
     &   + nflux2(mden,i,j  ,k)*(gp(i,j,k)-gp2a(i,j  ,k))
     &  )/(x2a(j+1)-x2a(j)) 
     &      )

          b1(i,j,k) = b1(i,j,k) 
     & +dt*(
     &  (- nflux1(mbm1,i+1,j,k)
     &   + nflux1(mbm1,i  ,j,k))/(x1a(i+1)-x1a(i)) 
     & +(- nflux2(mbm1,i,j+1,k)
     &   + nflux2(mbm1,i,j  ,k))/(x2a(j+1)-x2a(j)) 
     &      )

          b2(i,j,k) = b2(i,j,k) 
     & +dt*(
     &  (- nflux1(mbm2,i+1,j,k)
     &   + nflux1(mbm2,i  ,j,k))/(x1a(i+1)-x1a(i)) 
     & +(- nflux2(mbm2,i,j+1,k)
     &   + nflux2(mbm2,i,j  ,k))/(x2a(j+1)-x2a(j)) 
     &      )

          b3(i,j,k) = b3(i,j,k) 
     & +dt*(
     &  (- nflux1(mbm3,i+1,j,k)
     &   + nflux1(mbm3,i  ,j,k))/(x1a(i+1)-x1a(i)) 
     & +(- nflux2(mbm3,i,j+1,k)
     &   + nflux2(mbm3,i,j  ,k))/(x2a(j+1)-x2a(j)) 
     &      )

          bp(i,j,k) = bp(i,j,k) 
     & +dt*(
     &  (- nflux1(mbps,i+1,j,k)
     &   + nflux1(mbps,i  ,j,k))/(x1a(i+1)-x1a(i)) 
     & +(- nflux2(mbps,i,j+1,k)
     &   + nflux2(mbps,i,j  ,k))/(x2a(j+1)-x2a(j)) 
     &      )

!          write(6,*) i,j,k,bp(i,j,k)
      enddo
      enddo
      enddo


      return
      end subroutine UpdateConsv

      subroutine EvaulateCh
      use commons
      use fluxmod
      implicit none
      integer :: i,j,k,n
      real(8),parameter:: drate=0.1d0 ! 
! local variable
      real(8):: dh1l,dh2l,dh3l,dhl,dhd
      real(8):: ch1l,ch2l,ch3l,chl,chd
      real(8):: cts,css,cms
      real(8),parameter:: huge=1.0d90 

      chd = 0.0d0
      ch1l = 0.0d0; ch2l = 0.0d0; ch3l = 0.0d0
      dhd = huge
      dh1l =  huge; dh2l =  huge; dh3l =  huge

      do k=ks,ke
      do j=js,je
      do i=is,ie
            css  = (gam*svc(npre,i,j,k)/ svc(nden,i,j,k))
            cts  = css ! cs^2+c_a^2
     &          + (svc(nbm1,i,j,k)**2+svc(nbm2,i,j,k)**2+svc(nbm3,i,j,k)**2)/svc(nden,i,j,k)
            cms  = sqrt((cts +sqrt(cts**2
     &      -4.0d0*css*svc(nbm1,i,j,k)**2/svc(nden,i,j,k)))/2.0d0)
            ch1l = ( abs(svc(nve1,i,j,k)) + cms )
            dh1l =  (x1a(i+1)-x1a(i))

            cms  = sqrt((cts +sqrt(cts**2
     &      -4.0d0*css*svc(nbm2,i,j,k)**2/svc(nden,i,j,k)))/2.0d0)
            ch2l = ( abs(svc(nve2,i,j,k)) + cms )
            dh2l = (x2a(j+1)-x2a(j)) 

!         if (ldimen .eq. 3) then
!            cms  = sqrt((cts +sqrt(cts**2
!     &      -4.0d0*css*svc(nbm3,i,j,k)**2/svc(nden,i,j,k)))/2.0d0)
!            ch3l = ( abs(svc(nve3,i,j,k)) + cms )
!            dh3l = (x3a(k+1)-x3a(k))
!         endif

         chl     = max(ch1l,ch2l,ch3l)
         dhl     = min(dh1l,dh2l,dh3l)
         chd     = max(chl,chd)
         dhd     = min(dhl,dhd)
      enddo
      enddo
      enddo

      chg      = chd

      return
      end subroutine  EvaulateCh

      subroutine DampPsi
      use commons
      use fluxmod
      implicit none
      integer :: i,j,k,n
      real(8),parameter:: alphabp=0.1d0 !
      real(8):: taui
      real(8):: dhl,dh1l,dh2l,dh3l
      real(8),parameter:: huge=1.0d90 

      dh1l=huge
      dh2l=huge
      dh3l=huge

      do k=ks,ke
      do j=js,je
      do i=is,ie
            dh1l = x1a(i+1)-x1a(i)
            dh2l = x2a(j+1)-x2a(j)
!         if(ldimen .eq. 3) then
!            dh3l = dx3a(k)
!         endif

         dhl = min(dh1l,dh2l,dh3l)
         taui = alphabp * chg /dhl ! cm/s /cm => 1/s
         bp(i,j,k) = bp(i,j,k)*(1.0d0 - dt*taui) ! if dt = dtloc, damping by factor of (1.0-drate)
      enddo
      enddo
      enddo

      return
      end subroutine  DampPsi

      subroutine GravForce
      use commons
      use fluxmod
      implicit none
      integer :: i,j,k,n

      do k=ks,ke
      do j=js,je
      do i=is,ie+1
         gp1a(i  ,j,k) = gp(i,j,k)
     & - 0.5d0*(gp(i  ,j,k)-gp(i-1,j,k))

         gp1a(i+1,j,k) = gp(i,j,k)
     & + 0.5d0*(gp(i+1,j,k)-gp(i  ,j,k))

       grvsrc1(i,j,k) = (gp1a(i+1,j,k)-gp1a(i,j,k))/(x1a(i+1)-x1a(i))*d(i,j,k)

      enddo
      enddo
      enddo

      do k=ks,ke
      do i=is,ie
      do j=js,je+1
         gp2a(i  ,j,k) = gp(i,j,k)
     & - 0.5d0*(gp(i  ,j,k)-gp(i,j-1,k))

         gp2a(i,j+1,k) = gp(i,j,k)
     & + 0.5d0*(gp(i,j+1,k)-gp(i  ,j,k))

       grvsrc2(i,j,k) = (gp2a(i,j+1,k)-gp2a(i,j,k))/(x2a(j+1)-x2a(j))*d(i,j,k)

      enddo
      enddo
      enddo

       grvsrc3(:,:,:) = 0.0d0

      return
      end subroutine  GravForce

      subroutine Output(flag)
      use commons
      implicit none
      logical,intent(in)::flag
      integer::i,j,k
      character(20),parameter::dirname="xydata/"
      character(40)::filename
      real(8),save::tout
      data tout / 0.0d0 / 
      real(8),parameter:: dtout=1.0d-1
      integer::nout
      data nout / 1 /
      integer,parameter::unitout=13

      if(time .lt. tout+dtout .and. .not. flag ) return

      write(filename,'(a2,i5.5,a4)')"Sc",nout,".xss"

      filename = trim(dirname)//filename
      open(unitout,file=filename,status='replace',form='formatted') 

      write(unitout,'(1a,1(1x,F5.2))') "# ",time
      k=ks
!      do j=js,je
!      do i=is,ie
      do j=1,jn-1
      do i=1,in-1
         write(unitout,'(12(1x,E15.6e3))') x1b(i),x2b(j)              !  2 
     &                                 , d(i,j,k), p(i,j,k),gp(i,j,k) !  5
     &                                 ,v1(i,j,k),v2(i,j,k),v3(i,j,k) !  8
     &                                 ,b1(i,j,k),b2(i,j,k),b3(i,j,k) ! 11
     &                                 ,bp(i,j,k)                     ! 12
      enddo
         write(unitout,*)
      enddo
      close(unitout)

      write(6,*) "gp",gp(is,js,ks)
      write(6,*) "output:",nout,time

      nout=nout+1
      tout=time

      return
      end subroutine Output
