      module commons
      implicit none
      integer::nhy
      integer,parameter::nhymax=2000000
      real(8)::time,dt
      data time / 0.0d0 /
      real(8),parameter:: timemax=1.0d1
      real(8),parameter:: dtout=1.0d-1

      integer,parameter::ngrid=128
      integer,parameter::mgn=2
      integer,parameter::in=ngrid+2*mgn+1 &
     &                  ,jn=ngrid+2*mgn+1 &
     &                  ,kn=1
      integer,parameter::is=mgn+1 &
     &                  ,js=mgn+1 &
     &                  ,ks=1 
      integer,parameter::ie=ngrid+mgn &
     &                  ,je=ngrid+mgn &
     &                  ,ke=1

      real(8),parameter:: x1min=0.0d0,x1max=1.0d0
      real(8),parameter:: x2min=0.0d0,x2max=1.0d0
      real(8),dimension(in)::x1a,x1b
      real(8),dimension(jn)::x2a,x2b
      real(8),dimension(kn)::x3a,x3b

      real(8),dimension(in,jn,kn)::d,et,mv1,mv2,mv3
      real(8),dimension(in,jn,kn)::p,ei,v1,v2,v3,cs
      real(8),dimension(in,jn,kn)::b1,b2,b3,bp

      end module commons

      module eosmod
      implicit none
! adiabatic
      real(8),parameter::gam=5.0d0/3.0d0 !! adiabatic index
! isothermal
!      real(8)::csiso  !! isothemal sound speed
end module eosmod
      
      module fluxmod
      use commons, only : in,jn,kn
      implicit none
      real(8):: chg
      integer,parameter::nden=1,nve1=2,nve2=3,nve3=4,nene=5,npre=6,ncsp=7 &
     &                         ,nbm1=8,nbm2=9,nbm3=10,nbps=11
      integer,parameter::nhyd=11
      real(8),dimension(nhyd,in,jn,kn):: svc

      integer,parameter::mudn= 1,muvu= 2,muvv= 3,muvw= 4,muet= 5 &
     &                          ,mubu= 6,mubv= 7,mubw= 8,mubp= 9 &
     &                  ,mfdn=10,mfvu=11,mfvv=12,mfvw=13,mfet=14 &
     &                          ,mfbu=15,mfbv=16,mfbw=17,mfbp=18 &
     &                          ,mcsp=19,mvel=20,mpre=21
      integer,parameter:: mflx=9,madd=3

      integer,parameter:: mden=1,mrv1=2,mrv2=3,mrv3=4,meto=5   &
     &                          ,mrvu=muvu,mrvv=muvv,mrvw=muvw &
     &                          ,mbm1=6,mbm2=7,mbm3=8,mbps=9   &
     &                          ,mbmu=mubu,mbmv=mubv,mbmw=mubw
      real(8),dimension(mflx,in,jn,kn):: nflux1,nflux2,nflux3

      end module fluxmod

      program main
      use commons
      implicit none
      logical::is_final
      integer :: i,j,k
      data is_final /.false./
      real(8) :: divergenceB, error

      write(6,*) "setup grids and fiels"
      call GenerateGrid
      call GenerateProblem
      call ConsvVariable
      write(6,*) "entering main loop"
! main loop
                                  write(6,*)"step","time","dt"
      open(1,file="divB0.1_t.dat",action="write")
      mloop: do nhy=1,nhymax
         if(mod(nhy,100) .eq. 0 ) write(6,*)nhy,time,dt
         call TimestepControl
         call BoundaryCondition
         call StateVevtor
         call EvaulateCh
         call NumericalFlux1
         call NumericalFlux2
         call UpdateConsv
         call DampPsi
         call PrimVariable
         time=time+dt
         call Output(is_final)


         divergenceB = 0.d00
           do k=ks,ke
           do j=js,je
           do i=is,ie
              error = dabs( (b1(i+1,j,k) - b1(i-1,j,k) )/(x1a(i+1)-x1a(i-1))  &
                         + ( b2(i,j+1,k) - b2(i,j-1,k) )/(x2a(j+1)-x2a(j-1)) )
     !         error = dabs( ( Q(i+1,j,k,IB1) - Q(i-1,j,k,IB1) )/(x1a(i+1)-x1a(i-1))  &
     !                    + ( Q(i,j+1,k,IB2) - Q(i,j-1,k,IB2) )/(x2a(j+1)-x2a(j-1)) )/ &
     !                    dsqrt( Q(i,j,k,IB1)**2 + Q(i,j,k,IB2)**2 )*min(x1a(i+1) - x1a(i),x2a(j+1)-x2a(j))
     
              divergenceB = divergenceB + error
           enddo
           enddo
           enddo
           divergenceB = divergenceB/dble((ie-is+1)*(je-js+1)*(ke-ks+1))
           write(1,*) time, divergenceB


         if(time > timemax) exit mloop
      enddo mloop
      close(1)

      is_final = .true.
      call Output(is_final)

      write(6,*) "program has been finished"
      end program main

      subroutine GenerateGrid
      use commons
      use eosmod
      implicit none
      real(8)::dx,dy
      integer::i,j,k
! x coordinates
      dx=(x1max-x1min)/dble(ngrid)
      do i=1,in
         x1a(i) = dx*(i-(mgn+1))+x1min
      enddo
      do i=1,in-1
         x1b(i) = 0.5d0*(x1a(i+1)+x1a(i))
      enddo
 
! y coordinates
      dy=(x2max-x2min)/dble(ngrid)
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
      use eosmod
      implicit none
      integer::i,j,k
      real(8):: B0
      real(8):: pi,r0,rad
      pi=acos(-1.0d0)

      B0 = 1.0d0/gam
      r0 = 1.0d0/sqrt(8.0d0)

      d(:,:,:) = 1.0d0

      do k=ks,ke
      do j=js,je
      do i=is,ie
          rad = sqrt(x1b(i)**2 + x2b(j)**2)
          d(i,j,k) = 1.0d0
          p(i,j,k) = 6.0d0
         v1(i,j,k) = 1.0d0
         v2(i,j,k) = 1.0d0
         v3(i,j,k) = 0.0d0
           if( rad < r0 ) then
         b1(i,j,k) = ( (rad/r0)**8 - 2.0d0*(rad/r0)**4 + 1.0d0)/sqrt(4.0d0*pi)
     else
         b1(i,j,k) = 0.0d0
     endif
         b2(i,j,k) = 0.0d0
         b3(i,j,k) = 1.0d0/sqrt(4.0d0*pi)
         bp(i,j,k) = 0.0d0
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

      k=ks
      do j=1,jn-1
      do i=1,mgn
           d(i,j,k) =  d(ie-mgn+i,j,k)
          ei(i,j,k) = ei(ie-mgn+i,j,k)
          v1(i,j,k) = v1(ie-mgn+i,j,k)
          v2(i,j,k) = v2(ie-mgn+i,j,k)
          v3(i,j,k) = v3(ie-mgn+i,j,k)
          b1(i,j,k) = b1(ie-mgn+i,j,k)
          b2(i,j,k) = b2(ie-mgn+i,j,k)
          b3(i,j,k) = b3(ie-mgn+i,j,k)
          bp(i,j,k) = bp(ie-mgn+i,j,k)
      enddo
      enddo

      k=ks
      do j=1,jn-1
      do i=1,mgn
           d(ie+i,j,k) =  d(is+i-1,j,k)
          ei(ie+i,j,k) = ei(is+i-1,j,k)
          v1(ie+i,j,k) = v1(is+i-1,j,k)
          v2(ie+i,j,k) = v2(is+i-1,j,k)
          v3(ie+i,j,k) = v3(is+i-1,j,k)
          b1(ie+i,j,k) = b1(is+i-1,j,k)
          b2(ie+i,j,k) = b2(is+i-1,j,k)
          b3(ie+i,j,k) = b3(is+i-1,j,k)
          bp(ie+i,j,k) = bp(is+i-1,j,k)
      enddo
      enddo

      k=ks
      do i=1,in-1
      do j=1,mgn
           d(i,j,k) =  d(i,je-mgn+j,k)
          ei(i,j,k) = ei(i,je-mgn+j,k)
          v1(i,j,k) = v1(i,je-mgn+j,k)
          v2(i,j,k) = v2(i,je-mgn+j,k)
          v3(i,j,k) = v3(i,je-mgn+j,k)
          b1(i,j,k) = b1(i,je-mgn+j,k)
          b2(i,j,k) = b2(i,je-mgn+j,k)
          b3(i,j,k) = b3(i,je-mgn+j,k)
          bp(i,j,k) = bp(i,je-mgn+j,k)
      enddo
      enddo

      k=ks
      do i=1,in-1
      do j=1,mgn
           d(i,je+j,k) =  d(i,js+j-1,k)
          ei(i,je+j,k) = ei(i,js+j-1,k)
          v1(i,je+j,k) = v1(i,js+j-1,k)
          v2(i,je+j,k) = v2(i,js+j-1,k)
          v3(i,je+j,k) = v3(i,js+j-1,k)
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
          et(i,j,k) = 0.5d0*d(i,j,k)*(   &
     &                    +v1(i,j,k)**2  &
     &                    +v2(i,j,k)**2  &
     &                    +v3(i,j,k)**2) &
     &               +0.5d0*(            &
     &                    +b1(i,j,k)**2  &
     &                    +b2(i,j,k)**2  &
     &                    +b3(i,j,k)**2) &
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
      use eosmod
      implicit none
      integer::i,j,k
      do k=ks,ke
      do j=js,je
      do i=is,ie
          v1(i,j,k) = mv1(i,j,k)/d(i,j,k)
          v2(i,j,k) = mv2(i,j,k)/d(i,j,k)
          v3(i,j,k) = mv3(i,j,k)/d(i,j,k)

          ei(i,j,k) =  et(i,j,k)          &
     &          -0.5d0*d(i,j,k)*(         &
     &                    +v1(i,j,k)**2   &
     &                    +v2(i,j,k)**2   &
     &                    +v3(i,j,k)**2)  &
     &          -0.5d0*(                  &
     &                    +b1(i,j,k)**2   &
     &                    +b2(i,j,k)**2   &
     &                    +b3(i,j,k)**2)

! adiabatic
           p(i,j,k) =  ei(i,j,k)*(gam-1.0d0)
          cs(i,j,k) =  sqrt(gam*p(i,j,k)/d(i,j,k))
! isotermal
!           p(i,j,k) =  d(i,j,k)*csiso**2
!          cs(i,j,k) =  csiso
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
         ctot = sqrt(cs(i,j,k)**2 &
     &            +( b1(i,j,k)**2 &
     &              +b2(i,j,k)**2 &
     &              +b3(i,j,k)**2 &
     &                           )/d(i,j,k))
         dtl1 =(x1a(i+1)-x1a(i))/(abs(v1(i,j,k)) + ctot)
         dtl2 =(x2a(j+1)-x2a(j))/(abs(v2(i,j,k)) + ctot)
!         dtl3 =(x3a(k+1)-x3a(k))/(abs(v3(i,j,k)) + ctot)
         dtlocal = min (dtl1,dtl2)
         if(dtlocal .lt. dtmin) dtmin = dtlocal
      enddo
      enddo
      enddo

      dt = 0.05d0 * dtmin

      return
      end subroutine TimestepControl

      subroutine StateVevtor
      use commons
      use fluxmod
      use eosmod
      implicit none
      integer::i,j,k

      k=ks
      do j=1,jn-1
      do i=1,in-1
         svc(nden,i,j,k) =  d(i,j,k)
         svc(nve1,i,j,k) = v1(i,j,k)
         svc(nve2,i,j,k) = v2(i,j,k)
         svc(nve3,i,j,k) = v3(i,j,k)
         svc(nbm1,i,j,k) = b1(i,j,k)
         svc(nbm2,i,j,k) = b2(i,j,k)
         svc(nbm3,i,j,k) = b3(i,j,k)
         svc(nbps,i,j,k) = bp(i,j,k)
! adiabatic
         svc(nene,i,j,k) = ei(i,j,k)/d(i,j,k)
         svc(npre,i,j,k) = ei(i,j,k)*(gam-1.0d0)
         svc(ncsp,i,j,k) = sqrt(gam*(gam-1.0d0)*ei(i,j,k)/d(i,j,k))
! isotermal
!         svc(nene,i,j,k) = csiso**2
!         svc(npre,i,j,k) = d(i,j,k)*csiso**2
!         svc(ncsp,i,j,k) = csiso
         p(i,j,k) = svc(npre,i,j,k)  ! for output boundary 
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
         d(n) = sign(1.0d0,a(n))*max(0.0d0,min(abs(a(n)) &
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
         d(n) = sign(1.0d0,a(n))*max(0.0d0,min(abs(a(n))          &
     &                                  ,sign(1.0d0,a(n))*b(n)    &
     &                                  ,sign(1.0d0,a(n))*c(n)))
      enddo

      return
      end subroutine MClimiter

      subroutine NumericalFlux1
      use commons, only: is,ie,in,js,je,jn,ks,ke,kn
      use fluxmod
      implicit none
      integer::i,j,k
      real(8),dimension(nhyd):: dsvp,dsvm,dsvc,dsv
      real(8),dimension(nhyd,in,jn,kn):: leftpr,rigtpr
      real(8),dimension(2*mflx+madd,in,jn,kn):: leftco,rigtco
      real(8),dimension(2*mflx+madd):: leftst,rigtst
      real(8),dimension(mflx):: nflux
      real(8):: ptl,css,cts 

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

      do j=js,je
      do i=is,ie+1
!====================
! Left
!====================

! Consvative variables
         leftco(mudn,i,j,k)=leftpr(nden,i,j,k) ! rho
         leftco(muvu,i,j,k)=leftpr(nve1,i,j,k)*leftpr(nden,i,j,k)   ! rho v_x
         leftco(muvv,i,j,k)=leftpr(nve2,i,j,k)*leftpr(nden,i,j,k)   ! rho v_y
         leftco(muvw,i,j,k)=leftpr(nve3,i,j,k)*leftpr(nden,i,j,k)   ! rho v_z
         leftco(muet,i,j,k)=leftpr(nene,i,j,k)*leftpr(nden,i,j,k) & ! e_i
     &               +0.5d0*leftpr(nden,i,j,k)*(                  &
     &                     +leftpr(nve1,i,j,k)**2                 &
     &                     +leftpr(nve2,i,j,k)**2                 &
     &                     +leftpr(nve3,i,j,k)**2)                & ! + rho v^2/2
     &               +0.5d0*                    (                 &
     &                     +leftpr(nbm1,i,j,k)**2                 &
     &                     +leftpr(nbm2,i,j,k)**2                 &
     &                     +leftpr(nbm3,i,j,k)**2)                  ! + B^2/2

         leftco(mubu,i,j,k)=leftpr(nbm1,i,j,k)  ! b_x
         leftco(mubv,i,j,k)=leftpr(nbm2,i,j,k)  ! b_y
         leftco(mubw,i,j,k)=leftpr(nbm3,i,j,k)  ! b_z
         leftco(mubp,i,j,k)=leftpr(nbps,i,j,k)  ! psi

! Flux
         ptl = leftpr(npre,i,j,k) + ( leftpr(nbm1,i,j,k)**2        &
     &                               +leftpr(nbm2,i,j,k)**2        &
     &                               +leftpr(nbm3,i,j,k)**2)/2.0d0 

         leftco(mfdn,i,j,k)=leftpr(nden,i,j,k)                   *leftpr(nve1,i,j,k)
         leftco(mfvu,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve1,i,j,k)*leftpr(nve1,i,j,k) &
     &                     +ptl-leftpr(nbm1,i,j,k)**2
         leftco(mfvv,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve2,i,j,k)*leftpr(nve1,i,j,k) &
     &                                        -leftpr(nbm2,i,j,k)*leftpr(nbm1,i,j,k)
         leftco(mfvw,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve3,i,j,k)*leftpr(nve1,i,j,k) &
     &                                        -leftpr(nbm3,i,j,k)*leftpr(nbm1,i,j,k)
         leftco(mfet,i,j,k)=(leftco(muet,i,j,k)+ptl)*leftpr(nve1,i,j,k) &
     &                     -( leftpr(nbm1,i,j,k)*leftpr(nve1,i,j,k)     &
     &                       +leftpr(nbm2,i,j,k)*leftpr(nve2,i,j,k)     &
     &                       +leftpr(nbm3,i,j,k)*leftpr(nve3,i,j,k))*leftpr(nbm1,i,j,k)

         leftco(mfbu,i,j,k) =  0.0d0
         leftco(mfbv,i,j,k) =  leftpr(nbm2,i,j,k)*leftpr(nve1,i,j,k) &
     &                        -leftpr(nve2,i,j,k)*leftpr(nbm1,i,j,k)
         leftco(mfbw,i,j,k) =  leftpr(nbm3,i,j,k)*leftpr(nve1,i,j,k) &
     &                        -leftpr(nve3,i,j,k)*leftpr(nbm1,i,j,k)
         leftco(mfbp,i,j,k) = 0.0d0  ! psi
     
         css =leftpr(ncsp,i,j,k)**2
         cts =  css  & !c_s^2*c_a^2
     &                       +( leftpr(nbm1,i,j,k)**2  &
     &                         +leftpr(nbm2,i,j,k)**2  &
     &                         +leftpr(nbm3,i,j,k)**2)/leftpr(nden,i,j,k) 

         leftco(mcsp,i,j,k)= sqrt((cts +sqrt(cts**2                 &
     &                             -4.0d0*css*leftpr(nbm1,i,j,k)**2 &
     &                                       /leftpr(nden,i,j,k))   &
     &                            )/2.0d0)
         leftco(mvel,i,j,k)= leftpr(nve1,i,j,k)
         leftco(mpre,i,j,k)= ptl
!====================
! Right
!====================
! Consvative variables
         rigtco(mudn,i,j,k)=rigtpr(nden,i,j,k) ! rho
         rigtco(muvu,i,j,k)=rigtpr(nve1,i,j,k)*rigtpr(nden,i,j,k)   ! rho v_x
         rigtco(muvv,i,j,k)=rigtpr(nve2,i,j,k)*rigtpr(nden,i,j,k)   ! rho v_y
         rigtco(muvw,i,j,k)=rigtpr(nve3,i,j,k)*rigtpr(nden,i,j,k)   ! rho v_z
         rigtco(muet,i,j,k)=rigtpr(nene,i,j,k)*rigtpr(nden,i,j,k) & ! e_i
     &               +0.5d0*rigtpr(nden,i,j,k)*(                  &
     &                     +rigtpr(nve1,i,j,k)**2                 &
     &                     +rigtpr(nve2,i,j,k)**2                 &
     &                     +rigtpr(nve3,i,j,k)**2)                & ! + rho v^2/2
     &               +0.5d0*                    (                 &
     &                     +rigtpr(nbm1,i,j,k)**2                 &
     &                     +rigtpr(nbm2,i,j,k)**2                 &
     &                     +rigtpr(nbm3,i,j,k)**2)                 ! + B^2/2

         rigtco(mubu,i,j,k)=rigtpr(nbm1,i,j,k)  ! b_x
         rigtco(mubv,i,j,k)=rigtpr(nbm2,i,j,k)  ! b_y
         rigtco(mubw,i,j,k)=rigtpr(nbm3,i,j,k)  ! b_z
         rigtco(mubp,i,j,k)=rigtpr(nbps,i,j,k)  ! psi

! Flux
         ptl = rigtpr(npre,i,j,k) + ( rigtpr(nbm1,i,j,k)**2        &
     &                               +rigtpr(nbm2,i,j,k)**2        &
     &                               +rigtpr(nbm3,i,j,k)**2)/2.0d0

         rigtco(mfdn,i,j,k)=rigtpr(nden,i,j,k)                   *rigtpr(nve1,i,j,k)
         rigtco(mfvu,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve1,i,j,k)*rigtpr(nve1,i,j,k) &
     &                     +ptl-rigtpr(nbm1,i,j,k)**2
         rigtco(mfvv,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve2,i,j,k)*rigtpr(nve1,i,j,k) &
     &                                        -rigtpr(nbm2,i,j,k)*rigtpr(nbm1,i,j,k)
         rigtco(mfvw,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve3,i,j,k)*rigtpr(nve1,i,j,k) &
     &                                        -rigtpr(nbm3,i,j,k)*rigtpr(nbm1,i,j,k)
         rigtco(mfet,i,j,k)=(rigtco(muet,i,j,k)+ptl)*rigtpr(nve1,i,j,k) &
     &                     -( rigtpr(nbm1,i,j,k)*rigtpr(nve1,i,j,k)     &
     &                       +rigtpr(nbm2,i,j,k)*rigtpr(nve2,i,j,k)     &
     &                       +rigtpr(nbm3,i,j,k)*rigtpr(nve3,i,j,k))*rigtpr(nbm1,i,j,k)
     
         rigtco(mfbu,i,j,k) =  0.0d0
         rigtco(mfbv,i,j,k) =  rigtpr(nbm2,i,j,k)*rigtpr(nve1,i,j,k) &
     &                        -rigtpr(nve2,i,j,k)*rigtpr(nbm1,i,j,k)
         rigtco(mfbw,i,j,k) =  rigtpr(nbm3,i,j,k)*rigtpr(nve1,i,j,k) &
     &                        -rigtpr(nve3,i,j,k)*rigtpr(nbm1,i,j,k)
         rigtco(mfbp,i,j,k) = 0.0d0  ! b_z
         css =rigtpr(ncsp,i,j,k)**2
         cts =  css   &!c_s^2*c_a^2
     &                       +( rigtpr(nbm1,i,j,k)**2 &
     &                         +rigtpr(nbm2,i,j,k)**2 &
     &                         +rigtpr(nbm3,i,j,k)**2)/rigtpr(nden,i,j,k) 


         rigtco(mcsp,i,j,k)= sqrt((cts +sqrt(cts**2                 &
     &                             -4.0d0*css*rigtpr(nbm1,i,j,k)**2 &
     &                                       /rigtpr(nden,i,j,k))   &
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
!         call HLLC(leftst,rigtst,nflux)
         call HLLD(leftst,rigtst,nflux)
         nflux1(mden,i,j,k)=nflux(mden)
         nflux1(mrv1,i,j,k)=nflux(mrvu)
         nflux1(mrv2,i,j,k)=nflux(mrvv)
         nflux1(mrv3,i,j,k)=nflux(mrvw)
         nflux1(meto,i,j,k)=nflux(meto)
         nflux1(mbm1,i,j,k)=nflux(mbmu)
         nflux1(mbm2,i,j,k)=nflux(mbmv)
         nflux1(mbm3,i,j,k)=nflux(mbmw)

         nflux1(mbm1,i,j,k) =  0.5d0*(leftst(mubp)+rigtst(mubp)) &
     &                    -0.5d0*chg*(rigtst(mubu)-leftst(mubu))        ! finite volume
         nflux1(mbps,i,j,k) = (0.5d0*(leftst(mubu)+rigtst(mubu)) &
     &                    -0.5d0/chg*(rigtst(mubp)-leftst(mubp)))*chg**2 ! finite volume

!         write(6,*) "bpf1",nflux1(mbps,i,j,k)

      enddo
      enddo

      return
      end subroutine Numericalflux1

      subroutine NumericalFlux2
      use commons, only: is,ie,in,js,je,jn,ks,ke,kn
      use fluxmod
      implicit none
      integer::i,j,k
      real(8),dimension(nhyd):: dsvp,dsvm,dsvc,dsv
      real(8),dimension(nhyd,in,jn,kn):: leftpr,rigtpr
      real(8),dimension(2*mflx+madd,in,jn,kn):: leftco,rigtco
      real(8),dimension(2*mflx+madd):: leftst,rigtst
      real(8),dimension(mflx):: nflux
      real(8):: ptl,css,cts 

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

      do i=is,ie
      do j=js,je+1
         leftco(mudn,i,j,k)=leftpr(nden,i,j,k)
         leftco(muvw,i,j,k)=leftpr(nve1,i,j,k)*leftpr(nden,i,j,k)
         leftco(muvu,i,j,k)=leftpr(nve2,i,j,k)*leftpr(nden,i,j,k)   ! rho v
         leftco(muvv,i,j,k)=leftpr(nve3,i,j,k)*leftpr(nden,i,j,k)
         leftco(muet,i,j,k)=leftpr(nene,i,j,k)*leftpr(nden,i,j,k) & ! internal
     &               +0.5d0*leftpr(nden,i,j,k)*(                  &
     &                     +leftpr(nve1,i,j,k)**2                 &
     &                     +leftpr(nve2,i,j,k)**2                 &
     &                     +leftpr(nve3,i,j,k)**2)                & ! kinetic
     &               +0.5d0*                    (                 &
     &                     +leftpr(nbm1,i,j,k)**2                 &
     &                     +leftpr(nbm2,i,j,k)**2                 &
     &                     +leftpr(nbm3,i,j,k)**2) ! magnetic

         leftco(mubw,i,j,k)=leftpr(nbm1,i,j,k)  ! b_x
         leftco(mubu,i,j,k)=leftpr(nbm2,i,j,k)  ! b_y
         leftco(mubv,i,j,k)=leftpr(nbm3,i,j,k)  ! b_z
         leftco(mubp,i,j,k)=leftpr(nbps,i,j,k)  ! psi

         ptl = leftpr(npre,i,j,k) + ( leftpr(nbm1,i,j,k)**2        &
     &                               +leftpr(nbm2,i,j,k)**2        &
     &                               +leftpr(nbm3,i,j,k)**2)/2.0d0 

         leftco(mfdn,i,j,k)=leftpr(nden,i,j,k)                   *leftpr(nve2,i,j,k)   ! rho v
         leftco(mfvw,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve1,i,j,k)*leftpr(nve2,i,j,k) &
     &                                        -leftpr(nbm1,i,j,k)*leftpr(nbm2,i,j,k)
         leftco(mfvu,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve2,i,j,k)*leftpr(nve2,i,j,k) &
     &                     +ptl-leftpr(nbm2,i,j,k)**2
         leftco(mfvv,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve3,i,j,k)*leftpr(nve2,i,j,k) &
     &                                        -leftpr(nbm3,i,j,k)*leftpr(nbm2,i,j,k)
         leftco(mfet,i,j,k)=(leftco(muet,i,j,k)+ptl)*leftpr(nve2,i,j,k) &
     &                     -( leftpr(nbm1,i,j,k)*leftpr(nve1,i,j,k)     &
     &                       +leftpr(nbm2,i,j,k)*leftpr(nve2,i,j,k)     &
     &                       +leftpr(nbm3,i,j,k)*leftpr(nve3,i,j,k))*leftpr(nbm2,i,j,k)

         leftco(mfbw,i,j,k) =  leftpr(nbm1,i,j,k)*leftpr(nve2,i,j,k) &
     &                        -leftpr(nve1,i,j,k)*leftpr(nbm2,i,j,k)
         leftco(mfbu,i,j,k) =  0.0d0
         leftco(mfbv,i,j,k) =  leftpr(nbm3,i,j,k)*leftpr(nve2,i,j,k) &
     &                        -leftpr(nve3,i,j,k)*leftpr(nbm2,i,j,k)
         leftco(mfbp,i,j,k) = 0.0d0  ! psi
     
         css = leftpr(ncsp,i,j,k)**2
         cts =  css  & !c_s^2*c_a^2
     &                       +( leftpr(nbm1,i,j,k)**2  &
     &                         +leftpr(nbm2,i,j,k)**2  &
     &                         +leftpr(nbm3,i,j,k)**2)/leftpr(nden,i,j,k) 

         leftco(mcsp,i,j,k)= sqrt((cts +sqrt(cts**2                  &
     &                             -4.0d0*css*leftpr(nbm2,i,j,k)**2  &
     &                                          /leftpr(nden,i,j,k)) &
     &                            )/2.0d0)
         leftco(mvel,i,j,k)= leftpr(nve2,i,j,k)
         leftco(mpre,i,j,k)= ptl


         rigtco(mudn,i,j,k)=rigtpr(nden,i,j,k)
         rigtco(muvw,i,j,k)=rigtpr(nve1,i,j,k)*rigtpr(nden,i,j,k)
         rigtco(muvu,i,j,k)=rigtpr(nve2,i,j,k)*rigtpr(nden,i,j,k)   ! rho v
         rigtco(muvv,i,j,k)=rigtpr(nve3,i,j,k)*rigtpr(nden,i,j,k)
         rigtco(muet,i,j,k)=rigtpr(nene,i,j,k)*rigtpr(nden,i,j,k) & ! internal
     &               +0.5d0*rigtpr(nden,i,j,k)*(                  &
     &                     +rigtpr(nve1,i,j,k)**2                 &
     &                     +rigtpr(nve2,i,j,k)**2                 &
     &                     +rigtpr(nve3,i,j,k)**2)                & ! kinetic
     &               +0.5d0*                    (                 &
     &                     +rigtpr(nbm1,i,j,k)**2                 &
     &                     +rigtpr(nbm2,i,j,k)**2                 &
     &                     +rigtpr(nbm3,i,j,k)**2) ! magnetic

         rigtco(mubw,i,j,k)=rigtpr(nbm1,i,j,k)  ! b_x
         rigtco(mubu,i,j,k)=rigtpr(nbm2,i,j,k)  ! b_y
         rigtco(mubv,i,j,k)=rigtpr(nbm3,i,j,k)  ! b_z
         rigtco(mubp,i,j,k)=rigtpr(nbps,i,j,k)  ! psi

         ptl = rigtpr(npre,i,j,k) + ( rigtpr(nbm1,i,j,k)**2 &
     &                               +rigtpr(nbm2,i,j,k)**2 &
     &                               +rigtpr(nbm3,i,j,k)**2)/2.0d0 

         rigtco(mfdn,i,j,k)=rigtpr(nden,i,j,k)                   *rigtpr(nve2,i,j,k) ! rho v
         rigtco(mfvw,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve1,i,j,k)*rigtpr(nve2,i,j,k) &
     &                     -rigtpr(nbm1,i,j,k)*rigtpr(nbm2,i,j,k)
         rigtco(mfvu,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve2,i,j,k)*rigtpr(nve2,i,j,k) &
     &                     +ptl-rigtpr(nbm2,i,j,k)**2
         rigtco(mfvv,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve3,i,j,k)*rigtpr(nve2,i,j,k) &
     &                     -rigtpr(nbm3,i,j,k)*rigtpr(nbm2,i,j,k)
         rigtco(mfet,i,j,k)=(rigtco(muet,i,j,k)+ptl)*rigtpr(nve2,i,j,k) &
     &                     -( rigtpr(nbm1,i,j,k)*rigtpr(nve1,i,j,k)     &
     &                       +rigtpr(nbm2,i,j,k)*rigtpr(nve2,i,j,k)     &
     &                       +rigtpr(nbm3,i,j,k)*rigtpr(nve3,i,j,k))*rigtpr(nbm2,i,j,k)
         rigtco(mfbw,i,j,k) =  rigtpr(nbm1,i,j,k)*rigtpr(nve2,i,j,k) &
     &                        -rigtpr(nve1,i,j,k)*rigtpr(nbm2,i,j,k)
         rigtco(mfbu,i,j,k) =  0.0d0
         rigtco(mfbv,i,j,k) =  rigtpr(nbm3,i,j,k)*rigtpr(nve2,i,j,k) &
     &                        -rigtpr(nve3,i,j,k)*rigtpr(nbm2,i,j,k)
         rigtco(mfbp,i,j,k) = 0.0d0  ! psi
     
         css =rigtpr(nene,i,j,k)**2
         cts =  css  & !c_s^2*c_a^2
     &                       +( rigtpr(nbm1,i,j,k)**2 &
     &                         +rigtpr(nbm2,i,j,k)**2 &
     &                         +rigtpr(nbm3,i,j,k)**2)/rigtpr(nden,i,j,k) 

         rigtco(mcsp,i,j,k)= sqrt((cts +sqrt(cts**2                  &
     &                             -4.0d0*css*rigtpr(nbm2,i,j,k)**2  &
     &                                          /rigtpr(nden,i,j,k)) &
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
!         call HLLC(leftst,rigtst,nflux)
         call HLLD(leftst,rigtst,nflux)

         nflux2(mden,i,j,k)=nflux(mden)
         nflux2(mrv1,i,j,k)=nflux(mrvw)
         nflux2(mrv2,i,j,k)=nflux(mrvu) ! mrv2=3, mrvu=2
         nflux2(mrv3,i,j,k)=nflux(mrvv)
         nflux2(meto,i,j,k)=nflux(meto)
         nflux2(mbm1,i,j,k)=nflux(mbmw)
         nflux2(mbm2,i,j,k)=nflux(mbmu)
         nflux2(mbm3,i,j,k)=nflux(mbmv)

         nflux2(mbm2,i,j,k) =  0.5d0*(leftst(mubp)+rigtst(mubp)) &
     &                    -0.5d0*chg*(rigtst(mubu)-leftst(mubu))        ! finite volume
         nflux2(mbps,i,j,k) = (0.5d0*(leftst(mubu)+rigtst(mubu)) &
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
      use fluxmod, only: mflx,madd                &
     &                 , mudn,muvu,muvv,muvw,muet &
     &                 , mfdn,mfvu,mfvv,mfvw,mfet &
     &                 , mcsp,mvel,mpre           &
     &                 , mden,mrvu,mrvv,mrvw,meto &
     &                 , mubu,mubv,mubw,mubp      &
     &                 , mfbu,mfbv,mfbw,mfbp      &
     &                 , mbmu,mbmv,mbmw

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
      real(8) ::     byl,bzl
      real(8) ::     byr,bzr
      real(8) :: ptst

!----- U* ----
! qqlst ::  left state
! qqrst :: right state
      real(8) :: rolst,vxlst,vylst,vzlst,eelst
      real(8) :: rorst,vxrst,vyrst,vzrst,eerst
      real(8) :: rxlst,rylst,rzlst
      real(8) :: rxrst,ryrst,rzrst
      real(8) ::       bylst,bzlst
      real(8) ::       byrst,bzrst

!----- flux ---
! fqql ::  left physical flux
! fqqr :: right physical flux
      real(8) :: frol,frxl,fryl,frzl,feel
      real(8) ::           fbyl,fbzl
      real(8) :: fror,frxr,fryr,frzr,feer
      real(8) ::           fbyr,fbzr

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
        byl = leftst(mubv)
        bzl = leftst(mubv)
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
        byr = rigtst(mubv)
        bzr = rigtst(mubv)
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
        fbyl = leftst(mfbv)
        fbzl = leftst(mfbw)    

! Right value
! Left value
        fror = rigtst(mfdn)
        feer = rigtst(mfet)
        frxr = rigtst(mfvu)
        fryr = rigtst(mfvv) 
        frzr = rigtst(mfvw)
        fbyr = rigtst(mfbv)
        fbzr = rigtst(mfbw)

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

        bylst = rolst/rol * byl
        bzlst = rolst/rol * bzl

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

        byrst = rolst/rol * byr
        bzrst = rolst/rol * bzr
           
        eerst = (sdr*eer - ptr*vxr  + ptst*sm  )*isdmr
              
!----- Step 6. ----------------------------------------------------------|
! compute flux
        sign1 = sign(1.0d0,sm)    ! 1 for sm>0, -1 for sm<0
        maxs1 =  max(0.0d0,sign1) ! 1 sm>0, 0 for sm<0
        mins1 = -min(0.0d0,sign1) ! 0 sm>0,-1 for sm<0

        msl   = min(sl  ,0.0d0)   ! 0 for sl > 0, sl for sl < 0
        msr   = max(sr  ,0.0d0)   ! S_R > 0

        nflux(mden) = (frol+msl*(rolst-rol))*maxs1 &
     &               +(fror+msr*(rorst-ror))*mins1
        nflux(meto) = (feel+msl*(eelst-eel))*maxs1 &
     &               +(feer+msr*(eerst-eer))*mins1
        nflux(mrvu) = (frxl+msl*(rxlst-rxl))*maxs1 &
     &               +(frxr+msr*(rxrst-rxr))*mins1
        nflux(mrvv) = (fryl+msl*(rylst-ryl))*maxs1 &
     &               +(fryr+msr*(ryrst-ryr))*mins1
        nflux(mrvw) = (frzl+msl*(rzlst-rzl))*maxs1 &
     &               +(frzr+msr*(rzrst-rzr))*mins1
        nflux(mbmu) = 0.0d0
        nflux(mbmv) = (fbyl+msl*(bylst-byl))*maxs1 &
     &               +(fbyr+msr*(byrst-byr))*mins1
        nflux(mbmw) = (fbzl+msl*(bzlst-bzl))*maxs1 &
     &               +(fbzr+msr*(bzrst-bzr))*mins1

      return
      end subroutine HLLC

      subroutine HLLD(leftst,rigtst,nflux)
!=====================================================================
!
! HLLD Scheme
!
! Purpose
! Calculation of Numerical Flux by HLLD method
!
! Reference
!
! Input
! Output
!=====================================================================
      use fluxmod, only: mflx,madd                &
     &                 , mudn,muvu,muvv,muvw,muet &
     &                 , mfdn,mfvu,mfvv,mfvw,mfet &
     &                 , mcsp,mvel,mpre           &
     &                 , mden,mrvu,mrvv,mrvw,meto &
     &                 , mubu,mubv,mubw,mubp      &
     &                 , mfbu,mfbv,mfbw,mfbp      &
     &                 , mbmu,mbmv,mbmw

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
      real(8) :: bxs,byl,bzl
      real(8) ::     byr,bzr
      real(8) :: ptst

!----- U* ----
! qqlst ::  left state
! qqrst :: right state
      real(8) :: rolst,vxlst,vylst,vzlst,eelst
      real(8) :: rorst,vxrst,vyrst,vzrst,eerst
      real(8) :: rxlst,rylst,rzlst
      real(8) :: rxrst,ryrst,rzrst
      real(8) ::       bylst,bzlst
      real(8) ::       byrst,bzrst

!----- U** ----
! qqlst ::  left state
! qqrst :: right state
      real(8) :: vyldst,vzldst,eeldst
      real(8) :: vyrdst,vzrdst,eerdst
      real(8) :: ryldst,rzldst
      real(8) :: ryrdst,rzrdst
      real(8) ::       byldst,bzldst
      real(8) ::       byrdst,bzrdst

!----- flux ---
! fqql ::  left physical flux
! fqqr :: right physical flux
      real(8) :: frol,frxl,fryl,frzl,feel
      real(8) ::           fbyl,fbzl
      real(8) :: fror,frxr,fryr,frzr,feer
      real(8) ::           fbyr,fbzr

!----- wave speed ---
! sl ::  left-going fastest signal velocity
! sr :: right-going fastest signal velocity
! sm :: contact discontinuity velocity
! slst ::  left-going alfven velocity
! srst :: right-going alfven velocity
      real(8) :: sm,sl,sr,slst,srst

! cfl :: left-state Fast wave velocity
! cfr :: right-sate Fast wave velocity
      real(8) :: cfl,cfr

!--------------------
! temporary variables
      real(8) :: sdl,sdr,sdml,sdmr,isdml,isdmr,rosdl,rosdr
      real(8) :: temp
  
! no if
      real(8) :: sign1,maxs1,mins1
      real(8) :: msl,msr,mslst,msrst,temp1,invsumro,sqrtror,sqrtrol,abbx
      real(8) :: bxsq,temp_fst,eps,itf,vdbstl,vdbstr,signbx

!----- Step 0. ----------------------------------------------------------|
      eps = 1d-30
!---- Left state
        
        rol = leftst(mudn)
        eel = leftst(muet)
        rxl = leftst(muvu)
        ryl = leftst(muvv)
        rzl = leftst(muvw)
        vxl = leftst(muvu)/leftst(mudn)
        vyl = leftst(muvv)/leftst(mudn)
        vzl = leftst(muvw)/leftst(mudn)
        byl = leftst(mubv)
        bzl = leftst(mubw)
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
        byr = rigtst(mubv)
        bzr = rigtst(mubw)
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
        fbyl = leftst(mfbv)
        fbzl = leftst(mfbw)

! Right value
        fror = rigtst(mfdn)
        feer = rigtst(mfet)
        frxr = rigtst(mfvu)
        fryr = rigtst(mfvv) 
        frzr = rigtst(mfvw)
        fbyr = rigtst(mfbv)
        fbzr = rigtst(mfbw)


!----- Step 4. ----------------------------------------------------------|
! compute middle and alfven wave
!
        sdl = min(-1d-20,sl - vxl)
        sdr = max(1d-20,sr - vxr)
        rosdl = rol*sdl
        rosdr = ror*sdr

        temp = 1.0d0/(rosdr - rosdl)
! Eq. 45
        sm = (rosdr*vxr - rosdl*vxl - ptr + ptl)*temp
           
        sdml = min(-1d-20,sl - sm); isdml = 1.0d0/sdml
        sdmr = max(1d-20,sr - sm); isdmr = 1.0d0/sdmr

!----- Step 5. ----------------------------------------------------------|
! compute intermediate states
!
! Eq. 49
        ptst = (rosdr*ptl-rosdl*ptr+rosdl*rosdr*(vxr-vxl))*temp
		
!----- Step 5A. ----------------------------------------------------------|
! compute Ul*
!
           bxs = 0.5d0*(leftst(mubu)+rigtst(mubu))
           bxsq = bxs*bxs
           temp_fst = rosdl*sdml - bxsq
           sign1 = sign(1.0d0,abs(temp_fst)-eps)

           maxs1 = max(0.0d0,sign1)
           mins1 = min(0.0d0,sign1)

           itf = 1.0d0/(temp_fst+mins1)
           isdml = 1.0d0/sdml

           temp = bxs*(sdl-sdml)*itf
           rolst = maxs1*(rosdl*isdml) - mins1*rol
           vxlst = maxs1*sm - mins1*vxl
           rxlst = rolst*vxlst
           
           vylst = maxs1*(vyl-byl*temp) - mins1*vyl
           rylst = rolst*vylst
           vzlst = maxs1*(vzl-bzl*temp) - mins1*vzl
           rzlst = rolst*vzlst
           
           temp = (rosdl*sdl-bxsq)*itf
           bylst = maxs1*(byl*temp) - mins1*byl
           bzlst = maxs1*(bzl*temp) - mins1*bzl

           vdbstl = vxlst*bxs+vylst*bylst+vzlst*bzlst
           eelst = maxs1*(sdl*eel - ptl*vxl + ptst*sm +      &
     &          bxs*(vxl*bxs+vyl*byl+vzl*bzl-vdbstl))*isdml  &
     &          - mins1*eel		
           
!----- Step 5B. ----------------------------------------------------------|
! compute Ur*
!
           temp_fst = rosdr*sdmr - bxsq
           sign1 = sign(1.0d0,abs(temp_fst)-eps)
           maxs1 = max(0.0d0,sign1)
           mins1 = min(0.0d0,sign1)

           itf = 1.0d0/(temp_fst+mins1)
           isdmr = 1.0d0/sdmr
           
           temp = bxs*(sdr-sdmr)*itf
           rorst = maxs1*(rosdr*isdmr) - mins1*ror
           vxrst = maxs1*sm - mins1*vxr
           rxrst = rorst*vxrst
           
           vyrst = maxs1*(vyr-byr*temp) - mins1*vyr
           ryrst = rorst*vyrst
           vzrst = maxs1*(vzr-bzr*temp) - mins1*vzr
           rzrst = rorst*vzrst
           
           temp = (rosdr*sdr-bxsq)*itf
           byrst = maxs1*(byr*temp) - mins1*byr
           bzrst = maxs1*(bzr*temp) - mins1*bzr
				
           vdbstr = vxrst*bxs+vyrst*byrst+vzrst*bzrst
           eerst = maxs1*((sdr*eer - ptr*vxr  + ptst*sm)*isdmr +   & 
     &          bxs*(vxr*bxs+vyr*byr+vzr*bzr-vdbstr)*isdmr)        &
     &          - mins1*eer

!----- Step 5C. ----------------------------------------------------------|
! compute Ul** and Ur**
!
           sqrtrol = sqrt(rolst)
           sqrtror = sqrt(rorst)

           abbx = abs(bxs)
           signbx = sign(1.0d0,bxs)           
           sign1 = sign(1.0d0,abbx-eps)

           maxs1 = max(0d0,sign1)
           mins1 = -min(0d0,sign1)
           invsumro = maxs1/(sqrtrol + sqrtror)

           temp = invsumro*(sqrtrol*vylst + sqrtror*vyrst  &
     &          + signbx*(byrst-bylst))
           vyldst = vylst*mins1 + temp
           vyrdst = vyrst*mins1 + temp
           ryldst = rylst*mins1 + rolst * temp
           ryrdst = ryrst*mins1 + rorst * temp

           temp = invsumro*(sqrtrol*vzlst + sqrtror*vzrst  &
     &          + signbx*(bzrst-bzlst))
           vzldst = vzlst*mins1 + temp
           vzrdst = vzrst*mins1 + temp
           rzldst = rzlst*mins1 + rolst * temp
           rzrdst = rzrst*mins1 + rorst * temp

           temp = invsumro*(sqrtrol*byrst + sqrtror*bylst  &
     &          + signbx*sqrtrol*sqrtror*(vyrst-vylst))
           byldst = bylst*mins1 + temp
           byrdst = byrst*mins1 + temp
              
           temp = invsumro*(sqrtrol*bzrst + sqrtror*bzlst  &
     &           + signbx*sqrtrol*sqrtror*(vzrst-vzlst))
           bzldst = bzlst*mins1 + temp
           bzrdst = bzrst*mins1 + temp
              
           temp = sm*bxs + vyldst*byldst + vzldst*bzldst
           eeldst = eelst - sqrtrol*signbx*(vdbstl - temp)*maxs1
           eerdst = eerst + sqrtror*signbx*(vdbstr - temp)*maxs1
              
!----- Step 6. ----------------------------------------------------------|
! compute flux
           slst = (sm - abbx/sqrtrol)*maxs1
           srst = (sm + abbx/sqrtror)*maxs1

           sign1 = sign(1.0d0,sm)
           maxs1 = max(0.0d0,sign1)
           mins1 = -min(0.0d0,sign1)

           msl = min(sl,0.0d0)
           msr = max(sr,0.0d0)
           mslst = min(slst,0.0d0)
           msrst = max(srst,0.0d0)

           temp = mslst-msl
           temp1 = msrst-msr

           nflux(mden) = (frol+(rolst-rol)*msl)*maxs1   &
     &           +(fror+(rorst-ror)*msr)*mins1
           nflux(meto) = (feel+(eelst-eel)*msl+(eeldst-eelst)*mslst)*maxs1  & 
     &           +(feer+(eerst-eer)*msr+(eerdst-eerst)*msrst)*mins1
           nflux(mrvu) = (frxl+(rxlst-rxl)*msl)*maxs1   &
     &           +(frxr+(rxrst-rxr)*msr)*mins1
           nflux(mrvv) = (fryl+(rylst-ryl)*msl+(ryldst-rylst)*mslst)*maxs1  & 
     &           +(fryr+(ryrst-ryr)*msr+(ryrdst-ryrst)*msrst)*mins1
           nflux(mrvw) = (frzl+(rzlst-rzl)*msl+(rzldst-rzlst)*mslst)*maxs1  &
     &           +(frzr+(rzrst-rzr)*msr+(rzrdst-rzrst)*msrst)*mins1
           nflux(mbmu) = 0d0
           nflux(mbmv) = (fbyl+(bylst-byl)*msl+(byldst-bylst)*mslst)*maxs1  &
     &           +(fbyr+(byrst-byr)*msr+(byrdst-byrst)*msrst)*mins1
           nflux(mbmw) = (fbzl+(bzlst-bzl)*msl+(bzldst-bzlst)*mslst)*maxs1  &
     &           +(fbzr+(bzrst-bzr)*msr+(bzrdst-bzrst)*msrst)*mins1

      return
      end subroutine HLLD

      subroutine UpdateConsv
      use commons
      use fluxmod
      implicit none
      integer::i,j,k

      do k=ks,ke
      do j=js,je
      do i=is,ie
         
         d(i,j,k) = d(i,j,k)  &
     & +dt*(                  &
     &  (- nflux1(mden,i+1,j,k) &
     &   + nflux1(mden,i  ,j,k))/(x1a(i+1)-x1a(i)) & 
     & +(- nflux2(mden,i,j+1,k) &
     &   + nflux2(mden,i,j  ,k))/(x2a(j+1)-x2a(j)) &
     &      )

         mv1(i,j,k) = mv1(i,j,k) &
     & +dt*( &
     &  (- nflux1(mrv1,i+1,j,k) &
     &   + nflux1(mrv1,i  ,j,k))/(x1a(i+1)-x1a(i)) & 
     & +(- nflux2(mrv1,i,j+1,k) &
     &   + nflux2(mrv1,i,j  ,k))/(x2a(j+1)-x2a(j)) & 
     &      )
         mv2(i,j,k) = mv2(i,j,k) &
     & +dt*( &
     &  (- nflux1(mrv2,i+1,j,k) &
     &   + nflux1(mrv2,i  ,j,k))/(x1a(i+1)-x1a(i)) & 
     & +(- nflux2(mrv2,i,j+1,k) &
     &   + nflux2(mrv2,i,j  ,k))/(x2a(j+1)-x2a(j)) & 
     &      )

         mv3(i,j,k) = mv3(i,j,k)  &
     & +dt*( &
     &  (- nflux1(mrv3,i+1,j,k) &
     &   + nflux1(mrv3,i  ,j,k))/(x1a(i+1)-x1a(i)) & 
     & +(- nflux2(mrv3,i,j+1,k) &
     &   + nflux2(mrv3,i,j  ,k))/(x2a(j+1)-x2a(j)) & 
     &      )

          et(i,j,k) = et(i,j,k) &
     & +dt*( &
     &  (- nflux1(meto,i+1,j,k) &
     &   + nflux1(meto,i  ,j,k))/(x1a(i+1)-x1a(i)) &
     & +(- nflux2(meto,i,j+1,k) &
     &   + nflux2(meto,i,j  ,k))/(x2a(j+1)-x2a(j)) &
     &      )

          b1(i,j,k) = b1(i,j,k) &
     & +dt*( &
     &  (- nflux1(mbm1,i+1,j,k) &
     &   + nflux1(mbm1,i  ,j,k))/(x1a(i+1)-x1a(i)) & 
     & +(- nflux2(mbm1,i,j+1,k) &
     &   + nflux2(mbm1,i,j  ,k))/(x2a(j+1)-x2a(j)) & 
     &      )

          b2(i,j,k) = b2(i,j,k) &
     & +dt*( &
     &  (- nflux1(mbm2,i+1,j,k) &
     &   + nflux1(mbm2,i  ,j,k))/(x1a(i+1)-x1a(i)) & 
     & +(- nflux2(mbm2,i,j+1,k) &
     &   + nflux2(mbm2,i,j  ,k))/(x2a(j+1)-x2a(j)) & 
     &      )

          b3(i,j,k) = b3(i,j,k) &
     & +dt*( &
     &  (- nflux1(mbm3,i+1,j,k) &
     &   + nflux1(mbm3,i  ,j,k))/(x1a(i+1)-x1a(i)) & 
     & +(- nflux2(mbm3,i,j+1,k) &
     &   + nflux2(mbm3,i,j  ,k))/(x2a(j+1)-x2a(j)) & 
     &      )

          bp(i,j,k) = bp(i,j,k) &
     & +dt*( &
     &  (- nflux1(mbps,i+1,j,k) &
     &   + nflux1(mbps,i  ,j,k))/(x1a(i+1)-x1a(i)) & 
     & +(- nflux2(mbps,i,j+1,k) &
     &   + nflux2(mbps,i,j  ,k))/(x2a(j+1)-x2a(j)) &
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
            css  = svc(ncsp,i,j,k)/ svc(nden,i,j,k)**2
            cts  = css  &! cs^2+c_a^2
     &          + (svc(nbm1,i,j,k)**2+svc(nbm2,i,j,k)**2+svc(nbm3,i,j,k)**2)/svc(nden,i,j,k)
            cms  = sqrt((cts +sqrt(cts**2 &
     &      -4.0d0*css*svc(nbm1,i,j,k)**2/svc(nden,i,j,k)))/2.0d0)
            ch1l = ( abs(svc(nve1,i,j,k)) + cms )
            dh1l =  (x1a(i+1)-x1a(i))

            cms  = sqrt((cts +sqrt(cts**2 &
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


      subroutine Output
      use commons
      implicit none
      integer::i,j,k
      character(20),parameter::dirname="meddata/"
      character(40)::filename
      real(8),save::tout
      data tout / 0.0d0 / 
      integer::nout
      data nout / 1 /
      integer,parameter::unitout=13

      if(time .lt. tout+dtout) return

      write(filename,'(a2,i5.5,a4)')"Sc",nout,".xss"

      filename = trim(dirname)//filename
      open(unitout,file=filename,status='replace',form='formatted') 

      write(unitout,*) "# ",time
      k=ks
      do j=js,je
      do i=is,ie
         write(unitout,'(11(1x,E12.3))') x1b(i),x2b(j)                & !  2 
     &                                 , d(i,j,k), p(i,j,k),bp(i,j,k) & !  5
     &                                 ,v1(i,j,k),v2(i,j,k),v3(i,j,k) & !  8
     &                                 ,b1(i,j,k),b2(i,j,k),b3(i,j,k)   ! 11
      enddo
         write(unitout,*)
      enddo
      close(unitout)

      write(6,*) "output:",nout,time

      nout=nout+1
      tout=time

      return
      end subroutine Output

      subroutine makedirs(outdir)
        implicit none
        character(len=*), intent(in) :: outdir
        character(len=256) command
        write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
        write(*, *) trim(command)
        call system(command)
      end subroutine makedirs
