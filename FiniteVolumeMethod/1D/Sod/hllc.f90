subroutine HLLC(Ql,Qr,flx)
implicit none
real(8),intent(in)::Ql(:), Qr(:)
real(8),intent(out) :: flx(:)
integer :: i, n
real(8):: Ul(NVAR), Ur(NVAR), Ulst(NVAR), Urst(NVAR)
real(8):: Fl(NVAR), Fr(NVAR)
real(8):: csl,csr
real(8):: sl, sr, cl, cr, sm, pst

    csl = dsqrt(gam*Ql(IPR)/Ql(IDN))
    csr = dsqrt(gam*Qr(IPR)/Qr(IDN))
    sl = min(Ql(IVX),Qr(IVX)) - max(csl,csr)
    sr = max(Ql(IVX),Qr(IVX)) + max(csl,csr)

    cl = Ql(IDN)*(sl - Ql(IVX))
    cr = Qr(IDN)*(sr - Qr(IVX))

    sm = ( -cl*Ql(IVX) + cr*Qr(IVX) - Qr(IPR) + Ql(IPR) )/(-cl + cr)

    pst = Ql(IPR) + cl*(sm - Ql(IVX))
    if( abs( pst /( Qr(IPR) + cr*(sm - Qr(IVX)) ) - 1.0d0 ) > 1.0d-10 ) then
    print*,"error", pst
    stop
    endif

    ! conserved variables in the left and right states
    Ul(IDN) = Ql(IDN)
    Ur(IDN) = Qr(IDN)
    Ulst(IDN) = cl/(sl - sm)
    Urst(IDN) = cr/(sr - sm)

!    print*,Ul(IDN),Ulst(IDN)

    Ul(IMX) = Ql(IDN)*Ql(IVX)
    Ur(IMX) = Qr(IDN)*Qr(IVX)
    Ulst(IMX) = Ulst(IDN)*sm
    Urst(IMX) = Urst(IDN)*sm

    Ul(IEN) = 0.5d0*Ql(IDN)*Ql(IVX)**2 + Ql(IPR)/(gam - 1.0d0)
    Ur(IEN) = 0.5d0*Qr(IDN)*Qr(IVX)**2 + Qr(IPR)/(gam - 1.0d0)
    Ulst(IEN) = ( (sl - Ql(IVX))*Ul(IEN) - Ql(IPR)*Ql(IVX) + pst*sm )/(sl - sm)
    Urst(IEN) = ( (sr - Qr(IVX))*Ur(IEN) - Qr(IPR)*Qr(IVX) + pst*sm )/(sr - sm)

    ! flux in the left and right states
    Fl(IDN) = Ul(IMX)
    Fr(IDN) = Ur(IMX)

    Fl(IMX) = Ql(IPR) + Ql(IDN)*Ql(IVX)**2 
    Fr(IMX) = Qr(IPR) + Qr(IDN)*Qr(IVX)**2 

    Fl(IEN) = ( gam*Ql(IPR)/(gam - 1.0d0) + 0.5d0*Ql(IDN)*Ql(IVX)**2)*Ql(IVX)
    Fr(IEN) = ( gam*Qr(IPR)/(gam - 1.0d0) + 0.5d0*Qr(IDN)*Qr(IVX)**2)*Qr(IVX)


    if( sl > 0.0d0 ) then 
        do n=1,NVAR
            flx(n) = Fl(n)
        enddo
    else if (sr < 0.0d0 ) then
        do n=1,NVAR
             flx(n) = Fr(n)
        enddo
    else if (sm > 0.0d0 ) then
        do n=1,NVAR
            flx(n)  = Fl(n) + sl*(Ulst(n) - Ul(n))
        enddo
    else 
        do n=1,NVAR
            flx(n)  = Fr(n) + sr*(Urst(n) - Ur(n))
        enddo
    endif

return
end subroutine HLLC
