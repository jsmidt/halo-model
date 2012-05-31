module halo
use massfunc
implicit none


contains

! Get 1-Halo Term
real(dl) function P1h(kg)
    real(dl), intent(in) :: kg
    real(dl) :: rombint,tol2
    CALL qromb(P1hi,dlog(minval(hm%nu_m)),dlog(maxval(hm%nu_m)),P1h,kg)
end function P1h
real(dl) function P1hi(lnnu,k)
    real(dl),intent(in) :: lnnu
    real(dl) :: m,k,inu_fnu,omegam,omegal,omegamz,nu,iukm
    omegam = hm%Params%omegab + hm%Params%omegac
    omegal = hm%Params%omegav
    omegamz = omegam*(1.0+hm%z)**3/(omegam*(1.0+hm%z)**3+omegal)
    nu = exp(lnnu)
    m = interpf(dlog(hm%nu_m),hm%m,lnnu)
  call splint(dlog(hm%nu_m),hm%nu_fnum,hm%sp_nu_fnum,size(hm%nu_m),lnnu,inu_fnu)
call splin2(dlog(hm%k),dlog(hm%m),hm%ukm,hm%sp_ukm,hm%mpts,hm%mpts,dlog(k),dlog(m),iukm)
    P1hi = m/hm%rho_c/omegamz*inu_fnu*iukm**2
end function P1hi

! Get 2-Halo Term
real(dl) function P2h(kg)
    real(dl), intent(in) :: kg
    real(dl) :: rombint
    CALL qromb(P2hi,log(minval(hm%nu_m)),dlog(maxval(hm%nu_m)),P2h,kg)
    P2h = P2h**2*interpf(dlog(hm%k),dble(hm%Pk),log(kg))
end function P2h
real(dl) function P2hi(lnnu,k)
    real(dl), intent(in) :: lnnu,k
    real(dl) :: m,inu_fnu,ibias_1,nu,iukm
    nu = exp(lnnu)
    m = interpf(dlog(hm%nu_m),hm%m,lnnu)
  call splint(dlog(hm%nu_m),hm%nu_fnum,hm%sp_nu_fnum,size(hm%nu_m),lnnu,inu_fnu)
  call splint(dlog(hm%nu_m),hm%bias_1,hm%sp_bias_1,size(hm%nu_m),lnnu,ibias_1)
call splin2(dlog(hm%k),dlog(hm%m),hm%ukm,hm%sp_ukm,hm%mpts,hm%mpts,dlog(k),dlog(m),iukm)
    P2hi = ibias_1*inu_fnu*iukm
end function P2hi


end module halo
